module beam_utils_exact
  implicit none
  
  contains

  subroutine shape_funcs(xi, L, N_u, dN_u, N_v, dN_v, ddN_v)
    real(8), intent(in) :: xi, L
    real(8), intent(out), dimension(2) :: N_u, dN_u
    real(8), intent(out), dimension(4) :: N_v, dN_v, ddN_v
    
    real(8) :: one, two, three, four, half, quarter
    
    one = 1.0d0
    two = 2.0d0
    three = 3.0d0
    four = 4.0d0
    half = 0.5d0
    quarter = 0.25d0
    
    ! Linear (u, theta_x)
    N_u(1) = (one - xi) * half
    N_u(2) = (one + xi) * half
    
    dN_u(1) = -half * (two / L)
    dN_u(2) =  half * (two / L)

    ! Hermite Cubic (v, w)
    N_v(1) = quarter * (two - three*xi + xi**3)
    N_v(2) = quarter * (one - xi - xi**2 + xi**3) * (L*half)
    N_v(3) = quarter * (two + three*xi - xi**3)
    N_v(4) = quarter * (-one - xi + xi**2 + xi**3) * (L*half)

    dN_v(1) = quarter * (-three + three*xi**2) * (two/L)
    dN_v(2) = quarter * (-one - two*xi + three*xi**2) 
    dN_v(3) = quarter * (three - three*xi**2) * (two/L)
    dN_v(4) = quarter * (-one + two*xi + three*xi**2)
    
    ddN_v(1) = (1.5d0*xi) * (two/L)**2
    ddN_v(2) = quarter * (-two + 6.0d0*xi) * (two/L)
    ddN_v(3) = (-1.5d0*xi) * (two/L)**2
    ddN_v(4) = quarter * (two + 6.0d0*xi) * (two/L)
    
  end subroutine shape_funcs

  subroutine calculate_element_exact(props, coords, u_vec, R_res, K_mat)
      real(8), intent(in) :: props(:) 
      real(8), intent(in) :: coords(3,2)
      real(8), intent(in) :: u_vec(12)
      real(8), intent(out) :: R_res(12)
      real(8), intent(out) :: K_mat(12,12)
      
      real(8) :: E_mod, nu, Area, Iy, Iz, J_tor, G_mod, mu_neo
      real(8) :: L_elem
      real(8) :: gauss_pts(2), gauss_wts(2)
      integer :: i, j, k
      
      real(8) :: N_u(2), dN_u(2)
      real(8) :: N_v(4), dN_v(4), ddN_v(4)
      real(8) :: u_ax, du_ax
      real(8) :: v_lat, dv_lat, ddv_lat
      real(8) :: w_lat, dw_lat, ddw_lat
      real(8) :: phi_tor, dphi_tor
      
      real(8) :: lambda, lambda2, inv_lambdainv3
      real(8) :: stress_neo, moment_y, moment_z, torque
      real(8) :: v_dofs(4), w_dofs(4)
      real(8) :: F_ax, F_mz, F_my, F_t
      real(8) :: w_det
      
      ! Derivatives
      real(8) :: H_scalar
      real(8) :: Term1
      real(8) :: vec_A(12)
      
      ! Local index arrays
      integer :: dofs_v_idxs(4), dofs_w_idxs(4)
      real(8) :: s_i, s_j
      
      ! Manually assign index arrays to variables instead of DATA statements in execution
      dofs_v_idxs(1) = 2
      dofs_v_idxs(2) = 6
      dofs_v_idxs(3) = 8
      dofs_v_idxs(4) = 12
      
      dofs_w_idxs(1) = 3
      dofs_w_idxs(2) = 5
      dofs_w_idxs(3) = 9
      dofs_w_idxs(4) = 11

      ! Properties
      E_mod = props(1)
      nu    = props(2)
      Area  = props(3)
      Iy    = props(4)
      Iz    = props(5)
      J_tor = props(6)
      
      G_mod = E_mod / (2.0d0 * (1.0d0 + nu))
      mu_neo = G_mod 
      
      ! Geometry
      L_elem = sqrt( (coords(1,2)-coords(1,1))**2 + &
                     (coords(2,2)-coords(2,1))**2 + &
                     (coords(3,2)-coords(3,1))**2 )
                     
      ! Gauss Points
      gauss_pts(1) = -0.577350269189626d0
      gauss_pts(2) =  0.577350269189626d0
      gauss_wts(1) =  1.0d0
      gauss_wts(2) =  1.0d0
      
      do i=1,12
          R_res(i) = 0.0d0
          do j=1,12
              K_mat(i,j) = 0.0d0
          end do
      end do
      
      v_dofs(1) = u_vec(2)
      v_dofs(2) = u_vec(6)
      v_dofs(3) = u_vec(8)
      v_dofs(4) = u_vec(12)
      
      w_dofs(1) = u_vec(3)
      w_dofs(2) = -u_vec(5)
      w_dofs(3) = u_vec(9)
      w_dofs(4) = -u_vec(11)
      
      do k = 1, 2
          call shape_funcs(gauss_pts(k), L_elem, N_u, dN_u, N_v, dN_v, ddN_v)
          
          ! Interpolations
          du_ax   = u_vec(1)*dN_u(1) + u_vec(7)*dN_u(2)
          dphi_tor = u_vec(4)*dN_u(1) + u_vec(10)*dN_u(2)
          
          dv_lat  = v_dofs(1)*dN_v(1)+ v_dofs(2)*dN_v(2)+ v_dofs(3)*dN_v(3)+ v_dofs(4)*dN_v(4)
          ddv_lat = v_dofs(1)*ddN_v(1)+v_dofs(2)*ddN_v(2)+v_dofs(3)*ddN_v(3)+v_dofs(4)*ddN_v(4)

          dw_lat  = w_dofs(1)*dN_v(1)+ w_dofs(2)*dN_v(2)+ w_dofs(3)*dN_v(3)+ w_dofs(4)*dN_v(4)
          ddw_lat = w_dofs(1)*ddN_v(1)+w_dofs(2)*ddN_v(2)+w_dofs(3)*ddN_v(3)+w_dofs(4)*ddN_v(4)
          
          ! Strains
          lambda2 = (1.0d0 + du_ax)**2 + dv_lat**2 + dw_lat**2
          lambda  = sqrt(lambda2)
          
          w_det = gauss_wts(k) * (L_elem / 2.0d0)
          
          F_ax = mu_neo * Area * (1.0d0 - lambda**(-3)) * w_det
          F_mz = (E_mod * Iz) * ddv_lat * w_det
          F_my = (E_mod * Iy) * ddw_lat * w_det
          F_t  = (G_mod * J_tor) * dphi_tor * w_det
          
          ! --- 1. Residual Accumulation ---
          
          vec_A = 0.0d0
          vec_A(1) = (1.0d0+du_ax)*dN_u(1); vec_A(7) = (1.0d0+du_ax)*dN_u(2)
          vec_A(2) = dv_lat*dN_v(1); vec_A(6) = dv_lat*dN_v(2); vec_A(8) = dv_lat*dN_v(3); vec_A(12) = dv_lat*dN_v(4)
          vec_A(3) = dw_lat*dN_v(1); vec_A(5) = dw_lat*dN_v(2)*(-1.0d0); vec_A(9) = dw_lat*dN_v(3); vec_A(11) = dw_lat*dN_v(4)*(-1.0d0)
          
          do i=1,12
             R_res(i) = R_res(i) + F_ax * vec_A(i)
          end do
          
          ! Linear Bending/Torsion Terms
          R_res(2) = R_res(2) + F_mz * ddN_v(1)
          R_res(6) = R_res(6) + F_mz * ddN_v(2)
          R_res(8) = R_res(8) + F_mz * ddN_v(3)
          R_res(12)= R_res(12)+ F_mz * ddN_v(4)
          
          R_res(3) = R_res(3) + F_my * ddN_v(1)
          R_res(5) = R_res(5) + F_my * (ddN_v(2) * (-1.0d0))
          R_res(9) = R_res(9) + F_my * ddN_v(3)
          R_res(11)= R_res(11)+ F_my * (ddN_v(4) * (-1.0d0))
          
          R_res(4) = R_res(4) + F_t * dN_u(1)
          R_res(10)= R_res(10)+ F_t * dN_u(2)
          
          ! --- 2. Stiffness Calculation (Exact) ---
          
          ! 2a. Material Stiffness (for Bending/Torsion - Linear parts)
          
          ! Bending Z (v DOFs: 2, 6, 8, 12)
          do i=1,4
             do j=1,4
                K_mat(dofs_v_idxs(i), dofs_v_idxs(j)) = K_mat(dofs_v_idxs(i), dofs_v_idxs(j)) + (E_mod*Iz*w_det) * ddN_v(i) * ddN_v(j)
             end do
          end do
          
          ! Bending Y (w DOFs: 3, 5, 9, 11)
          do i=1,4
              do j=1,4
                  s_i = 1.0d0; s_j = 1.0d0
                  if(i.eq.2 .or. i.eq.4) s_i = -1.0d0
                  if(j.eq.2 .or. j.eq.4) s_j = -1.0d0
                  
                  K_mat(dofs_w_idxs(i), dofs_w_idxs(j)) = K_mat(dofs_w_idxs(i), dofs_w_idxs(j)) + (E_mod*Iy*w_det) * (ddN_v(i)*s_i) * (ddN_v(j)*s_j)
              end do
          end do
          
          ! Torsion (phi: 4, 10)
          K_mat(4,4)   = K_mat(4,4)   + (G_mod*J_tor*w_det) * dN_u(1)*dN_u(1)
          K_mat(4,10)  = K_mat(4,10)  + (G_mod*J_tor*w_det) * dN_u(1)*dN_u(2)
          K_mat(10,4)  = K_mat(10,4)  + (G_mod*J_tor*w_det) * dN_u(2)*dN_u(1)
          K_mat(10,10) = K_mat(10,10) + (G_mod*J_tor*w_det) * dN_u(2)*dN_u(2)
          
          ! 2b. Geometric & Axial Stiffness (Nonlinear part)
          
          ! Axial Block (1,7)
          K_mat(1,1) = K_mat(1,1) + F_ax * dN_u(1)*dN_u(1)
          K_mat(1,7) = K_mat(1,7) + F_ax * dN_u(1)*dN_u(2)
          K_mat(7,1) = K_mat(7,1) + F_ax * dN_u(2)*dN_u(1)
          K_mat(7,7) = K_mat(7,7) + F_ax * dN_u(2)*dN_u(2)
          
          ! V-block
          do i=1,4
             do j=1,4
                K_mat(dofs_v_idxs(i), dofs_v_idxs(j)) = K_mat(dofs_v_idxs(i), dofs_v_idxs(j)) + F_ax * dN_v(i) * dN_v(j)
             end do
          end do
          
          ! W-block
          do i=1,4
              do j=1,4
                  s_i = 1.0d0; s_j = 1.0d0
                  if(i.eq.2 .or. i.eq.4) s_i = -1.0d0
                  if(j.eq.2 .or. j.eq.4) s_j = -1.0d0
                  
                  K_mat(dofs_w_idxs(i), dofs_w_idxs(j)) = K_mat(dofs_w_idxs(i), dofs_w_idxs(j)) + F_ax * (dN_v(i)*s_i) * (dN_v(j)*s_j)
              end do
          end do
          
          ! 2c. Material/Geometric Coupling (dF_ax/dq * vec_A)
          
          H_scalar = (mu_neo * Area * w_det) * (3.0d0 * lambda**(-4))
          Term1 = H_scalar / lambda
          
          do i=1,12
             do j=1,12
                K_mat(i,j) = K_mat(i,j) + Term1 * vec_A(i) * vec_A(j)
             end do
          end do

      end do
  end subroutine calculate_element_exact
  
end module beam_utils_exact

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS, &
     & PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME, &
     & KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF, &
     & LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
      
      USE beam_utils_exact
      IMPLICIT NONE
      
      ! Variables passed in
      REAL(8) :: RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*), &
     & SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL), &
     & DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),DTIME, &
     & PARAMS(*),ADLMAG(*),PREDEF(2,NPREDF,NNODE), &
     & DDLMAG(*),MDLOAD(MLVARX,*),PNEWDT,JPROPS(*),PERIOD
      
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP, &
     & KINC,JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS(*),MLVARX,NJPROP
     
      ! Local Variables
      REAL(8) :: R_loc(12), K_loc(12,12)
      REAL(8) :: ex(3), ey(3), ez(3), L_val
      REAL(8) :: T(12,12), T_trans(12,12)
      REAL(8) :: u_loc(12)
      REAL(8) :: R_glob(12), K_glob_temp(12,12), K_glob(12,12)
      INTEGER :: i, j

      ! Initialize
      R_loc = 0.0d0
      K_loc = 0.0d0
      RHS(1:12,1) = 0.0d0
      AMATRX(1:NDOFEL,1:NDOFEL) = 0.0d0
      
      ! 1. Calculate Rotation Matrix
      ex(1) = COORDS(1,2) - COORDS(1,1)
      ex(2) = COORDS(2,2) - COORDS(2,1)
      ex(3) = COORDS(3,2) - COORDS(3,1)
      
      L_val = sqrt(ex(1)**2 + ex(2)**2 + ex(3)**2)
      ex = ex / L_val
      
      ! Arbitrary ey (assume beam not vertical for now, or use standard safe vector)
      if (dabs(ex(1)) < 0.9d0) then
         ey(1) = 1.0d0; ey(2) = 0.0d0; ey(3) = 0.0d0
      else
         ey(1) = 0.0d0; ey(2) = 1.0d0; ey(3) = 0.0d0
      end if
      ! Orthonormalize ey w.r.t ex
      ! ey = ey - (ey.ex)ex
      L_val = ey(1)*ex(1) + ey(2)*ex(2) + ey(3)*ex(3)
      ey(1) = ey(1) - L_val*ex(1)
      ey(2) = ey(2) - L_val*ex(2)
      ey(3) = ey(3) - L_val*ex(3)
      L_val = sqrt(ey(1)**2 + ey(2)**2 + ey(3)**2)
      ey = ey / L_val
      
      ez(1) = ex(2)*ey(3) - ex(3)*ey(2)
      ez(2) = ex(3)*ey(1) - ex(1)*ey(3)
      ez(3) = ex(1)*ey(2) - ex(2)*ey(1)
      
      T = 0.0d0
      do i=1,4 
          do j=1,3
             T((i-1)*3 + j, (i-1)*3 + 1) = ex(j)
             T((i-1)*3 + j, (i-1)*3 + 2) = ey(j)
             T((i-1)*3 + j, (i-1)*3 + 3) = ez(j)
          end do
      end do
      T_trans = transpose(T)
      
      ! 2. Transform Displacements to Local System
      ! u_loc = T^T * U_glob
      u_loc = matmul(T_trans, U(1:12))
      
      ! 3. Call Exact Logic
      call calculate_element_exact(PROPS(1:NPROPS), COORDS(1:MCRD,1:NNODE), u_loc, R_loc, K_loc)
      
      ! 4. Transform Forces back to Global
      ! R_glob = T * R_loc
      R_glob = matmul(T, R_loc)
      
      ! 5. Transform Stiffness back to Global
      ! K_glob = T * K_loc * T^T
      K_glob_temp = matmul(K_loc, T_trans)
      K_glob = matmul(T, K_glob_temp)
      
      ! Assign to Abaqus Outputs
      do i=1,12
          RHS(i,1) = -R_glob(i)
          do j=1,12
              AMATRX(i,j) = K_glob(i,j)
          end do
      end do

      return
      end subroutine uel
