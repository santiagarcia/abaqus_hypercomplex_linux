      module beam_utils_fd
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

      subroutine calculate_element_residual(props, coords, u_vec, R_res)
          real(8), intent(in) :: props(:) 
          real(8), intent(in) :: coords(3,2)
          real(8), intent(in) :: u_vec(12)
          real(8), intent(out) :: R_res(12)
          
          real(8) :: E_mod, nu, Area, Iy, Iz, J_tor, G_mod, mu_neo
          real(8) :: L_elem
          real(8) :: gauss_pts(2), gauss_wts(2)
          integer :: i, k
          
          real(8) :: N_u(2), dN_u(2)
          real(8) :: N_v(4), dN_v(4), ddN_v(4)
          real(8) :: u_ax, du_ax
          real(8) :: v_lat, dv_lat, ddv_lat
          real(8) :: w_lat, dw_lat, ddw_lat
          real(8) :: phi_tor, dphi_tor
          
          real(8) :: lambda, lambda2, inv_lambda
          real(8) :: stress_neo, moment_y, moment_z, torque
          real(8) :: v_dofs(4), w_dofs(4)
          real(8) :: F_ax, F_mz, F_my, F_t
          real(8) :: w_det
          
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
              
              u_ax    = u_vec(1)*N_u(1)  + u_vec(7)*N_u(2)
              du_ax   = u_vec(1)*dN_u(1) + u_vec(7)*dN_u(2)
              
              phi_tor  = u_vec(4)*N_u(1)  + u_vec(10)*N_u(2)
              dphi_tor = u_vec(4)*dN_u(1) + u_vec(10)*dN_u(2)
              
              v_lat   = v_dofs(1)*N_v(1) + v_dofs(2)*N_v(2) + v_dofs(3)*N_v(3) + v_dofs(4)*N_v(4)
              dv_lat  = v_dofs(1)*dN_v(1)+ v_dofs(2)*dN_v(2)+ v_dofs(3)*dN_v(3)+ v_dofs(4)*dN_v(4)
              ddv_lat = v_dofs(1)*ddN_v(1)+v_dofs(2)*ddN_v(2)+v_dofs(3)*ddN_v(3)+v_dofs(4)*ddN_v(4)

              w_lat   = w_dofs(1)*N_v(1) + w_dofs(2)*N_v(2) + w_dofs(3)*N_v(3) + w_dofs(4)*N_v(4)
              dw_lat  = w_dofs(1)*dN_v(1)+ w_dofs(2)*dN_v(2)+ w_dofs(3)*dN_v(3)+ w_dofs(4)*dN_v(4)
              ddw_lat = w_dofs(1)*ddN_v(1)+w_dofs(2)*ddN_v(2)+w_dofs(3)*ddN_v(3)+w_dofs(4)*ddN_v(4)
              
              ! Strains
              lambda2 = (1.0d0 + du_ax)**2 + dv_lat**2 + dw_lat**2
              lambda  = sqrt(lambda2)
              
              ! Stresses
              inv_lambda = 1.0d0/lambda
              stress_neo = (mu_neo/2.0d0) * (2.0d0*lambda - 2.0d0*(inv_lambda**2)) * Area
              moment_z = (E_mod * Iz) * ddv_lat
              moment_y = (E_mod * Iy) * ddw_lat
              torque   = (G_mod * J_tor) * dphi_tor
              
              w_det = gauss_wts(k) * (L_elem / 2.0d0)
              
              ! Axial Term
              F_ax = (stress_neo * inv_lambda) * w_det
              R_res(1) = R_res(1) + F_ax * (1.0d0 + du_ax) * dN_u(1)
              R_res(7) = R_res(7) + F_ax * (1.0d0 + du_ax) * dN_u(2)
              R_res(2)  = R_res(2)  + F_ax * dv_lat * dN_v(1)
              R_res(6)  = R_res(6)  + F_ax * dv_lat * dN_v(2)
              R_res(8)  = R_res(8)  + F_ax * dv_lat * dN_v(3)
              R_res(12) = R_res(12) + F_ax * dv_lat * dN_v(4)
              R_res(3)  = R_res(3)  + F_ax * dw_lat * dN_v(1)
              R_res(5)  = R_res(5)  + F_ax * dw_lat * (dN_v(2) * -1.0d0)
              R_res(9)  = R_res(9)  + F_ax * dw_lat * dN_v(3)
              R_res(11) = R_res(11) + F_ax * dw_lat * (dN_v(4) * -1.0d0)
              
              ! Bending Z
              F_mz = moment_z * w_det
              R_res(2)  = R_res(2)  + F_mz * ddN_v(1)
              R_res(6)  = R_res(6)  + F_mz * ddN_v(2)
              R_res(8)  = R_res(8)  + F_mz * ddN_v(3)
              R_res(12) = R_res(12) + F_mz * ddN_v(4)
              
              ! Bending Y
              F_my = moment_y * w_det
              R_res(3)  = R_res(3)  + F_my * ddN_v(1)
              R_res(5)  = R_res(5)  + F_my * (ddN_v(2) * -1.0d0)
              R_res(9)  = R_res(9)  + F_my * ddN_v(3)
              R_res(11) = R_res(11) + F_my * (ddN_v(4) * -1.0d0)
              
              ! Torsion
              F_t = torque * w_det
              R_res(4)  = R_res(4)  + F_t * dN_u(1)
              R_res(10) = R_res(10) + F_t * dN_u(2)

          end do
      
      end subroutine calculate_element_residual

      end module beam_utils_fd

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS, &
               PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE, &
               TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG, &
               PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS, &
               NJPROP,PERIOD)
      
      use beam_utils_fd
      implicit none
      
      ! Arguments with Explicit Verification - Standard Style
      INTEGER :: NDOFEL, NRHS, NSVARS, NPROPS, MCRD, NNODE, JTYPE, KSTEP, KINC, JELEM, NDLOAD, NPREDF, MLVARX, MDLOAD, NJPROP
      INTEGER :: LFLAGS(*), JDLTYP(MDLOAD,*), JPROPS(*)
      REAL(8)  :: RHS(MLVARX,NRHS), AMATRX(NDOFEL,NDOFEL)
      REAL(8)  :: SVARS(NSVARS), ENERGY(8)
      REAL(8)  :: PROPS(NPROPS)
      REAL(8)  :: COORDS(MCRD,NNODE)
      REAL(8)  :: U(NDOFEL), V(NDOFEL), A(NDOFEL)
      REAL(8)  :: DU(MLVARX,NRHS)
      REAL(8)  :: TIME(2), DTIME, PARAMS(*), ADLMAG(MDLOAD,*), DDLMAG(MDLOAD,*)
      REAL(8)  :: PREDEF(2,NPREDF,NNODE)
      REAL(8)  :: PNEWDT, PERIOD

      ! Local variables
      integer :: i, j
      
      real(8) :: dx, dy, dz, L_val
      real(8) :: ex(3), ey(3), ez(3) 
      real(8) :: T(12,12)
      
      real(8) :: u_tmp_loc(12), u_pert_loc(12)
      real(8) :: R_loc_vec(12), R_glob_vec(12), R_pert_vec(12)
      real(8) :: K_mat_local_val(12,12)
      real(8) :: K_temp(12,12)
      real(8) :: perturbation, epsilon

      logical, save :: first_call = .true.

      ! 1. Guards & Checks FIRST
      ! Initialization of locals to avoid garbage
      T = 0.0d0
      K_mat_local_val = 0.0d0
      K_temp = 0.0d0
      R_loc_vec = 0.0d0
      R_pert_vec = 0.0d0
      R_glob_vec = 0.0d0

      ! Initialization of outputs
      RHS(:,1) = 0.0d0
      AMATRX(:,:) = 0.0d0
      
      ! Debug Print
      if (first_call) then
         write(6,*) '---------------------------------------------------'
         write(6,*) ' FD UEL ENTERED: JELEM=',JELEM
         write(6,*) ' THIS IS THE NON-HYPERDUAL VERSION'
         write(6,*) '---------------------------------------------------'
         call flush(6)
         first_call = .false.
      endif

      ! 5. Local Coordinates System
      dx = COORDS(1,2) - COORDS(1,1)
      dy = COORDS(2,2) - COORDS(2,1)
      dz = COORDS(3,2) - COORDS(3,1)
      L_val = sqrt(dx*dx + dy*dy + dz*dz)
      
      if (L_val < 1.0d-8) then
         write(6,*) 'ERROR: Zero length element'
         call xit
      endif
      
      ex(1) = dx / L_val
      ex(2) = dy / L_val
      ex(3) = dz / L_val
      
      if (abs(ex(3)) .lt. 0.9d0) then
         ey(1) = -ex(2)
         ey(2) =  ex(1)
         ey(3) =  0.0d0
      else
         ey(1) =  0.0d0
         ey(2) = -ex(3)
         ey(3) =  ex(2)
      end if
      L_val = sqrt(ey(1)**2 + ey(2)**2 + ey(3)**2)
      ey(1) = ey(1) / L_val
      ey(2) = ey(2) / L_val
      ey(3) = ey(3) / L_val
      
      ez(1) = ex(2)*ey(3) - ex(3)*ey(2)
      ez(2) = ex(3)*ey(1) - ex(1)*ey(3)
      ez(3) = ex(1)*ey(2) - ex(2)*ey(1)
      
      do i=1,4 
          do j=1,3
             T((i-1)*3 + j, (i-1)*3 + 1) = ex(j)
             T((i-1)*3 + j, (i-1)*3 + 2) = ey(j)
             T((i-1)*3 + j, (i-1)*3 + 3) = ez(j)
          end do
      end do
      
      ! 6. Residual Computation (Value)
      u_tmp_loc = matmul(transpose(T), U(1:12))
      
      call calculate_element_residual(PROPS(1:NPROPS), COORDS(1:MCRD,1:NNODE), u_tmp_loc, R_loc_vec)
      
      R_glob_vec = matmul(T, R_loc_vec)
      
      do i=1,12
          RHS(i,1) = -R_glob_vec(i)
      end do
      
      ! 7. Stiffness Computation (Finite Difference - Forward)
      ! Numerical Tangent (Lower accuracy O(h), prone to convergence issues)
      epsilon = 1.0d-8
      
      do j=1,12
          u_pert_loc = u_tmp_loc
          u_pert_loc(j) = u_pert_loc(j) + epsilon
          
          call calculate_element_residual(PROPS(1:NPROPS), COORDS(1:MCRD,1:NNODE), u_pert_loc, R_pert_vec)
          
          do i=1,12
              K_mat_local_val(i,j) = (R_pert_vec(i) - R_loc_vec(i)) / epsilon
          end do
      end do
      
      K_temp = matmul(K_mat_local_val, transpose(T))
      AMATRX(1:NDOFEL, 1:NDOFEL) = matmul(T, K_temp)
      
      return
      END SUBROUTINE UEL
