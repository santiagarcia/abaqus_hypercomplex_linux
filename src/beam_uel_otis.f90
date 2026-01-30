
      include 'master_parameters.f90'
      include 'real_utils.f90'
      include 'otim6n1.f90'

      module beam_utils_otis
      use OTIM6N1
      implicit none
      
      contains

      subroutine shape_funcs(xi, L, N_u, dN_u, N_v, dN_v, ddN_v)
        real(8), intent(in) :: xi, L
        type(ONUMM6N1), intent(out), dimension(2) :: N_u, dN_u
        type(ONUMM6N1), intent(out), dimension(4) :: N_v, dN_v, ddN_v
        
        type(ONUMM6N1) :: xi_hd
        type(ONUMM6N1) :: one, two, three, four, half, quarter
        
        xi_hd = xi
        one = 1.0d0
        two = 2.0d0
        three = 3.0d0
        four = 4.0d0
        half = 0.5d0
        quarter = 0.25d0
        
        ! Linear (u, theta_x)
        N_u(1) = (one - xi_hd) * half
        N_u(2) = (one + xi_hd) * half
        
        dN_u(1) = -half * (two / L)
        dN_u(2) =  half * (two / L)

        ! Hermite Cubic (v, w)
        N_v(1) = quarter * (two - three*xi_hd + xi_hd**3.0d0)
        N_v(2) = quarter * (one - xi_hd - xi_hd**2.0d0 + xi_hd**3.0d0) * (L*half)
        N_v(3) = quarter * (two + three*xi_hd - xi_hd**3.0d0)
        N_v(4) = quarter * (-one - xi_hd + xi_hd**2.0d0 + xi_hd**3.0d0) * (L*half)

        dN_v(1) = quarter * (-three + three*xi_hd**2.0d0) * (two/L)
        dN_v(2) = quarter * (-one - two*xi_hd + three*xi_hd**2.0d0) 
        dN_v(3) = quarter * (three - three*xi_hd**2.0d0) * (two/L)
        dN_v(4) = quarter * (-one + two*xi_hd + three*xi_hd**2.0d0)
        
        ddN_v(1) = (1.5d0*xi_hd) * (two/L)**2.0d0
        ddN_v(2) = quarter * (-two + 6.0d0*xi_hd) * (two/L)
        ddN_v(3) = (-1.5d0*xi_hd) * (two/L)**2.0d0
        ddN_v(4) = quarter * (two + 6.0d0*xi_hd) * (two/L)
        
      end subroutine shape_funcs

      subroutine calculate_element_residual(props, coords, u_hd, R_res)
          real(8), intent(in) :: props(:) 
          real(8), intent(in) :: coords(3,2)
          type(ONUMM6N1), intent(in) :: u_hd(12)
          type(ONUMM6N1), intent(out) :: R_res(12)
          
          real(8) :: E_mod, nu, Area, Iy, Iz, J_tor, G_mod, mu_neo
          real(8) :: L_elem
          real(8) :: gauss_pts(2), gauss_wts(2)
          integer :: i, k
          
          type(ONUMM6N1) :: N_u(2), dN_u(2)
          type(ONUMM6N1) :: N_v(4), dN_v(4), ddN_v(4)
          type(ONUMM6N1) :: u_ax, du_ax
          type(ONUMM6N1) :: v_lat, dv_lat, ddv_lat
          type(ONUMM6N1) :: w_lat, dw_lat, ddw_lat
          type(ONUMM6N1) :: phi_tor, dphi_tor
          
          type(ONUMM6N1) :: lambda, lambda2, inv_lambda
          type(ONUMM6N1) :: stress_neo, moment_y, moment_z, torque
          type(ONUMM6N1) :: v_dofs(4), w_dofs(4)
          type(ONUMM6N1) :: F_ax, F_mz, F_my, F_t
          real(8) :: w_det
          type(ONUMM6N1) :: half, one, two, inv_b2
          
          one = 1.0d0
          two = 2.0d0
          half = 0.5d0

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
          
          v_dofs(1) = u_hd(2)
          v_dofs(2) = u_hd(6)
          v_dofs(3) = u_hd(8)
          v_dofs(4) = u_hd(12)
          
          w_dofs(1) = u_hd(3)
          w_dofs(2) = -u_hd(5)
          w_dofs(3) = u_hd(9)
          w_dofs(4) = -u_hd(11)
          
          do k = 1, 2
              call shape_funcs(gauss_pts(k), L_elem, N_u, dN_u, N_v, dN_v, ddN_v)
              
              u_ax    = u_hd(1)*N_u(1)  + u_hd(7)*N_u(2)
              du_ax   = u_hd(1)*dN_u(1) + u_hd(7)*dN_u(2)
              
              phi_tor  = u_hd(4)*N_u(1)  + u_hd(10)*N_u(2)
              dphi_tor = u_hd(4)*dN_u(1) + u_hd(10)*dN_u(2)
              
              v_lat   = v_dofs(1)*N_v(1) + v_dofs(2)*N_v(2) + v_dofs(3)*N_v(3) + v_dofs(4)*N_v(4)
              dv_lat  = v_dofs(1)*dN_v(1)+ v_dofs(2)*dN_v(2)+ v_dofs(3)*dN_v(3)+ v_dofs(4)*dN_v(4)
              ddv_lat = v_dofs(1)*ddN_v(1)+v_dofs(2)*ddN_v(2)+v_dofs(3)*ddN_v(3)+v_dofs(4)*ddN_v(4)

              w_lat   = w_dofs(1)*N_v(1) + w_dofs(2)*N_v(2) + w_dofs(3)*N_v(3) + w_dofs(4)*N_v(4)
              dw_lat  = w_dofs(1)*dN_v(1)+ w_dofs(2)*dN_v(2)+ w_dofs(3)*dN_v(3)+ w_dofs(4)*dN_v(4)
              ddw_lat = w_dofs(1)*ddN_v(1)+w_dofs(2)*ddN_v(2)+w_dofs(3)*ddN_v(3)+w_dofs(4)*ddN_v(4)
              
              ! Strains
              lambda2 = (one + du_ax)**2.0d0 + dv_lat**2.0d0 + dw_lat**2.0d0
              lambda  = sqrt(lambda2)
              
              ! Stresses
              inv_lambda = one/lambda
              stress_neo = (mu_neo*half) * (two*lambda - two*(inv_lambda**2.0d0)) * Area
              moment_z = (E_mod * Iz) * ddv_lat
              moment_y = (E_mod * Iy) * ddw_lat
              torque   = (G_mod * J_tor) * dphi_tor
              
              w_det = gauss_wts(k) * (L_elem / 2.0d0)
              
              ! Axial Term
              F_ax = (stress_neo * inv_lambda) * w_det
              R_res(1) = R_res(1) + F_ax * (one + du_ax) * dN_u(1)
              R_res(7) = R_res(7) + F_ax * (one + du_ax) * dN_u(2)
              R_res(2)  = R_res(2)  + F_ax * dv_lat * dN_v(1)
              R_res(6)  = R_res(6)  + F_ax * dv_lat * dN_v(2)
              R_res(8)  = R_res(8)  + F_ax * dv_lat * dN_v(3)
              R_res(12) = R_res(12) + F_ax * dv_lat * dN_v(4)
              R_res(3)  = R_res(3)  + F_ax * dw_lat * dN_v(1)
              R_res(5)  = R_res(5)  + F_ax * dw_lat * (dN_v(2) * (-1.0d0))
              R_res(9)  = R_res(9)  + F_ax * dw_lat * dN_v(3)
              R_res(11) = R_res(11) + F_ax * dw_lat * (dN_v(4) * (-1.0d0))
              
              ! Bending Z
              F_mz = moment_z * w_det
              R_res(2)  = R_res(2)  + F_mz * ddN_v(1)
              R_res(6)  = R_res(6)  + F_mz * ddN_v(2)
              R_res(8)  = R_res(8)  + F_mz * ddN_v(3)
              R_res(12) = R_res(12) + F_mz * ddN_v(4)
              
              ! Bending Y
              F_my = moment_y * w_det
              R_res(3)  = R_res(3)  + F_my * ddN_v(1)
              R_res(5)  = R_res(5)  + F_my * (ddN_v(2) * (-1.0d0))
              R_res(9)  = R_res(9)  + F_my * ddN_v(3)
              R_res(11) = R_res(11) + F_my * (ddN_v(4) * (-1.0d0))
              
              ! Torsion
              F_t = torque * w_det
              R_res(4)  = R_res(4)  + F_t * dN_u(1)
              R_res(10) = R_res(10) + F_t * dN_u(2)

          end do
      
      end subroutine calculate_element_residual

      end module beam_utils_otis
      
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS, &
               PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE, &
               TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG, &
               PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS, &
               NJPROP,PERIOD)
      
      use OTIM6N1
      use beam_utils_otis
      implicit none
      
      ! Arguments
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
      type(ONUMM6N1) :: u_hd(12), R_hd(12)
      integer :: i, j
      
      real(8) :: dx, dy, dz, L_val
      real(8) :: ex(3), ey(3), ez(3) 
      real(8) :: T(12,12)
      
      real(8) :: u_tmp_loc(12)
      real(8) :: R_loc_vec(12), R_glob_vec(12)
      real(8) :: K_mat_local_val(12,12)
      real(8) :: K_temp(12,12)

      ! Initialize
      T = 0.0d0
      K_mat_local_val = 0.0d0
      K_temp = 0.0d0
      R_loc_vec = 0.0d0
      R_glob_vec = 0.0d0

      RHS(:,1) = 0.0d0
      AMATRX(:,:) = 0.0d0
      
      ! Coordinate rotation logic
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
      
      ! 6. Residual Computation
      u_tmp_loc = matmul(transpose(T), U(1:12))
      
      do i=1,12
         u_hd(i) = u_tmp_loc(i) 
      end do
      
      call calculate_element_residual(PROPS(1:NPROPS), COORDS(1:MCRD,1:NNODE), u_hd, R_hd)
      
      do i=1,12
          R_loc_vec(i) = R_hd(i)%R
      end do
      R_glob_vec = matmul(T, R_loc_vec)
      
      do i=1,12
          RHS(i,1) = -R_glob_vec(i)
      end do
      
      ! 7. Stiffness Computation (Tangent) with Otis (6 Directions)
      ! Pass 1: Cols 1-6
      ! Pass 2: Cols 7-12
      
      do j=1,12,6
          ! Initialize
          do i=1,12
              u_hd(i) = u_tmp_loc(i)
          end do
          
          ! Perturb 6 directions
                 u_hd(j)%E1 = 1.0d0
          if(j+1<=12) u_hd(j+1)%E2 = 1.0d0
          if(j+2<=12) u_hd(j+2)%E3 = 1.0d0
          if(j+3<=12) u_hd(j+3)%E4 = 1.0d0
          if(j+4<=12) u_hd(j+4)%E5 = 1.0d0
          if(j+5<=12) u_hd(j+5)%E6 = 1.0d0
          
          call calculate_element_residual(PROPS(1:NPROPS), COORDS(1:MCRD,1:NNODE), u_hd, R_hd)
          
          do i=1,12
                     K_mat_local_val(i,j)   = R_hd(i)%E1
             if(j+1<=12) K_mat_local_val(i,j+1) = R_hd(i)%E2
             if(j+2<=12) K_mat_local_val(i,j+2) = R_hd(i)%E3
             if(j+3<=12) K_mat_local_val(i,j+3) = R_hd(i)%E4
             if(j+4<=12) K_mat_local_val(i,j+4) = R_hd(i)%E5
             if(j+5<=12) K_mat_local_val(i,j+5) = R_hd(i)%E6
          end do
      end do
      
      K_temp = matmul(K_mat_local_val, transpose(T))
      AMATRX(1:NDOFEL, 1:NDOFEL) = matmul(T, K_temp)
      
      return
      END SUBROUTINE UEL
