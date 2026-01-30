      MODULE complex_utils
      IMPLICIT NONE
      CONTAINS

      subroutine calculate_element_residual_c(props, coords, u_c, R_res)
          real(8), intent(in) :: props(:) 
          real(8), intent(in) :: coords(3,2)
          complex(8), intent(in) :: u_c(12)
          complex(8), intent(out) :: R_res(12)
          
          real(8) :: E_mod, nu, Area, Iy, Iz, J_tor, G_mod, mu_neo
          real(8) :: L_elem
          real(8) :: gauss_pts(2), gauss_wts(2)
          integer :: i, k
          
          complex(8) :: N_u(2), dN_u(2)
          complex(8) :: N_v(4), dN_v(4), ddN_v(4)
          complex(8) :: u_ax, du_ax
          complex(8) :: v_lat, dv_lat, ddv_lat
          complex(8) :: w_lat, dw_lat, ddw_lat
          complex(8) :: phi_tor, dphi_tor
          
          complex(8) :: lambda, lambda2, inv_lambda
          complex(8) :: stress_neo, moment_y, moment_z, torque
          complex(8) :: v_dofs(4), w_dofs(4)
          complex(8) :: F_ax, F_mz, F_my, F_t
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
          
          R_res = cmplx(0.0d0, 0.0d0, 8)
          
          v_dofs(1) = u_c(2)
          v_dofs(2) = u_c(6)
          v_dofs(3) = u_c(8)
          v_dofs(4) = u_c(12)
          
          w_dofs(1) = u_c(3)
          w_dofs(2) = -u_c(5)
          w_dofs(3) = u_c(9)
          w_dofs(4) = -u_c(11)
          
          do k = 1, 2
              call shape_funcs_c(gauss_pts(k), L_elem, N_u, dN_u, N_v, dN_v, ddN_v)
              
              u_ax    = u_c(1)*N_u(1)  + u_c(7)*N_u(2)
              du_ax   = u_c(1)*dN_u(1) + u_c(7)*dN_u(2)
              
              phi_tor  = u_c(4)*N_u(1)  + u_c(10)*N_u(2)
              dphi_tor = u_c(4)*dN_u(1) + u_c(10)*dN_u(2)
              
              v_lat   = v_dofs(1)*N_v(1) + v_dofs(2)*N_v(2) + v_dofs(3)*N_v(3) + v_dofs(4)*N_v(4)
              dv_lat  = v_dofs(1)*dN_v(1)+ v_dofs(2)*dN_v(2)+ v_dofs(3)*dN_v(3)+ v_dofs(4)*dN_v(4)
              ddv_lat = v_dofs(1)*ddN_v(1)+v_dofs(2)*ddN_v(2)+v_dofs(3)*ddN_v(3)+v_dofs(4)*ddN_v(4)

              w_lat   = w_dofs(1)*N_v(1) + w_dofs(2)*N_v(2) + w_dofs(3)*N_v(3) + w_dofs(4)*N_v(4)
              dw_lat  = w_dofs(1)*dN_v(1)+ w_dofs(2)*dN_v(2)+ w_dofs(3)*dN_v(3)+ w_dofs(4)*dN_v(4)
              ddw_lat = w_dofs(1)*ddN_v(1)+w_dofs(2)*ddN_v(2)+w_dofs(3)*ddN_v(3)+w_dofs(4)*ddN_v(4)
              
              ! Strains
              lambda2 = (cmplx(1.0d0,0.0d0,8) + du_ax)**2 + dv_lat**2 + dw_lat**2
              lambda  = sqrt(lambda2)
              
              ! Stresses
              inv_lambda = cmplx(1.0d0,0.0d0,8)/lambda
              stress_neo = cmplx(mu_neo/2.0d0,0.0d0,8) * (cmplx(2.0d0,0.0d0,8)*lambda - cmplx(2.0d0,0.0d0,8)*(inv_lambda**2)) * cmplx(Area,0.0d0,8)
              moment_z = cmplx(E_mod * Iz, 0.0d0, 8) * ddv_lat
              moment_y = cmplx(E_mod * Iy, 0.0d0, 8) * ddw_lat
              torque   = cmplx(G_mod * J_tor, 0.0d0, 8) * dphi_tor
              
              w_det = gauss_wts(k) * (L_elem / 2.0d0)
              
              ! Axial Term
              F_ax = (stress_neo * inv_lambda) * cmplx(w_det,0.0d0,8)
              
              R_res(1) = R_res(1) + F_ax * (cmplx(1.0d0,0.0d0,8) + du_ax) * dN_u(1)
              R_res(7) = R_res(7) + F_ax * (cmplx(1.0d0,0.0d0,8) + du_ax) * dN_u(2)
              R_res(2)  = R_res(2)  + F_ax * dv_lat * dN_v(1)
              R_res(6)  = R_res(6)  + F_ax * dv_lat * dN_v(2)
              R_res(8)  = R_res(8)  + F_ax * dv_lat * dN_v(3)
              R_res(12) = R_res(12) + F_ax * dv_lat * dN_v(4)
              R_res(3)  = R_res(3)  + F_ax * dw_lat * dN_v(1)
              R_res(5)  = R_res(5)  + F_ax * dw_lat * (dN_v(2) * cmplx(-1.0d0,0.0d0,8))
              R_res(9)  = R_res(9)  + F_ax * dw_lat * dN_v(3)
              R_res(11) = R_res(11) + F_ax * dw_lat * (dN_v(4) * cmplx(-1.0d0,0.0d0,8))
              
              ! Bending Z
              F_mz = moment_z * cmplx(w_det,0.0d0,8)
              R_res(2)  = R_res(2)  + F_mz * ddN_v(1)
              R_res(6)  = R_res(6)  + F_mz * ddN_v(2)
              R_res(8)  = R_res(8)  + F_mz * ddN_v(3)
              R_res(12) = R_res(12) + F_mz * ddN_v(4)
              
              ! Bending Y
              F_my = moment_y * cmplx(w_det,0.0d0,8)
              R_res(3)  = R_res(3)  + F_my * ddN_v(1)
              R_res(5)  = R_res(5)  + F_my * (ddN_v(2) * cmplx(-1.0d0,0.0d0,8))
              R_res(9)  = R_res(9)  + F_my * ddN_v(3)
              R_res(11) = R_res(11) + F_my * (ddN_v(4) * cmplx(-1.0d0,0.0d0,8))
              
              ! Torsion
              F_t = torque * cmplx(w_det,0.0d0,8)
              R_res(4)  = R_res(4)  + F_t * dN_u(1)
              R_res(10) = R_res(10) + F_t * dN_u(2)

          end do
      end subroutine calculate_element_residual_c

      subroutine shape_funcs_c(xi, L, N_u, dN_u, N_v, dN_v, ddN_v)
        real(8), intent(in) :: xi, L
        complex(8), intent(out) :: N_u(2), dN_u(2)
        complex(8), intent(out) :: N_v(4), dN_v(4), ddN_v(4)
        
        complex(8) :: xi_c
        complex(8) :: one, two, three, four, half, quarter, L_c
        
        xi_c = cmplx(xi, 0.0d0, 8)
        L_c  = cmplx(L, 0.0d0, 8)
        one = cmplx(1.0d0, 0.0d0, 8)
        two = cmplx(2.0d0, 0.0d0, 8)
        three = cmplx(3.0d0, 0.0d0, 8)
        four = cmplx(4.0d0, 0.0d0, 8)
        half = cmplx(0.5d0, 0.0d0, 8)
        quarter = cmplx(0.25d0, 0.0d0, 8)
        
        ! Linear (u, theta_x)
        N_u(1) = (one - xi_c) * half
        N_u(2) = (one + xi_c) * half
        
        dN_u(1) = -half * (two / L_c)
        dN_u(2) =  half * (two / L_c)

        ! Hermite Cubic (v, w)
        N_v(1) = quarter * (two - three*xi_c + xi_c**3)
        N_v(2) = quarter * (one - xi_c - xi_c**2 + xi_c**3) * (L_c*half)
        N_v(3) = quarter * (two + three*xi_c - xi_c**3)
        N_v(4) = quarter * (-one - xi_c + xi_c**2 + xi_c**3) * (L_c*half)

        dN_v(1) = quarter * (-three + three*xi_c**2) * (two/L_c)
        dN_v(2) = quarter * (-one - two*xi_c + three*xi_c**2) 
        dN_v(3) = quarter * (three - three*xi_c**2) * (two/L_c)
        dN_v(4) = quarter * (-one + two*xi_c + three*xi_c**2)
        
        ddN_v(1) = (cmplx(1.5d0,0.0d0,8)*xi_c) * (two/L_c)**2
        ddN_v(2) = quarter * (-two + cmplx(6.0d0,0.0d0,8)*xi_c) * (two/L_c)
        ddN_v(3) = (cmplx(-1.5d0,0.0d0,8)*xi_c) * (two/L_c)**2
        ddN_v(4) = quarter * (two + cmplx(6.0d0,0.0d0,8)*xi_c) * (two/L_c)
        
      end subroutine shape_funcs_c

      END MODULE complex_utils

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS, &
               PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE, &
               TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG, &
               PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS, &
               NJPROP,PERIOD)
      
      USE complex_utils
      ! Hypercomplex (Complex Step) UEL Implementation
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
      complex(8) :: u_c(12), R_c(12) 
      integer :: i, j
      
      real(8) :: dx, dy, dz, L_val
      real(8) :: ex(3), ey(3), ez(3) 
      real(8) :: T(12,12)
      
      real(8) :: u_tmp_loc(12)
      real(8) :: R_loc_vec(12), R_glob_vec(12)
      real(8) :: K_mat_local_val(12,12)
      real(8) :: K_temp(12,12)
      real(8), parameter :: h_step = 1.0d-20
      
      logical, save :: first_call = .true.

      ! Initialization
      T = 0.0d0
      K_mat_local_val = 0.0d0
      K_temp = 0.0d0
      R_loc_vec = 0.0d0
      R_glob_vec = 0.0d0
      RHS(:,1) = 0.0d0
      AMATRX(:,:) = 0.0d0
      
      if (first_call) then
         write(6,*) '---------------------------------------------------'
         write(6,*) ' HYPERCOMPLEX UEL STARTED'
         write(6,*) ' Method: Complex Step Differentiation'
         write(6,*) ' Step Size:', h_step
         write(6,*) '---------------------------------------------------'
         call flush(6)
         first_call = .false.
      endif

      ! Coordinate System (Same as other UELs)
      dx = COORDS(1,2) - COORDS(1,1)
      dy = COORDS(2,2) - COORDS(2,1)
      dz = COORDS(3,2) - COORDS(3,1)
      L_val = sqrt(dx*dx + dy*dy + dz*dz)
      
      if (L_val < 1.0d-8) then
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
      ey = ey / L_val
      
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
      
      u_tmp_loc = matmul(transpose(T), U(1:12))
      
      ! 1. Residual (Value)
      do i=1,12
         u_c(i) = cmplx(u_tmp_loc(i), 0.0d0, 8)
      end do
      
      call calculate_element_residual_c(PROPS(1:NPROPS), COORDS(1:MCRD,1:NNODE), u_c, R_c)
      
      do i=1,12
          R_loc_vec(i) = dreal(R_c(i))
      end do
      R_glob_vec = matmul(T, R_loc_vec)
      
      do i=1,12
          RHS(i,1) = -R_glob_vec(i)
      end do
      
      ! 2. Stiffness (Complex Step)
      do j=1,12
          u_c = cmplx(u_tmp_loc, 0.0d0, 8) ! Reset
          u_c(j) = u_c(j) + cmplx(0.0d0, h_step, 8) ! Perturb
          
          call calculate_element_residual_c(PROPS(1:NPROPS), COORDS(1:MCRD,1:NNODE), u_c, R_c)
          
          do i=1,12
              K_mat_local_val(i,j) = aimag(R_c(i)) / h_step
          end do
      end do
      
      K_temp = matmul(K_mat_local_val, transpose(T))
      AMATRX(1:NDOFEL, 1:NDOFEL) = matmul(T, K_temp)
      
      return
      END SUBROUTINE UEL
