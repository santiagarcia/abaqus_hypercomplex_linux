      module hyperdual_mod
      implicit none
      private

      ! Type definition for Dual Vector numbers (2 Channels)
      ! Optimized for calculating two gradients simultaneously
      ! Dropped second-derivative (e12) cross term for performance
      type, public :: hd
          real(8) :: val
          real(8) :: e1
          real(8) :: e2
      end type hd

      ! Public interfaces
      public :: assignment(=)
      public :: operator(+), operator(-), operator(*), operator(/)
      public :: operator(**)
      public :: sqrt, abs, exp, log, sin, cos

      ! Constructors
      interface hd
        module procedure real_to_hd
      end interface

      ! Assignment
      interface assignment(=)
          module procedure assign_real_to_hd
      end interface

      ! Operators
      interface operator(+)
          module procedure add_hd_hd, add_hd_r, add_r_hd
      end interface

      interface operator(-)
          module procedure sub_hd_hd, sub_hd_r, sub_r_hd, neg_hd
      end interface

      interface operator(*)
          module procedure mul_hd_hd, mul_hd_r, mul_r_hd
      end interface

      interface operator(/)
          module procedure div_hd_hd, div_hd_r, div_r_hd
      end interface

      interface operator(**)
          module procedure pow_hd_i
      end interface

      ! Functions
      interface sqrt
          module procedure sqrt_hd
      end interface
      
      interface abs
          module procedure abs_hd
      end interface

      interface exp
          module procedure exp_hd
      end interface

      interface log
          module procedure log_hd
      end interface

      interface sin
          module procedure sin_hd
      end interface

      interface cos
           module procedure cos_hd
      end interface

      contains

      ! Constructor
      function real_to_hd(x) result(res)
          real(8), intent(in) :: x
          type(hd) :: res
          res%val = x
          res%e1  = 0.0d0
          res%e2  = 0.0d0
      end function real_to_hd

      ! Assignment
      subroutine assign_real_to_hd(res, x)
          type(hd), intent(out) :: res
          real(8), intent(in) :: x
          res%val = x
          res%e1  = 0.0d0
          res%e2  = 0.0d0
      end subroutine assign_real_to_hd

      ! Addition
      function add_hd_hd(a, b) result(res)
          type(hd), intent(in) :: a, b
          type(hd) :: res
          res%val = a%val + b%val
          res%e1  = a%e1  + b%e1
          res%e2  = a%e2  + b%e2
      end function add_hd_hd

      function add_hd_r(a, b) result(res)
          type(hd), intent(in) :: a
          real(8), intent(in) :: b
          type(hd) :: res
          res%val = a%val + b
          res%e1  = a%e1
          res%e2  = a%e2
      end function add_hd_r

      function add_r_hd(a, b) result(res)
          real(8), intent(in) :: a
          type(hd), intent(in) :: b
          type(hd) :: res
          res%val = a + b%val
          res%e1  = b%e1
          res%e2  = b%e2
      end function add_r_hd

      ! Subtraction
      function sub_hd_hd(a, b) result(res)
          type(hd), intent(in) :: a, b
          type(hd) :: res
          res%val = a%val - b%val
          res%e1  = a%e1  - b%e1
          res%e2  = a%e2  - b%e2
      end function sub_hd_hd

      function sub_hd_r(a, b) result(res)
          type(hd), intent(in) :: a
          real(8), intent(in) :: b
          type(hd) :: res
          res%val = a%val - b
          res%e1  = a%e1
          res%e2  = a%e2
      end function sub_hd_r

      function sub_r_hd(a, b) result(res)
          real(8), intent(in) :: a
          type(hd), intent(in) :: b
          type(hd) :: res
          res%val = a - b%val
          res%e1  = -b%e1
          res%e2  = -b%e2
      end function sub_r_hd

      function neg_hd(a) result(res)
          type(hd), intent(in) :: a
          type(hd) :: res
          res%val = -a%val
          res%e1  = -a%e1
          res%e2  = -a%e2
      end function neg_hd

      ! Multiplication
      function mul_hd_hd(a, b) result(res)
          type(hd), intent(in) :: a, b
          type(hd) :: res
          res%val = a%val * b%val
          res%e1  = a%val * b%e1 + a%e1 * b%val
          res%e2  = a%val * b%e2 + a%e2 * b%val
      end function mul_hd_hd

      function mul_hd_r(a, b) result(res)
          type(hd), intent(in) :: a
          real(8), intent(in) :: b
          type(hd) :: res
          res%val = a%val * b
          res%e1  = a%e1  * b
          res%e2  = a%e2  * b
      end function mul_hd_r

      function mul_r_hd(a, b) result(res)
          real(8), intent(in) :: a
          type(hd), intent(in) :: b
          type(hd) :: res
          res%val = a * b%val
          res%e1  = a * b%e1
          res%e2  = a * b%e2
      end function mul_r_hd

      ! Division
      function div_hd_hd(a, b) result(res)
          type(hd), intent(in) :: a, b
          type(hd) :: res
          real(8) :: inv_b
          inv_b = 1.0d0 / b%val
          res%val = a%val * inv_b
          res%e1  = (a%e1 - res%val * b%e1) * inv_b
          res%e2  = (a%e2 - res%val * b%e2) * inv_b
      end function div_hd_hd

      function div_hd_r(a, b) result(res)
          type(hd), intent(in) :: a
          real(8), intent(in) :: b
          type(hd) :: res
          real(8) :: inv_b
          inv_b = 1.0d0 / b
          res%val = a%val * inv_b
          res%e1  = a%e1  * inv_b
          res%e2  = a%e2  * inv_b
      end function div_hd_r

      function div_r_hd(a, b) result(res)
          real(8), intent(in) :: a
          type(hd), intent(in) :: b
          type(hd) :: res
          real(8) :: inv_b, inv_b2
          inv_b = 1.0d0 / b%val
          inv_b2 = inv_b * inv_b
          
          res%val = a * inv_b
          res%e1  = - a * b%e1 * inv_b2
          res%e2  = - a * b%e2 * inv_b2
      end function div_r_hd

      ! Power integer
      function pow_hd_i(a, n) result(res)
          type(hd), intent(in) :: a
          integer, intent(in) :: n
          type(hd) :: res
          integer :: i
          
          if (n == 0) then
              res%val = 1.0d0
              res%e1 = 0.0d0
              res%e2 = 0.0d0
              return
          end if
          if (n == 1) then
              res = a
              return
          end if
          if (n == 2) then
              res = a * a
              return
          end if
          
          res = a
          do i = 2, n
              res = res * a
          end do
      end function pow_hd_i

      ! Sqrt
      function sqrt_hd(a) result(res)
          type(hd), intent(in) :: a
          type(hd) :: res
          real(8) :: s
          s = sqrt(a%val)
          res%val = s
          res%e1  = 0.5d0 * a%e1 / s
          res%e2  = 0.5d0 * a%e2 / s
      end function sqrt_hd

      ! Abs
      function abs_hd(a) result(res)
          type(hd), intent(in) :: a
          type(hd) :: res
          real(8) :: sgn
          if (a%val >= 0.0d0) then
              sgn = 1.0d0
          else
              sgn = -1.0d0
          endif
          res%val = sgn * a%val
          res%e1  = sgn * a%e1
          res%e2  = sgn * a%e2
      end function abs_hd

      ! Exp
      function exp_hd(a) result(res)
          type(hd), intent(in) :: a
          type(hd) :: res
          real(8) :: ev
          ev = exp(a%val)
          res%val = ev
          res%e1  = ev * a%e1
          res%e2  = ev * a%e2
      end function exp_hd

      ! Log
      function log_hd(a) result(res)
          type(hd), intent(in) :: a
          type(hd) :: res
          real(8) :: inv_a
          inv_a = 1.0d0 / a%val
          res%val = log(a%val)
          res%e1  = a%e1 * inv_a
          res%e2  = a%e2 * inv_a
      end function log_hd

      ! Sin
      function sin_hd(a) result(res)
            type(hd), intent(in) :: a
            type(hd) :: res
            real(8) :: s, c
            s = sin(a%val)
            c = cos(a%val)
            res%val = s
            res%e1  = c * a%e1
            res%e2  = c * a%e2
      end function sin_hd

        ! Cos
        function cos_hd(a) result(res)
              type(hd), intent(in) :: a
              type(hd) :: res
              real(8) :: s, c
              s = sin(a%val)
              c = cos(a%val)
              res%val = c
              res%e1  = -s * a%e1
              res%e2  = -s * a%e2
        end function cos_hd

      end module hyperdual_mod
      
      module beam_utils
      use hyperdual_mod
      implicit none
      
      contains

      subroutine shape_funcs(xi, L, N_u, dN_u, N_v, dN_v, ddN_v)
        real(8), intent(in) :: xi, L
        type(hd), intent(out), dimension(2) :: N_u, dN_u
        type(hd), intent(out), dimension(4) :: N_v, dN_v, ddN_v
        
        type(hd) :: xi_hd
        type(hd) :: one, two, three, four, half, quarter
        
        xi_hd = hd(xi)
        one = hd(1.0d0)
        two = hd(2.0d0)
        three = hd(3.0d0)
        four = hd(4.0d0)
        half = hd(0.5d0)
        quarter = hd(0.25d0)
        
        ! Linear (u, theta_x)
        N_u(1) = (one - xi_hd) * half
        N_u(2) = (one + xi_hd) * half
        
        dN_u(1) = -half * (two / hd(L))
        dN_u(2) =  half * (two / hd(L))

        ! Hermite Cubic (v, w)
        N_v(1) = quarter * (two - three*xi_hd + xi_hd**3)
        N_v(2) = quarter * (one - xi_hd - xi_hd**2 + xi_hd**3) * (hd(L)*half)
        N_v(3) = quarter * (two + three*xi_hd - xi_hd**3)
        N_v(4) = quarter * (-one - xi_hd + xi_hd**2 + xi_hd**3) * (hd(L)*half)

        dN_v(1) = quarter * (-three + three*xi_hd**2) * (two/hd(L))
        dN_v(2) = quarter * (-one - two*xi_hd + three*xi_hd**2) 
        dN_v(3) = quarter * (three - three*xi_hd**2) * (two/hd(L))
        dN_v(4) = quarter * (-one + two*xi_hd + three*xi_hd**2)
        
        ddN_v(1) = (hd(1.5d0)*xi_hd) * (two/hd(L))**2
        ddN_v(2) = quarter * (-two + hd(6.0d0)*xi_hd) * (two/hd(L))
        ddN_v(3) = (hd(-1.5d0)*xi_hd) * (two/hd(L))**2
        ddN_v(4) = quarter * (two + hd(6.0d0)*xi_hd) * (two/hd(L))
        
      end subroutine shape_funcs

      subroutine calculate_element_residual(props, coords, u_hd, R_res)
          real(8), intent(in) :: props(:) 
          real(8), intent(in) :: coords(3,2)
          type(hd), intent(in) :: u_hd(12)
          type(hd), intent(out) :: R_res(12)
          
          real(8) :: E_mod, nu, Area, Iy, Iz, J_tor, G_mod, mu_neo
          real(8) :: L_elem
          real(8) :: gauss_pts(2), gauss_wts(2)
          integer :: i, k
          
          type(hd) :: N_u(2), dN_u(2)
          type(hd) :: N_v(4), dN_v(4), ddN_v(4)
          type(hd) :: u_ax, du_ax
          type(hd) :: v_lat, dv_lat, ddv_lat
          type(hd) :: w_lat, dw_lat, ddw_lat
          type(hd) :: phi_tor, dphi_tor
          
          type(hd) :: lambda, lambda2, inv_lambda
          type(hd) :: stress_neo, moment_y, moment_z, torque
          type(hd) :: v_dofs(4), w_dofs(4)
          type(hd) :: F_ax, F_mz, F_my, F_t
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
              R_res(i) = hd(0.0d0)
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
              lambda2 = (hd(1.0d0) + du_ax)**2 + dv_lat**2 + dw_lat**2
              lambda  = sqrt(lambda2)
              
              ! Stresses
              inv_lambda = hd(1.0d0)/lambda
              stress_neo = hd(mu_neo/2.0d0) * (hd(2.0d0)*lambda - hd(2.0d0)*(inv_lambda**2)) * hd(Area)
              moment_z = hd(E_mod * Iz) * ddv_lat
              moment_y = hd(E_mod * Iy) * ddw_lat
              torque   = hd(G_mod * J_tor) * dphi_tor
              
              w_det = gauss_wts(k) * (L_elem / 2.0d0)
              
              ! Axial Term
              F_ax = (stress_neo * inv_lambda) * hd(w_det)
              R_res(1) = R_res(1) + F_ax * (hd(1.0d0) + du_ax) * dN_u(1)
              R_res(7) = R_res(7) + F_ax * (hd(1.0d0) + du_ax) * dN_u(2)
              R_res(2)  = R_res(2)  + F_ax * dv_lat * dN_v(1)
              R_res(6)  = R_res(6)  + F_ax * dv_lat * dN_v(2)
              R_res(8)  = R_res(8)  + F_ax * dv_lat * dN_v(3)
              R_res(12) = R_res(12) + F_ax * dv_lat * dN_v(4)
              R_res(3)  = R_res(3)  + F_ax * dw_lat * dN_v(1)
              R_res(5)  = R_res(5)  + F_ax * dw_lat * (dN_v(2) * hd(-1.0d0))
              R_res(9)  = R_res(9)  + F_ax * dw_lat * dN_v(3)
              R_res(11) = R_res(11) + F_ax * dw_lat * (dN_v(4) * hd(-1.0d0))
              
              ! Bending Z
              F_mz = moment_z * hd(w_det)
              R_res(2)  = R_res(2)  + F_mz * ddN_v(1)
              R_res(6)  = R_res(6)  + F_mz * ddN_v(2)
              R_res(8)  = R_res(8)  + F_mz * ddN_v(3)
              R_res(12) = R_res(12) + F_mz * ddN_v(4)
              
              ! Bending Y
              F_my = moment_y * hd(w_det)
              R_res(3)  = R_res(3)  + F_my * ddN_v(1)
              R_res(5)  = R_res(5)  + F_my * (ddN_v(2) * hd(-1.0d0))
              R_res(9)  = R_res(9)  + F_my * ddN_v(3)
              R_res(11) = R_res(11) + F_my * (ddN_v(4) * hd(-1.0d0))
              
              ! Torsion
              F_t = torque * hd(w_det)
              R_res(4)  = R_res(4)  + F_t * dN_u(1)
              R_res(10) = R_res(10) + F_t * dN_u(2)

          end do
      
      end subroutine calculate_element_residual

      end module beam_utils

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS, &
               PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE, &
               TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG, &
               PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS, &
               NJPROP,PERIOD)
      
      use hyperdual_mod
      use beam_utils
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
      type(hd) :: u_hd(12), R_hd(12)
      integer :: i, j
      
      real(8) :: dx, dy, dz, L_val
      real(8) :: ex(3), ey(3), ez(3) 
      real(8) :: T(12,12)
      
      real(8) :: u_tmp_loc(12)
      real(8) :: R_loc_vec(12), R_glob_vec(12)
      real(8) :: K_mat_local_val(12,12)
      real(8) :: K_temp(12,12)
      
      logical, save :: first_call = .true.
      logical, parameter :: SMOKE_TEST = .false. 

      ! 1. Guards & Checks FIRST
      ! Initialization of locals to avoid garbage
      T = 0.0d0
      K_mat_local_val = 0.0d0
      K_temp = 0.0d0
      R_loc_vec = 0.0d0
      R_glob_vec = 0.0d0

      ! Initialization of outputs
      RHS(:,1) = 0.0d0
      AMATRX(:,:) = 0.0d0
      
      ! Debug Print
      if (first_call) then
         write(6,*) '---------------------------------------------------'
         write(6,*) ' UEL ENTERED: JELEM=',JELEM,' NDOFEL=',NDOFEL,' NNODE=',NNODE,' NRHS=',NRHS,' NSVARS=',NSVARS
         write(6,*) ' Check: MCRD  =', MCRD
         write(6,*) ' Check: NPROPS=', NPROPS
         write(6,*) ' Check: LFLAGS=', LFLAGS(1:6)
         write(6,*) ' SMOKE_TEST =', SMOKE_TEST
         write(6,*) '---------------------------------------------------'
         call flush(6)
         first_call = .false.
      endif

      if (NDOFEL .ne. 12) then
         write(6,*) 'CRITICAL ERROR: NDOFEL must be 12. Got', NDOFEL
         call xit
      endif
      if (NNODE .ne. 2) then
         write(6,*) 'CRITICAL ERROR: NCOORDS/NNODE must be 2. Got', NNODE
         call xit
      endif
      if (MCRD .ne. 3) then
         write(6,*) 'CRITICAL ERROR: MCRD must be 3. Got', MCRD
         call xit
      endif
      if (NPROPS .lt. 6) then
         write(6,*) 'CRITICAL ERROR: NPROPS too small. Need 6. Got', NPROPS
         call xit
      endif
      
      ! 3. Smoke Test Mode
      if (SMOKE_TEST) then
          do i=1,12
             AMATRX(i,i) = 1.0d8 ! Simple stiffness
             RHS(i,1) = -1.0d4 * U(i) ! Restore force
          end do
          return
      endif
      
      ! 4. Operational Flag Check
      ! LFLAGS(1): 1=Static, 2=Dynamic. 
      ! LFLAGS(3): 1=Normal, 100=Perturbation? 
      ! Standard execution: just proceed unless specific skips needed.
      
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
      
      do i=1,12
         u_hd(i) = hd(u_tmp_loc(i)) ! e1=0
      end do
      
      call calculate_element_residual(PROPS(1:NPROPS), COORDS(1:MCRD,1:NNODE), u_hd, R_hd)
      
      do i=1,12
          R_loc_vec(i) = R_hd(i)%val
      end do
      R_glob_vec = matmul(T, R_loc_vec)
      
      do i=1,12
          RHS(i,1) = -R_glob_vec(i)
      end do
      
      ! 7. Stiffness Computation (Tangent)
      ! OPTIMIZED: Vectorized using Hyperdual channels (e1, e2)
      ! Computes 2 columns of Stiffness Matrix per pass.
      ! Moved u_tmp_loc calculation outside loops.
      
      u_tmp_loc = matmul(transpose(T), U(1:12))

      do j=1,12,2
          ! Initialize Hyperduals with base state (and zero derivatives)
          do i=1,12
              u_hd(i)    = hd(u_tmp_loc(i))
          end do
          
          ! Perturb j-th DOF on e1 channel
          u_hd(j)%e1 = 1.0d0
          
          ! Perturb (j+1)-th DOF on e2 channel (safe for NDOF=12)
          u_hd(j+1)%e2 = 1.0d0
          
          call calculate_element_residual(PROPS(1:NPROPS), COORDS(1:MCRD,1:NNODE), u_hd, R_hd)
          
          do i=1,12
              ! Column j from e1 partials
              K_mat_local_val(i,j)   = R_hd(i)%e1
              ! Column j+1 from e2 partials
              K_mat_local_val(i,j+1) = R_hd(i)%e2
          end do
      end do
      
      K_temp = matmul(K_mat_local_val, transpose(T))
      AMATRX(1:NDOFEL, 1:NDOFEL) = matmul(T, K_temp)
      
      return
      END SUBROUTINE UEL
