
import sympy as sp
from sympy.utilities.codegen import codegen

def generate_fortran_code():
    # Define Symbols
    u = sp.MatrixSymbol('u', 12, 1) # Local DOFs
    
    # Properties
    E_mod, nu, Area, Iy, Iz, J_tor = sp.symbols('E_mod nu Area Iy Iz J_tor')
    L_elem = sp.symbols('L_elem')
    w_det = sp.symbols('w_det') # Weight * Jacobian
    
    G_mod = E_mod / (2 * (1 + nu))
    mu_neo = G_mod
    
    # Shape function values at a gauss point (passed as inputs to the subroutine)
    # We will assume the subroutine loops over gauss points and accumulates K
    N_u1, N_u2 = sp.symbols('N_u(1) N_u(2)')
    dN_u1, dN_u2 = sp.symbols('dN_u(1) dN_u(2)')
    
    N_v1, N_v2, N_v3, N_v4 = sp.symbols('N_v(1) N_v(2) N_v(3) N_v(4)')
    dN_v1, dN_v2, dN_v3, dN_v4 = sp.symbols('dN_v(1) dN_v(2) dN_v(3) dN_v(4)')
    ddN_v1, ddN_v2, ddN_v3, ddN_v4 = sp.symbols('ddN_v(1) ddN_v(2) ddN_v(3) ddN_v(4)')
    
    # Map DOFs to Field Variables
    # u_ax
    du_ax = u[0]*dN_u1 + u[6]*dN_u2
    
    # torsion
    dphi_tor = u[3]*dN_u1 + u[9]*dN_u2
    
    # v_lat (using u2, u6, u8, u12)
    dv_lat = u[1]*dN_v1 + u[5]*dN_v2 + u[7]*dN_v3 + u[11]*dN_v4
    ddv_lat = u[1]*ddN_v1 + u[5]*ddN_v2 + u[7]*ddN_v3 + u[11]*ddN_v4
    
    # w_lat (using u3, -u5, u9, -u11) 
    # NOTE: The FD code uses:
    # w_dofs(1) = u_vec(3)
    # w_dofs(2) = -u_vec(5)
    # w_dofs(3) = u_vec(9)
    # w_dofs(4) = -u_vec(11)
    
    # Be careful with indices (0-based in sympy, 1-based in fortran description)
    # u_vec(1) -> u[0]
    # ...
    # u_vec(5) is Rotation about Y (DOF 5). 
    # Standard beam: theta_y = -dw/dx. 
    # Abaqus UEL Manual or typical expectation: DOF 5 is RotY.
    
    dw_lat = u[2]*dN_v1 - u[4]*dN_v2 + u[8]*dN_v3 - u[10]*dN_v4
    ddw_lat = u[2]*ddN_v1 - u[4]*ddN_v2 + u[8]*ddN_v3 - u[10]*ddN_v4
    
    # Strains
    # lambda2 = (1.0d0 + du_ax)**2 + dv_lat**2 + dw_lat**2
    lambda2 = (1 + du_ax)**2 + dv_lat**2 + dw_lat**2
    lam = sp.sqrt(lambda2)
    inv_lambda = 1/lam
    
    # Stresses (Internal Forces intensities)
    # stress_neo = (mu_neo/2.0d0) * (2.0d0*lambda - 2.0d0*(inv_lambda**2)) * Area
    stress_neo = (mu_neo/2) * (2*lam - 2*(inv_lambda**2)) * Area
    
    moment_z = (E_mod * Iz) * ddv_lat
    moment_y = (E_mod * Iy) * ddw_lat
    torque   = (G_mod * J_tor) * dphi_tor
    
    # Multiplied by determinant (integration weight)
    # F_ax = (stress_neo * inv_lambda) * w_det
    F_mt = stress_neo * inv_lambda * w_det
    
    # Residual Vector Construction
    R = sp.zeros(12, 1)
    
    # Axial / Geometric Stiffness Parts
    # R_res(1) = R_res(1) + F_ax * (1.0d0 + du_ax) * dN_u(1)
    term_geo = F_mt * (1 + du_ax)
    R[0] += term_geo * dN_u1
    R[6] += term_geo * dN_u2
    
    # R_res(2)  = R_res(2)  + F_ax * dv_lat * dN_v(1)
    # ...
    term_geo_v = F_mt * dv_lat
    R[1] += term_geo_v * dN_v1
    R[5] += term_geo_v * dN_v2
    R[7] += term_geo_v * dN_v3
    R[11]+= term_geo_v * dN_v4
    
    # R_res(3)  = R_res(3)  + F_ax * dw_lat * dN_v(1)
    term_geo_w = F_mt * dw_lat
    R[2] += term_geo_w * dN_v1
    R[4] -= term_geo_w * dN_v2 # Note Negative sign from w_dofs definition above
    R[8] += term_geo_w * dN_v3
    R[10]-= term_geo_w * dN_v4
    
    # Bending Z
    F_mz = moment_z * w_det
    R[1] += F_mz * ddN_v1
    R[5] += F_mz * ddN_v2
    R[7] += F_mz * ddN_v3
    R[11]+= F_mz * ddN_v4
    
    # Bending Y
    F_my = moment_y * w_det
    R[2] += F_my * ddN_v1
    R[4] -= F_my * ddN_v2
    R[8] += F_my * ddN_v3
    R[10]-= F_my * ddN_v4
    
    # Torsion
    F_t = torque * w_det
    R[3] += F_t * dN_u1
    R[9] += F_t * dN_u2
    
    # Compute Jacobian (Stiffness Matrix)
    K = R.jacobian(u)
    
    # Generate Code
    print("Generating Fortran code...")
    from sympy.printing.fortran import fcode
    
    # Use CSE to optimize efficiently
    replacements, reduced_K = sp.cse(K)
    
    with open('stiffness_code.inc', 'w') as f:
        f.write("      ! Generated Analytical Stiffness\n")
        f.write("      ! Common Subexpressions\n")
        for sym, expr in replacements:
            f.write(f"      real(8) :: {sym}\n")
            f.write(f"      {sym} = {fcode(expr, source_format='free')}\n")
            
        f.write("\n      ! Stiffness Matrix Accumulation\n")
        # reduced_K is [Matrix(...)]
        matrix_res = reduced_K[0]
        rows, cols = matrix_res.shape # 12, 12
        for i in range(rows):
            for j in range(cols):
                val = matrix_res[i, j]
                if val != 0:
                    code = fcode(val, source_format='free')
                    # Remember K_mat_local_val is accumulated over Gauss points
                    # i+1, j+1 for 1-based indexing
                    f.write(f"      K_mat_local_val({i+1},{j+1}) = K_mat_local_val({i+1},{j+1}) + {code}\n")

if __name__ == "__main__":
    generate_fortran_code()
