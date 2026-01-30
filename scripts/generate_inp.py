import sys
import numpy as np

def generate_inp(filename, num_elements, force_y, num_increments=10, length=10.0):
    num_nodes = num_elements + 1
    dx = length / num_elements
    
    # Calculate time increment to achieve target number of increments
    # Total time is 1.0
    dt = 1.0 / float(num_increments)
    
    with open(filename, 'w') as f:
        f.write("*HEADING\n")
        f.write(f"Cantilever Beam {num_elements} Elements Force {force_y}\n")
        
        # Nodes
        f.write("*NODE\n")
        for i in range(num_nodes):
            x = i * dx
            # ID, x, y, z
            f.write(f"{i+1}, {x:.4f}, 0.0, 0.0\n")
            
        # User Element Definition
        # Note: VARIABLES=12 (Number of State variables, separate from DOFs)
        f.write("*USER ELEMENT, NODES=2, TYPE=U1, PROPERTIES=6, COORDINATES=3, VARIABLES=12\n")
        # DOF list: 1-6. Abaqus applies this to all nodes if only one line is given.
        f.write("1, 2, 3, 4, 5, 6\n")
        
        # Elements
        f.write("*ELEMENT, TYPE=U1, ELSET=BEAM\n")
        for i in range(num_elements):
            # ElemID, Node1, Node2
            f.write(f"{i+1}, {i+1}, {i+2}\n")
            
        # Properties
        f.write("*UEL PROPERTY, ELSET=BEAM\n")
        # E, nu, Area, Iy, Iz, J
        f.write("210.e9, 0.3, 0.01, 1.e-5, 1.e-5, 2.e-5\n")
        
        # Boundary Conditions (Fix Node 1 completely)
        f.write("*BOUNDARY\n")
        f.write("1, 1, 6, 0.0\n")
        
        # Step
        f.write("*STEP, NLGEOM=YES, INC=10000\n")
        f.write("*STATIC\n")
        # Initial, Total, Min, Max
        f.write(f"{dt}, 1.0, 1e-10, {dt}\n")
        
        # Load at Tip (Last Node)
        f.write("*CLOAD\n")
        # NodeID, DOF, Mag
        f.write(f"{num_nodes}, 2, {force_y}\n")
        
        # Output
        f.write("*NODE PRINT\n")
        f.write("U,\n")
        f.write("*OUTPUT, FIELD\n")
        f.write("*NODE OUTPUT\n")
        f.write("U\n")
        f.write("*END STEP\n")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python3 generate_inp.py <filename> <num_elements> <force_y> [num_increments]")
        sys.exit(1)
        
    fname = sys.argv[1]
    nel = int(sys.argv[2])
    f_val = float(sys.argv[3])
    
    n_inc = 10
    if len(sys.argv) > 4:
        n_inc = int(sys.argv[4])
    
    generate_inp(fname, nel, f_val, n_inc)
