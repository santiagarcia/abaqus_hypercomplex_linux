import pandas as pd
import numpy as np
import plotly.graph_objects as go
import re
import os

# Configuration
DAT_FILE = os.path.join(os.path.dirname(__file__), '../sim/cantilever_uel.dat')
OUTPUT_HTML = os.path.join(os.path.dirname(__file__), '../sim/beam_animation.html')
NUM_POINTS = 50  # Points along the beam for smooth curve

def parse_dat_file(filepath):
    """Parses Abaqus .dat file for Node Output."""
    if not os.path.exists(filepath):
        print(f"Error: {filepath} not found.")
        return []

    with open(filepath, 'r') as f:
        lines = f.readlines()

    increments = []
    current_increment = {}
    
    # Regex to capture Increment Time
    # "INCREMENT     X SUMMARY"
    re_inc_summary = re.compile(r'INCREMENT\s+(\d+)\s+SUMMARY')
    re_time_comp = re.compile(r'TOTAL TIME COMPLETED\s+([\d\.E\+\-]+)')
    
    # Regex for Node Output table lines: "Node U1 U2 U3 UR1 UR2 UR3"
    # Example:
    # 1     -7.3459036E-36  2.5112420E-37  0.0000000E+00  0.0000000E+00  0.0000000E+00  5.8924313E-36 
    re_node_line = re.compile(r'^\s*(\d+)\s+([\d\.E\+\-]+)\s+([\d\.E\+\-]+)\s+([\d\.E\+\-]+)\s+([\d\.E\+\-]+)\s+([\d\.E\+\-]+)\s+([\d\.E\+\-]+)')

    reading_nodes = False
    
    inc_data = None 
    
    for line in lines:
        if "INCREMENT" in line and "SUMMARY" in line:
            if inc_data and 'nodes' in inc_data and len(inc_data['nodes']) > 0:
                 increments.append(inc_data)
            
            m = re_inc_summary.search(line)
            if m:
                inc_num = int(m.group(1))
                inc_data = {'id': inc_num, 'nodes': {}, 'time': 0.0}
        
        if inc_data is not None:
             if "TOTAL TIME COMPLETED" in line:
                 m = re_time_comp.search(line)
                 if m:
                     try:
                         inc_data['time'] = float(m.group(1))
                     except:
                         pass
            
             if "THE FOLLOWING TABLE IS PRINTED FOR ALL NODES" in line:
                 reading_nodes = True
                 continue
             
             if reading_nodes:
                 m = re_node_line.match(line)
                 if m:
                     node_id = int(m.group(1))
                     vals = [float(x) for x in m.groups()[1:]]
                     # Store [U1, U2, U3, UR1, UR2, UR3]
                     inc_data['nodes'][node_id] = vals
                 elif "MAXIMUM" in line:
                     reading_nodes = False
    
    # Append last
    if inc_data and 'nodes' in inc_data and len(inc_data['nodes']) > 0:
         increments.append(inc_data)
         
    return increments

def hermite_shape(xi, L):
    """
    Returns shape functions for v (transverse)
    """
    one = 1.0
    two = 2.0
    three = 3.0
    quart = 0.25
    half_L = L / 2.0
    
    Nv = np.zeros(4)
    Nv[0] = quart * (two - three*xi + xi**3)
    Nv[1] = quart * (one - xi - xi**2 + xi**3) * half_L
    Nv[2] = quart * (two + three*xi - xi**3)
    Nv[3] = quart * (-one - xi + xi**2 + xi**3) * half_L
    
    return Nv

def linear_shape(xi):
    # N1 = (1 - xi)/2
    # N2 = (1 + xi)/2
    Nu = np.zeros(2)
    Nu[0] = 0.5 * (1.0 - xi)
    Nu[1] = 0.5 * (1.0 + xi)
    return Nu

def interpolate_beam(inc_data, L=10.0, num_points=50):
    """
    Returns x, y, z arrays for the deformed beam at this increment.
    """
    # Extract DOFs. Expect Nodes 1 and 2.
    if 1 not in inc_data['nodes'] or 2 not in inc_data['nodes']:
        # Fallback if nodes missing
        return [], [], []

    n1 = inc_data['nodes'][1]
    n2 = inc_data['nodes'][2]
    
    # Map to element DOFs
    u1_dofs = [n1[0], n2[0]]
    
    # v_lat (U2): N1_u2, N1_ur3, N2_u2, N2_ur3
    v_dofs = [n1[1], n1[5], n2[1], n2[5]]
    
    # w_lat (U3): N1_u3, -N1_ur2, N2_u3, -N2_ur2
    w_dofs = [n1[2], -n1[4], n2[2], -n2[4]]
    
    xi_list = np.linspace(-1, 1, num_points)
    
    x_coords = []
    y_coords = []
    z_coords = []
    
    for xi in xi_list:
        # Map xi [-1,1] to x [0, L]
        x0 = (xi + 1.0) * L / 2.0
        
        Nu = linear_shape(xi)
        Nv = hermite_shape(xi, L)
        
        # Displacements
        du = Nu[0]*u1_dofs[0] + Nu[1]*u1_dofs[1]
        dv = np.dot(Nv, v_dofs)
        dw = np.dot(Nv, w_dofs)
        
        x_coords.append(x0 + du)
        y_coords.append(0.0 + dv)
        z_coords.append(0.0 + dw)
        
    return x_coords, y_coords, z_coords

def main():
    print(f"Reading {DAT_FILE}...")
    increments = parse_dat_file(DAT_FILE)
    print(f"Found {len(increments)} increments.")
    
    if not increments:
        print("No increments found. Check dat file format.")
        return

    # Create AnimationFrames
    frames = []
    
    # Determine bounds for plotting
    all_x, all_y, all_z = [], [], []
    for inc in increments:
        x, y, z = interpolate_beam(inc)
        all_x.extend(x)
        all_y.extend(y)
        all_z.extend(z)
        
    x_min, x_max = min(all_x), max(all_x)
    y_min, y_max = min(all_y), max(all_y)
    z_min, z_max = min(all_z), max(all_z) 
    
    # Pad limits
    pad_x = (x_max - x_min)*0.1 if (x_max!=x_min) else 1.0
    pad_y = (y_max - y_min)*0.1 if (y_max!=y_min) else 1.0
    
    x_range = [x_min-pad_x, x_max+pad_x]
    y_range = [y_min-pad_y, y_max+pad_y]
    z_range = [-1, 1] # Fixed small range for 2D problem
    
    
    for i, inc in enumerate(increments):
        x, y, z = interpolate_beam(inc)
        
        frame = go.Frame(
            data=[go.Scatter3d(
                x=x, y=y, z=z,
                mode='lines+markers',
                line=dict(color='blue', width=8),
                marker=dict(size=4),
                name='Beam'
            )],
            name=f"Inc{i}"
        )
        frames.append(frame)
        
    # Initial Data
    x0, y0, z0 = interpolate_beam(increments[0])
    
    fig = go.Figure(
        data=[go.Scatter3d(
            x=x0, y=y0, z=z0,
            mode='lines+markers',
            line=dict(color='blue', width=8),
            marker=dict(size=4),
            name='Beam'
        )],
        layout=go.Layout(
            title="Cantilever Beam Deflection (UEL)",
            scene=dict(
                xaxis=dict(range=x_range, title='X (Axial)'),
                yaxis=dict(range=y_range, title='Y (Deflection)'),
                zaxis=dict(range=z_range, title='Z'), 
                aspectmode='manual',
                aspectratio=dict(x=4, y=1, z=0.5)
            ),
            updatemenus=[{
                "buttons": [
                    {
                        "args": [None, {"frame": {"duration": 100, "redraw": True},
                                        "fromcurrent": True, "transition": {"duration": 0}}],
                        "label": "Play",
                        "method": "animate"
                    },
                    {
                        "args": [[None], {"frame": {"duration": 0, "redraw": True},
                                          "mode": "immediate",
                                          "transition": {"duration": 0}}],
                        "label": "Pause",
                        "method": "animate"
                    }
                ],
                "direction": "left",
                "pad": {"r": 10, "t": 87},
                "showactive": False,
                "type": "buttons",
                "x": 0.1,
                "xanchor": "right",
                "y": 0,
                "yanchor": "top"
            }],
            sliders= [{
                "active": 0,
                "yanchor": "top",
                "xanchor": "left",
                "currentvalue": {
                    "font": {"size": 20},
                    "prefix": "Increment:",
                    "visible": True,
                    "xanchor": "right"
                },
                "transition": {"duration": 300, "easing": "cubic-in-out"},
                "pad": {"b": 10, "t": 50},
                "len": 0.9,
                "x": 0.1,
                "y": 0,
                "steps": [
                    {
                        "args": [
                            [f"Inc{k}"],
                            {"frame": {"duration": 300, "redraw": True},
                             "mode": "immediate",
                             "transition": {"duration": 300}}
                        ],
                        "label": f"{inc['id']}",
                        "method": "animate"
                    }
                    for k, inc in enumerate(increments)
                ]
            }]
        ),
        frames=frames
    )

    print(f"Writing animation to {OUTPUT_HTML}")
    fig.write_html(OUTPUT_HTML)
    print("Done.")

if __name__ == "__main__":
    main()