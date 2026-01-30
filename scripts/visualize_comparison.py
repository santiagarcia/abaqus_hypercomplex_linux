import plotly.graph_objects as go
from plotly.subplots import make_subplots
import re
import numpy as np
import sys
import os

def interpolate_beam(inc_data):
    """
    Given increment data with nodal coordinates and displacements,
    reconstruct the full beam shape.
    """
    nodes = inc_data['nodes']
    sorted_node_ids = sorted(nodes.keys())
    
    x_full, y_full, z_full = [], [], []
    
    # We need to handle the case where we don't have explicit coords, 
    # but we know it's a cantilever along X.
    # In parse_all_simulation_data we reconstruct 'coord' based on node count.
    
    for i in range(len(sorted_node_ids) - 1):
        n1 = nodes[sorted_node_ids[i]]
        n2 = nodes[sorted_node_ids[i+1]]
        
        # Original Coords
        X1, Y1, Z1 = n1['coord']
        X2, Y2, Z2 = n2['coord']
        
        # Displacements
        u1, v1, w1 = n1['disp'][0], n1['disp'][1], n1['disp'][2]
        u2, v2, w2 = n2['disp'][0], n2['disp'][1], n2['disp'][2]
        
        # Deformed Coords
        def_x1 = X1 + u1
        def_y1 = Y1 + v1
        def_z1 = Z1 + w1
        
        def_x2 = X2 + u2
        def_y2 = Y2 + v2
        def_z2 = Z2 + w2
        
        x_full.append(def_x1)
        y_full.append(def_y1)
        z_full.append(def_z1)
        
    # Add last point
    if sorted_node_ids:
        last_n = nodes[sorted_node_ids[-1]]
        x_full.append(last_n['coord'][0] + last_n['disp'][0])
        y_full.append(last_n['coord'][1] + last_n['disp'][1])
        z_full.append(last_n['coord'][2] + last_n['disp'][2])
        
    return x_full, y_full, z_full

def parse_dat_file(dat_file):
    """
    Parses .dat file for Nodal Displacements.
    Returns list of dicts: [{'id': inc_num, 'nodes': {node_id: {disp:[], coord:[]}}}]
    """
    if not os.path.exists(dat_file):
        print(f"Warning: {dat_file} not found")
        return []
        
    with open(dat_file, 'r') as f:
        lines = f.readlines()
        
    all_increments = []
    
    current_inc_id = None
    current_nodes = {}
    
    # State flags
    in_table = False
    
    for line in lines:
        stripped = line.strip()
        
        # 1. Detect Increment Header
        # "INCREMENT     5 SUMMARY"
        if "INCREMENT" in line and "SUMMARY" in line:
            # If we had a previous increment active, save it.
            # BUT: In Abaqus .dat, the "INCREMENT x SUMMARY" usually comes *before* 
            # the "NODE OUTPUT" table for that increment (based on user grep).
            # Wait, user grep showed:
            # INCREMENT 5 SUMMARY ... NODE OUTPUT ... TABLE
            # So detecting "INCREMENT x SUMMARY" is the start of a block.
            
            if current_inc_id is not None and current_nodes:
                all_increments.append({'id': current_inc_id, 'nodes': current_nodes.copy()})
                
            m = re.search(r'INCREMENT\s+(\d+)\s+SUMMARY', line)
            if m:
                current_inc_id = int(m.group(1))
                current_nodes = {}
                in_table = False
                
        # 2. Detect Table Start
        if "THE FOLLOWING TABLE IS PRINTED" in line:
            in_table = True
            continue
            
        # 3. Parse Table Data
        if in_table:
            # Check for end of table (empty line or text like "MAXIMUM")
            if not stripped:
                continue # Skip empty lines inside table? Or is empty line the terminator?
                # Abaqus often puts empty lines. Let's look for non-digit start.
            
            tokens = stripped.split()
            if not tokens: 
                continue
                
            # Check if first token is an integer
            try:
                node_id = int(tokens[0])
            except ValueError:
                # Not an integer (e.g. "NODE", "MAXIMUM", "MINIMUM")
                # This breaks the table reading
                if "NODE" in tokens[0] or "FOOT-" in tokens[0] or "NOTE" in tokens[0]:
                    continue # Header lines
                else:
                    in_table = False # End of table logic (e.g. MAXIMUM)
                    continue

            # If we are here, we have a node_id.
            # Expecting 6 floats after it: U1, U2, U3, UR1, UR2, UR3
            if len(tokens) >= 7:
                try:
                    disps = [float(x) for x in tokens[1:7]]
                    current_nodes[node_id] = {'disp': disps, 'coord': [0.0,0.0,0.0]}
                except ValueError:
                    pass
    
    # End of file: save last
    if current_inc_id is not None and current_nodes:
        all_increments.append({'id': current_inc_id, 'nodes': current_nodes.copy()})
        
    # Post-processing: Add Coords
    # We assume a cantilever of Length=10 along X
    if all_increments:
        # Determine dx based on number of nodes in the last valid increment
        last_inc = all_increments[-1]
        node_ids = sorted(last_inc['nodes'].keys())
        if not node_ids: return []
        
        max_id = node_ids[-1]
        L = 10.0
        if max_id > 1:
            dx = L / (max_id - 1)
        else:
            dx = 10.0
            
        for inc in all_increments:
            for nid, data in inc['nodes'].items():
                x_coord = (nid - 1) * dx
                data['coord'] = [x_coord, 0.0, 0.0]
                
    return all_increments

def parse_msg_file(msg_file):
    """
    Parses .msg file for iterations.
    """
    iteration_history = []
    if not os.path.exists(msg_file):
        return []
        
    with open(msg_file, 'r') as f:
        lines = f.readlines()
        
    current_inc = 0
    current_iters = 0
    
    for line in lines:
        if "INCREMENT" in line and "STARTS" in line:
            m = re.search(r'INCREMENT\s+(\d+)\s+STARTS', line)
            if m:
                current_inc = int(m.group(1))
                
        if "ITERATION" in line and "SUMMARY" not in line:
             m = re.search(r'ITERATION\s+(\d+)', line)
             if m:
                 it = int(m.group(1))
                 if it > current_iters:
                     current_iters = it

        if "INCREMENT" in line and "COMPLETED" in line:
             iteration_history.append({'increment': current_inc, 'iterations': current_iters})
             current_iters = 0
             
    return iteration_history

def main():
    print("Reading simulation data...")
    inc_uel = parse_dat_file("../sim/cantilever_uel.dat")
    iter_uel = parse_msg_file("../sim/cantilever_uel.msg")
    
    inc_fd = parse_dat_file("../sim/cantilever_fd.dat")
    iter_fd = parse_msg_file("../sim/cantilever_fd.msg")
    
    inc_hc = parse_dat_file("../sim/cantilever_complex.dat")
    iter_hc = parse_msg_file("../sim/cantilever_complex.msg")
    
    inc_ex = parse_dat_file("../sim/cantilever_exact.dat")
    iter_ex = parse_msg_file("../sim/cantilever_exact.msg")
    
    print(f"Loaded Hyperdual: {len(inc_uel)} increments")
    print(f"Loaded FD: {len(inc_fd)} increments")
    print(f"Loaded Complex: {len(inc_hc)} increments")
    print(f"Loaded Exact: {len(inc_ex)} increments")
    
    has_hc = len(inc_hc) > 0
    has_ex = len(inc_ex) > 0
    
    # Determine frames (using smallest available count to sync)
    lens = [len(inc_uel), len(inc_fd)]
    if has_hc: lens.append(len(inc_hc))
    if has_ex: lens.append(len(inc_ex))
    num_frames = min(lens)
    
    if num_frames == 0:
        print("Error: No valid frames found.")
        sys.exit(0)
        
    # Create Figure
    fig = make_subplots(
        rows=2, cols=1,
        row_heights=[0.7, 0.3],
        specs=[[{'type': 'scene'}], [{'type': 'xy'}]],
        subplot_titles=("Beam Deflection (Blue:HD, Red:FD, Green:CS)", "Tip Deflection Error (vs Exact)")
    )
    
    # Calculate Axis Ranges
    all_x = []
    # Just gather last frame data for bounds
    x_l, y_l, _ = interpolate_beam(inc_uel[-1])
    all_x.extend(x_l)
    
    x_min, x_max = 0, 10
    if all_x: x_max = max(all_x)
    
    # Y-range typically -1 to 0 for cantilever down, or +1
    y_min, y_max = -2.0, 0.5 
    
    # Calculate Errors
    err_hd = []
    err_fd = []
    err_cs = []
    
    inc_indices = list(range(num_frames))
    
    min_err_floor = 1e-18
    
    for k in range(num_frames):
        # We need the Exact solution to compare against
        if not has_ex:
            err_hd.append(None)
            err_fd.append(None)
            err_cs.append(None)
            continue
            
        # Get Exact Tip Displacement (U2) for this increment
        ex_nodes = inc_ex[k]['nodes']
        if not ex_nodes:
             err_hd.append(None); err_fd.append(None); err_cs.append(None); continue
             
        tip_id = max(ex_nodes.keys())
        u2_ex = ex_nodes[tip_id]['disp'][1]
        
        # HD
        hd_nodes = inc_uel[k]['nodes']
        if tip_id in hd_nodes:
            u2_hd = hd_nodes[tip_id]['disp'][1]
            diff = abs(u2_hd - u2_ex)
            if diff < min_err_floor: diff = min_err_floor
            err_hd.append(diff)
        else:
            err_hd.append(None)
            
        # FD
        fd_nodes = inc_fd[k]['nodes']
        if tip_id in fd_nodes:
            u2_fd = fd_nodes[tip_id]['disp'][1]
            diff = abs(u2_fd - u2_ex)
            if diff < min_err_floor: diff = min_err_floor
            err_fd.append(diff)
        else:
            err_fd.append(None)
            
        # CS
        if has_hc:
            hc_nodes = inc_hc[k]['nodes']
            if tip_id in hc_nodes:
                u2_hc = hc_nodes[tip_id]['disp'][1]
                diff = abs(u2_hc - u2_ex)
                if diff < min_err_floor: diff = min_err_floor
                err_cs.append(diff)
            else:
                err_cs.append(None)
    
    # Initial Data (k=0)
    x0_u, y0_u, z0_u = interpolate_beam(inc_uel[0])
    x0_f, y0_f, z0_f = interpolate_beam(inc_fd[0])
    
    trace_indices = [0, 1]
    fig.add_trace(go.Scatter3d(x=x0_u, y=y0_u, z=z0_u, mode='lines+markers', line=dict(color='blue', width=6), marker=dict(size=3), name='Hyperdual'), row=1, col=1)
    fig.add_trace(go.Scatter3d(x=x0_f, y=y0_f, z=z0_f, mode='lines+markers', line=dict(color='red', width=6, dash='longdash'), marker=dict(size=3, symbol='circle-open'), name='Finite Diff'), row=1, col=1)
    
    if has_hc:
        x0_h, y0_h, z0_h = interpolate_beam(inc_hc[0])
        fig.add_trace(go.Scatter3d(x=x0_h, y=y0_h, z=z0_h, mode='lines+markers', line=dict(color='green', width=4, dash='dot'), marker=dict(size=3, symbol='diamond'), name='Complex'), row=1, col=1)
        trace_indices.append(len(trace_indices)) # Auto append

    if has_ex:
        x0_e, y0_e, z0_e = interpolate_beam(inc_ex[0])
        fig.add_trace(go.Scatter3d(x=x0_e, y=y0_e, z=z0_e, mode='lines+markers', line=dict(color='orange', width=4), marker=dict(size=3, symbol='cross'), name='Exact'), row=1, col=1)
        trace_indices.append(len(trace_indices))

    # Error Plots
    # Use log scale for Y because errors can vary by orders of magnitude (1e-16 vs 1e-6)
    fig.add_trace(go.Scatter(x=inc_indices, y=err_hd, mode='lines+markers', name='HD Error', line=dict(color='blue')), row=2, col=1)
    fig.add_trace(go.Scatter(x=inc_indices, y=err_fd, mode='lines+markers', name='FD Error', line=dict(color='red')), row=2, col=1)
    if has_hc:
        fig.add_trace(go.Scatter(x=inc_indices, y=err_cs, mode='lines+markers', name='CS Error', line=dict(color='green')), row=2, col=1)
        
    fig.update_yaxes(type="log", title="Tip Error (Log Scale)", row=2, col=1)
    fig.update_xaxes(title="Increment", row=2, col=1)
    
    # Frames for Animation
    frames = []
    for k in range(num_frames):
        x_u, y_u, z_u = interpolate_beam(inc_uel[k])
        x_f, y_f, z_f = interpolate_beam(inc_fd[k])
        
        frame_data = [
            go.Scatter3d(x=x_u, y=y_u, z=z_u),
            go.Scatter3d(x=x_f, y=y_f, z=z_f)
        ]
        
        if has_hc:
            x_h, y_h, z_h = interpolate_beam(inc_hc[k])
            frame_data.append(go.Scatter3d(x=x_h, y=y_h, z=z_h))

        if has_ex:
            x_e, y_e, z_e = interpolate_beam(inc_ex[k])
            frame_data.append(go.Scatter3d(x=x_e, y=y_e, z=z_e))
            
        frames.append(go.Frame(data=frame_data, name=f"Inc{k}", traces=trace_indices))
        
    fig.frames = frames

    # Sliders and Buttons
    sliders = [dict(
        steps=[dict(
            method='animate',
            args=[[f'Inc{k}'], dict(mode='immediate', frame=dict(duration=0, redraw=True), transition=dict(duration=0))],
            label=f'{k}'
        ) for k in range(num_frames)],
        transition=dict(duration=0),
        x=0.1,
        y=0,
        currentvalue=dict(font=dict(size=12), prefix='Increment: ', visible=True, xanchor='center'),
        len=0.9
    )]

    fig.update_layout(
        scene=dict(
            xaxis=dict(title='X (Length)', range=[0, 11]),
            yaxis=dict(title='U2 (Deflection)', range=[-4.0, 1.0]),
            zaxis=dict(title='Z', range=[-1.0, 1.0]),
            aspectmode='manual',
            aspectratio=dict(x=4, y=2, z=0.5)
        ),
        updatemenus=[{
            'type': 'buttons',
            'showactive': False,
            'y': 0,
            'x': 0,
            'xanchor': 'right',
            'yanchor': 'top',
            'pad': dict(t=0, r=10),
            'buttons': [{
                'label': 'Play',
                'method': 'animate',
                'args': [None, dict(frame=dict(duration=100, redraw=True), fromcurrent=True, transition=dict(duration=0))]
            }, {
                'label': 'Pause',
                'method': 'animate',
                'args': [[None], dict(mode='immediate', frame=dict(duration=0, redraw=False), transition=dict(duration=0))]
            }]
        }],
        sliders=sliders,
        title="Beam Deflection History"
    )
    
    print("Writing output...")
    fig.write_html("../sim/beam_comparison.html")
    print("Done.")

if __name__ == "__main__":
    main()
de
