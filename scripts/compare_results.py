import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
import re
import os

# Configuration
DAT_UEL = os.path.join(os.path.dirname(__file__), '../sim/cantilever_uel.dat')
DAT_FD  = os.path.join(os.path.dirname(__file__), '../sim/cantilever_fd.dat')
OUTPUT_HTML = os.path.join(os.path.dirname(__file__), '../sim/comparison.html')

def parse_dat_file(filepath):
    """Parses Abaqus .dat file for Node Output."""
    if not os.path.exists(filepath):
        print(f"Error: {filepath} not found.")
        return [], 0.0

    with open(filepath, 'r') as f:
        lines = f.readlines()

    increments = []
    current_increment = {}
    
    cpu_time = 0.0
    
    re_inc_summary = re.compile(r'INCREMENT\s+(\d+)\s+SUMMARY')
    re_time_comp = re.compile(r'TOTAL TIME COMPLETED\s+([\d\.E\+\-]+)')
    re_node_line = re.compile(r'^\s*(\d+)\s+([\d\.E\+\-]+)\s+([\d\.E\+\-]+)\s+([\d\.E\+\-]+)\s+([\d\.E\+\-]+)\s+([\d\.E\+\-]+)\s+([\d\.E\+\-]+)')
    re_cpu = re.compile(r'TOTAL CPU TIME \(SEC\)\s+=\s+([\d\.E\+\-]+)')

    reading_nodes = False
    
    inc_data = None 
    
    for line in lines:
        if "TOTAL CPU TIME" in line:
            m = re_cpu.search(line)
            if m:
                # Capture the last one found (Job Time Summary is usually at end)
                try:
                    cpu_time = float(m.group(1))
                except:
                    pass

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
                     inc_data['nodes'][node_id] = vals
                 elif "MAXIMUM" in line:
                     reading_nodes = False
    
    if inc_data and 'nodes' in inc_data and len(inc_data['nodes']) > 0:
         increments.append(inc_data)
         
    return increments, cpu_time

def get_tip_displacement(increments):
    # Get DOF 2 (v) at Node 2 for each increment
    times = []
    disps = []
    for inc in increments:
        if 2 in inc['nodes']:
            times.append(inc['time'])
            disps.append(inc['nodes'][2][1]) # U2
    return times, disps

def main():
    print("Parsing UEL Results...")
    inc_uel, time_uel = parse_dat_file(DAT_UEL)
    print("Parsing FD Results...")
    inc_fd, time_fd = parse_dat_file(DAT_FD)
    
    t_uel, u_uel = get_tip_displacement(inc_uel)
    t_fd, u_fd = get_tip_displacement(inc_fd)
    
    print(f"UEL CPU Time: {time_uel} s")
    print(f"FD CPU Time:  {time_fd} s")
    
    fig = go.Figure()
    
    fig.add_trace(go.Scatter(x=t_uel, y=u_uel, mode='lines+markers', name=f'Hyperdual (CPU: {time_uel}s)'))
    fig.add_trace(go.Scatter(x=t_fd, y=u_fd, mode='lines+markers', name=f'Finite Diff (CPU: {time_fd}s)', line=dict(dash='dash')))
    
    fig.update_layout(title='Tip Deflection Comparison: Hyperdual vs FD',
                      xaxis_title='Time',
                      yaxis_title='Tip Displacement (U2)')
    
    print(f"Writing comparison to {OUTPUT_HTML}")
    fig.write_html(OUTPUT_HTML)

if __name__ == "__main__":
    main()
