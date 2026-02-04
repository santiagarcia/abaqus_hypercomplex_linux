import sys
import os
import csv
import numpy as np
from odbAccess import openOdb

def extract_mse(run_dir, tag):
    # Dict of methods mapping to ODB filenames
    methods = {
        'exact': 'cantilever_{}_exact.odb'.format(tag),
        'hd': 'cantilever_{}_uel.odb'.format(tag),
        'fd': 'cantilever_{}_fd.odb'.format(tag),
        'cs': 'cantilever_{}_complex.odb'.format(tag),
        'otis': 'cantilever_{}_otis.odb'.format(tag)
    }
    
    # Store displacements for all nodes
    # Structure: data[method] = { node_label: displacement_u2 }
    data = {}
    
    # We will compute MSE at the last frame
    
    for method, filename in methods.items():
        odb_path = os.path.join(run_dir, filename)
        if not os.path.exists(odb_path):
            continue
            
        try:
            odb = openOdb(odb_path, readOnly=True)
            if len(odb.steps) == 0:
                odb.close()
                continue
                
            last_step_key = list(odb.steps.keys())[-1]
            last_step = odb.steps[last_step_key]
            
            if len(last_step.frames) == 0:
                odb.close()
                continue
            
            last_frame = last_step.frames[-1]
            
            # Get Field Output U (Displacement)
            if 'U' not in last_frame.fieldOutputs:
                odb.close()
                continue
                
            u_field = last_frame.fieldOutputs['U']
            
            node_data = {}
            for v in u_field.values:
                # v.data is tuple (Ux, Uy, Uz) for 3D or 2D. 
                # Cantilever input defines 3 coordinates.
                # Displacement matches.
                # We can store the full vector and compute distance, or just U2.
                # User asked for "positions", implying vector magnitude or just components.
                # Let's use the magnitude of the difference vector for safety, 
                # or just U2 since it's a cantilever loaded in Y.
                # Let's use U2 to be consistent with previous tip history which focused on U2. 
                # Or even better: MSE of the scalar displacement field U2.
                
                val = v.data[1] # U2
                node_data[v.nodeLabel] = val
                
            data[method] = node_data
            odb.close()
            
        except Exception as e:
            print("Error processing {}: {}".format(method, e))
            
    # Calculate MSE relative to Exact
    if 'exact' not in data:
        print("Exact solution not found in {}".format(run_dir))
        return
        
    exact_data = data['exact']
    results = {}
    
    # Output file
    csv_path = os.path.join(run_dir, 'mse_error.csv')
    
    with open(csv_path, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['method', 'mse', 'max_error'])
        
        for method in ['hd', 'fd', 'cs', 'otis']:
            if method not in data:
                 writer.writerow([method, 'NaN', 'NaN'])
                 continue
                
            method_data = data[method]
            sq_diff_sum = 0.0
            max_diff = 0.0
            count = 0
            
            for label, u_exact in exact_data.items():
                if label in method_data:
                    u_method = method_data[label]
                    diff = u_method - u_exact
                    sq_diff_sum += diff * diff
                    max_diff = max(max_diff, abs(diff))
                    count += 1
            
            if count > 0:
                mse = sq_diff_sum / count
                results[method] = mse
                writer.writerow([method, mse, max_diff])
            else:
                writer.writerow([method, 'NaN', 'NaN'])
                
    print("MSE calculation complete for {}".format(tag))

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: abaqus python extract_mse_error.py <run_dir> <tag>")
        sys.exit(1)
        
    run_dir = sys.argv[1]
    tag = sys.argv[2]
    extract_mse(run_dir, tag)
