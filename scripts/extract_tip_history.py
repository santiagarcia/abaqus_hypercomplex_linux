#!/usr/bin/env python3
import sys
import os
import csv
import glob
import numpy as np
from odbAccess import openOdb

def extract_for_tag(run_dir, tag):
    """
    Extracts tip displacement history for all methods in the run_dir with the given tag.
    Writes to 'tip_history.csv' in the run_dir.
    """
    
    # Define expected ODB files
    methods = {
        'exact': 'cantilever_{}_exact.odb'.format(tag),
        'hd': 'cantilever_{}_uel.odb'.format(tag),
        'fd': 'cantilever_{}_fd.odb'.format(tag),
        'cs': 'cantilever_{}_complex.odb'.format(tag),
        'otis': 'cantilever_{}_otis.odb'.format(tag)
    }
    
    print("Extracting from " + run_dir)
    
    # Data Storage
    results = {}
    
    # Iterate over methods
    max_frames = 0
    
    for method, filename in methods.items():
        odb_path = os.path.join(run_dir, filename)
        
        # Initialize list for this method
        results[method] = []
        
        if not os.path.exists(odb_path):
            print("Warning: ODB not found: " + odb_path)
            continue
            
        try:
            odb = openOdb(odb_path, readOnly=True)
            
            # Assume first step
            if len(odb.steps) == 0:
                print("Warning: No steps in " + filename)
                odb.close()
                continue
                
            step_name = odb.steps.keys()[0]
            step = odb.steps[step_name]
            
            # Find Tip Node (Highest X)
            # We look at the first frame's instance
            root_assembly = odb.rootAssembly
            # Try to find the instance name. Usually 'PART-1-1' or similar.
            inst_name = root_assembly.instances.keys()[0]
            instance = root_assembly.instances[inst_name]
            
            # Find node with max X
            max_x = -1e20
            tip_label = -1
            
            # Scanning all nodes can be slow for large models, but here we have few elements
            for node in instance.nodes:
                if node.coordinates[0] > max_x:
                    max_x = node.coordinates[0]
                    tip_label = node.label
            
            # print("Tip Node for " + method + ": " + str(tip_label) + " at X=" + str(max_x))
            
            # Extract U2 history
            # We can use history output if defined, or field output. 
            # Field output is safer if history wasn't requested.
            
            frames = step.frames
            if len(frames) > max_frames:
                max_frames = len(frames)
                
            for i, frame in enumerate(frames):
                # if i == 0: continue # Skip frame 0 if desired, but usually we want it (t=0, u=0)
                
                # fieldOutputs['U'] might not exist in frame 0 if not requested, but usually zero-results are there or implied
                try:
                    u_field = frame.fieldOutputs['U']
                    # Get subset for tip node
                    # Note: getSubset returns a FieldOutput object
                    # We need to match the node label.
                    # region=instance.nodes[tip_label-1] might work if indices align, 
                    # but simpler to use getSubset(region=...)
                    
                    # Create a node set for the tip if possible? No, can't modify ODB.
                    # Just iterate values? 
                    # For performance on large models this is bad, but for beam elements it is fine.
                    
                    # Optimization: getSubset by nodeLabel
                    # In Abaqus Python, we can't easily select by label in getSubset without a Set.
                    # So we iterate values.
                    
                    val = None
                    # Attempt sub-setting is tricky without a nodeSet.
                    # Let's iterate. 
                    for v in u_field.values:
                         if v.nodeLabel == tip_label:
                             val = v
                             break
                    
                    if val:
                        u2 = val.data[1] # Y-displacement
                        results[method].append(u2)
                    else:
                        results[method].append(None)
                        
                except Exception as e:
                    # Might happen in frame 0 if U not defined
                    results[method].append(0.0)
            
            odb.close()
            
        except Exception as e:
            print("Error reading " + filename + ": " + str(e))
    
    # Write CSV
    csv_path = os.path.join(run_dir, 'tip_history.csv')
    with open(csv_path, 'w') as f:
        writer = csv.writer(f)
        header = ['frame', 'u2_exact', 'u2_hd', 'u2_fd', 'u2_cs', 'u2_otis', 'err_hd', 'err_fd', 'err_cs', 'err_otis']
        writer.writerow(header)
        
        for i in range(max_frames):
            row = [i]
            
            # Get Displacements
            u_ex = results.get('exact')[i] if 'exact' in results and i < len(results['exact']) else None
            u_hd = results.get('hd')[i] if 'hd' in results and i < len(results['hd']) else None
            u_fd = results.get('fd')[i] if 'fd' in results and i < len(results['fd']) else None
            u_cs = results.get('cs')[i] if 'cs' in results and i < len(results['cs']) else None
            u_ot = results.get('otis')[i] if 'otis' in results and i < len(results['otis']) else None
            
            row.extend([u_ex, u_hd, u_fd, u_cs, u_ot])
            
            # Calculate Errors
            if u_ex is not None:
                err_hd = abs(u_hd - u_ex) if u_hd is not None else None
                err_fd = abs(u_fd - u_ex) if u_fd is not None else None
                err_cs = abs(u_cs - u_ex) if u_cs is not None else None
                err_ot = abs(u_ot - u_ex) if u_ot is not None else None
            else:
                err_hd, err_fd, err_cs, err_ot = None, None, None, None
                
            row.extend([err_hd, err_fd, err_cs, err_ot])
            writer.writerow(row)
            
    print("Saved tip history to " + csv_path)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: abaqus python extract_tip_history.py <run_dir> <tag>")
        sys.exit(1)
        
    run_dir = sys.argv[1]
    tag = sys.argv[2]
    
    extract_for_tag(run_dir, tag)
