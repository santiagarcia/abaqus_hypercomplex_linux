import subprocess
import os
import sys
import time

def main():
    # Define Sweep Parameters
    nelems_list = [1, 2, 4, 8, 20]
    force_list = [-100.0]
    nincs_list = [10]
    
    # Or simplified for testing:
    # nelems_list = [1, 5]
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    run_script = os.path.join(script_dir, "run_comparison.sh")
    extract_script = os.path.join(script_dir, "extract_tip_history.py")
    viz_script = os.path.join(script_dir, "visualize_comparison.py")
    plot_script = os.path.join(script_dir, "plot_sweep_results.py")
    
    repo_root = os.path.dirname(script_dir)
    sim_dir = os.path.join(repo_root, "sim")
    
    print("Starting Parametric Sweep...")
    
    for ne in nelems_list:
        for f in force_list:
            for ninc in nincs_list:
                tag = f"ne{ne}_f{abs(int(f))}_inc{ninc}"
                
                print(f"--- Running Case: {tag} (Ne={ne}, F={f}, Ninc={ninc}) ---")
                
                # 1. Run Simulation Comparison
                # Arguments: NELEMS FORCE NUM_INCS TAG
                cmd_run = ["bash", run_script, str(ne), str(f), str(ninc), tag]
                
                try:
                    subprocess.run(cmd_run, check=True)
                except subprocess.CalledProcessError as e:
                    print(f"Error running simulation for {tag}: {e}")
                    continue
                
                # 2. Extract History (Abaqus Python)
                # Usage: abaqus python extract_tip_history.py <run_dir> <tag>
                run_dir = os.path.join(sim_dir, "runs", tag)
                cmd_extract = ["abaqus", "python", extract_script, run_dir, tag]
                
                try:
                    print(f"Extracting history for {tag}...")
                    subprocess.run(cmd_extract, check=True)
                except subprocess.CalledProcessError as e:
                     print(f"Error extracting history for {tag}: {e}")
                     continue

                # 3. Generate Visualization (Plotly HTML)
                # Usage: python3 visualize_comparison.py <run_dir> <tag>
                cmd_viz = ["python3", viz_script, run_dir, tag]
                try:
                    print(f"Generating visualization for {tag}...")
                    subprocess.run(cmd_viz, check=True)
                except subprocess.CalledProcessError as e:
                    print(f"Error generating visualization for {tag}: {e}")
                     
    print("--- Sweep Complete ---")
    
    # 3. Generate Summary Plots
    print("Generating Plots...")
    try:
        subprocess.run(["python3", plot_script], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error plotting results: {e}")
        
    print("Done.")

if __name__ == "__main__":
    main()
