import subprocess
import os
import sys

def main():
    # User requirements
    # Linear steps of 10 from 2 to 1000: 2, 12, 22, ...
    nelems_list = list(range(2, 1001, 10))
    
    force = -10000.0  # Force of 10000 (downward)
    ninc = 10         # Default number of increments
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    run_script = os.path.join(script_dir, "run_comparison.sh")
    extract_script = os.path.join(script_dir, "extract_tip_history.py")
    extract_mse_script = os.path.join(script_dir, "extract_mse_error.py")
    
    repo_root = os.path.dirname(script_dir)
    sim_dir = os.path.join(repo_root, "sim")
    
    print(f"Starting Custom Sweep for {len(nelems_list)} cases: {nelems_list}")
    
    for ne in nelems_list:
        # Construct a unique tag for this run
        tag = f"lin_ne{ne}_f{abs(int(force))}_inc{ninc}"
        
        print(f"--- Running Case: {tag} (Ne={ne}, F={force}) ---")
        
        # 1. Run Simulation Comparison
        # Arguments: NELEMS FORCE NUM_INCS TAG
        cmd_run = ["bash", run_script, str(ne), str(force), str(ninc), tag]
        
        # Check if run already exists to save time? 
        # User said "run ... again" implicitly by asking for new results, but maybe we can just re-extract if ODBs exist
        # However, to be safe and simple, let's re-run or assume run_script handles skipping if robust (it doesn't seem to check).
        # But wait, the previous run_custom_sweep tool call was cancelled/failed? 
        # No, the PREVIOUS one succeeded with 4 cases.
        # I can keep the run command.
        
        try:
            subprocess.run(cmd_run, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running simulation for {tag}: {e}")
            continue
        
        # 2. Extract History (Abaqus Python)
        run_dir = os.path.join(sim_dir, "runs", tag)
        cmd_extract = ["abaqus", "python", extract_script, run_dir, tag]
        
        try:
            print(f"Extracting history for {tag}...")
            # We capture output to avoid cluttering the terminal too much
            subprocess.run(cmd_extract, check=True, stdout=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
             print(f"Error extracting history for {tag}: {e}")
             continue

        # 3. Extract MSE Error (Abaqus Python)
        cmd_mse = ["abaqus", "python", extract_mse_script, run_dir, tag]
        try:
            print(f"Extracting MSE Error for {tag}...")
            subprocess.run(cmd_mse, check=True, stdout=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"Error extracting MSE for {tag}: {e}")
             
    print("--- Sweep Complete ---")

if __name__ == "__main__":
    main()
