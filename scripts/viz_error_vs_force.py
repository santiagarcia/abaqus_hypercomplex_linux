import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import re
import subprocess
import sys

def main():
    sim_dir = os.path.join(os.path.dirname(__file__), '../sim')
    runs_dir = os.path.join(sim_dir, 'runs')
    scripts_dir = os.path.dirname(os.path.abspath(__file__))
    extract_mse_script = os.path.join(scripts_dir, "extract_mse_error.py")
    
    plot_file = os.path.join(sim_dir, 'plots', 'error_vs_force.png')
    
    if not os.path.exists(os.path.dirname(plot_file)):
        os.makedirs(os.path.dirname(plot_file))

    # Pattern: ne32_f*_inc10
    # Note: excluding 'lin_' or 'exp_' prefixes for now, assuming standard runs have the data
    # The grep showed "ne32_f...." without prefix, and also "lin_ne32_f..."
    # I'll generally look for *ne32_f*_inc10 and parse carefully.
    
    pattern = os.path.join(runs_dir, '*ne32_f*_inc10')
    dirs = glob.glob(pattern)
    
    data = []
    
    print(f"Found {len(dirs)} directories matching pattern.")
    
    for d in dirs:
        dir_name = os.path.basename(d)
        
        # Parse Force
        # Expected formats: ne32_f10000_inc10, lin_ne32_f10000_inc10
        # Regex to find f{NUMBER}
        match = re.search(r'_f(\d+)_', dir_name)
        if match:
            force = float(match.group(1))
        else:
            continue
            
        # Filter: We want exactly ne32.
        # Check ne
        ne_match = re.search(r'ne(\d+)_', dir_name)
        if not ne_match or int(ne_match.group(1)) != 32:
            continue
            
        # Check inc
        if '_inc10' not in dir_name:
            continue

        mse_file = os.path.join(d, 'mse_error.csv')
        
        # Ensure MSE file exists
        if not os.path.exists(mse_file):
            print(f"Generating MSE for {dir_name}...")
            tag = dir_name # Folder name is the tag
            try:
                # Call abaqus python ...
                # We need to construct absolute path suitable for command line
                cmd = ["abaqus", "python", extract_mse_script, d, tag]
                subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
            except subprocess.CalledProcessError as e:
                print(f"Failed to extract MSE for {tag}: {e}")
                continue
        
        # Read Data
        if os.path.exists(mse_file):
            try:
                df = pd.read_csv(mse_file)
                if df.empty:
                    continue
                    
                row_data = {'force': force}
                
                for index, row in df.iterrows():
                    method = row['method']
                    mse = float(row['mse'])
                    if mse < 1e-13:
                        mse = 0.0
                    row_data[f'err_{method}'] = mse
                
                data.append(row_data)
            except Exception as e:
                print(f"Error reading {mse_file}: {e}")

    if not data:
        print("No data found.")
        return

    df_results = pd.DataFrame(data)
    df_results = df_results.sort_values('force')
    
    print("Data collected:")
    print(df_results.head())
    
    # Validation: Filter out duplicates (e.g. if lin_ne32 and ne32 both exist with same force)
    # Prefer non-prefixed if available? Or just take one.
    df_results = df_results.drop_duplicates(subset=['force'], keep='last')
    
    # Plotting
    plt.figure(figsize=(10, 6))
    
    if 'err_hd' in df_results.columns:
        plt.plot(df_results['force'], df_results['err_hd'], 'o-', label='HyperDual (UEL)')
    if 'err_fd' in df_results.columns:
        plt.plot(df_results['force'], df_results['err_fd'], 's-', label='Finite Difference')
    if 'err_cs' in df_results.columns:
        plt.plot(df_results['force'], df_results['err_cs'], '^-', label='Complex Step')
    if 'err_otis' in df_results.columns:
        plt.plot(df_results['force'], df_results['err_otis'], 'd-', label='Otis (HyperDual)')

    plt.xlabel('Force Magnitude', fontsize=12)
    plt.ylabel('MSE Error (Relative to Exact)', fontsize=12)
    plt.title('Error vs Force (Elements = 32)', fontsize=14)
    
    # Scale choice: Force spans 100 to 100,000. Login X is good.
    # Error might span orders.
    # Let's try log-log.
    plt.xscale('log')
    plt.yscale('symlog', linthresh=1e-13)
    
    plt.grid(True, which="major", ls="-", alpha=0.8)
    plt.grid(True, which="minor", ls=":", alpha=0.5)
    
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.legend(fontsize=10)
    
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {plot_file}")

if __name__ == "__main__":
    main()
