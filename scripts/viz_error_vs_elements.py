import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import re
import numpy as np

def main():
    sim_dir = os.path.join(os.path.dirname(__file__), '../sim')
    runs_dir = os.path.join(sim_dir, 'runs')
    plot_file = os.path.join(sim_dir, 'plots', 'error_vs_elements.png')
    
    if not os.path.exists(os.path.dirname(plot_file)):
        os.makedirs(os.path.dirname(plot_file))

    # Pattern to match: lin_ne{N}_f10000_inc10
    pattern = os.path.join(runs_dir, 'lin_ne*_f10000_inc10')
    dirs = glob.glob(pattern)
    
    data = []
    
    print(f"Found {len(dirs)} directories matching pattern.")
    
    for d in dirs:
        dir_name = os.path.basename(d)
        # Extract number of elements
        match = re.search(r'ne(\d+)_', dir_name)
        if match:
            nelems = int(match.group(1))
        else:
            continue
            
        # Old history file
        # hist_file = os.path.join(d, 'tip_history.csv')
        
        # New MSE file
        mse_file = os.path.join(d, 'mse_error.csv')
        
        if not os.path.exists(mse_file):
            print(f"Missing MSE file in {dir_name}")
            continue
            
        try:
            df = pd.read_csv(mse_file)
            if df.empty:
                continue
                
            # Structure: method, mse, max_error
            # Transpose to single row
            row_data = {'nelems': nelems}
            
            for index, row in df.iterrows():
                method = row['method']
                mse = float(row['mse'])
                if mse < 1e-13:
                    mse = 0.0
                row_data[f'err_{method}'] = mse # Using MSE as the error metric
            
            data.append(row_data)
        except Exception as e:
            print(f"Error reading {mse_file}: {e}")
            
    if not data:
        print("No data found.")
        return

    df_results = pd.DataFrame(data)
    df_results = df_results.sort_values('nelems')
    
    print("Data collected:")
    print(df_results.head())
    
    # Plotting
    plt.figure(figsize=(10, 6))
    
    # Plot each method present
    if 'err_hd' in df_results.columns:
        plt.plot(df_results['nelems'], df_results['err_hd'], 'o-', label='HyperDual (UEL)')
    if 'err_fd' in df_results.columns:
        plt.plot(df_results['nelems'], df_results['err_fd'], 's-', label='Finite Difference')
    if 'err_cs' in df_results.columns:
        plt.plot(df_results['nelems'], df_results['err_cs'], '^-', label='Complex Step')
    if 'err_otis' in df_results.columns:
        plt.plot(df_results['nelems'], df_results['err_otis'], 'd-', label='Otis (HyperDual)')

    plt.xlabel('Number of Elements', fontsize=12)
    plt.ylabel('Final Error (Relative to Exact)', fontsize=12)
    plt.title('Error vs Number of Elements (Force = 10000)', fontsize=14)
    
    # Use linear scale
    plt.yscale('linear')
    
    # Customize grid and ticks
    plt.grid(True, which="major", ls="-", alpha=0.8)
    plt.grid(True, which="minor", ls=":", alpha=0.5)
    
    # Increase tick label size
    plt.tick_params(axis='both', which='major', labelsize=10)
    
    plt.legend(fontsize=10)
    
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {plot_file}")

if __name__ == "__main__":
    main()
