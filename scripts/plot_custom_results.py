import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import re

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    sim_dir = os.path.join(script_dir, '../sim')
    runs_dir = os.path.join(sim_dir, 'runs')
    output_plot = os.path.join(sim_dir, 'plots', 'error_vs_nelems_custom.png')
    
    if not os.path.exists(os.path.dirname(output_plot)):
        os.makedirs(os.path.dirname(output_plot))

    data = []
    
    # Iterate through run directories
    # Looking for tags like: custom_ne{ne}_f{force}_inc{ninc}
    # But specifically force 10000
    
    if not os.path.exists(runs_dir):
        print(f"Runs directory {runs_dir} does not exist.")
        return

    print("collecting data...")
    for tag in os.listdir(runs_dir):
        run_path = os.path.join(runs_dir, tag)
        if not os.path.isdir(run_path):
            continue
            
        # Check if it matches our pattern
        # custom_ne100_f10000_inc10
        match = re.search(r'ne(\d+)_f10000_inc', tag) 
        # Note: Depending on my run script, I used "custom_" prefix or not. 
        # In run_custom_sweep.py I used: tag = f"custom_ne{ne}_f{abs(int(force))}_inc{ninc}"
        # So it should be custom_ne...
        
        if not match:
            continue
            
        nelems = int(match.group(1))
        
        hist_file = os.path.join(run_path, 'tip_history.csv')
        if not os.path.exists(hist_file):
            continue
            
        try:
            df = pd.read_csv(hist_file)
            # Use the last value (end of simulation) or max error?
            # Usually end of simulation for static analysis.
            row = df.iloc[-1]
            
            # Using max error over the history is also valid.
            # Let's stick effectively to max error as static load increases monotonically.
            # In plot_sweep_results.py it used df['err_hd'].max()
            
            err_hd = df['err_hd'].max()
            err_fd = df['err_fd'].max()
            err_cs = df['err_cs'].max()
            err_ot = df['err_otis'].max()
            
            data.append({
                'nelems': nelems,
                'HyperDual': err_hd,
                'FiniteDiff': err_fd,
                'ComplexStep': err_cs,
                'Otis': err_ot
            })
        except Exception as e:
            print(f"Error reading {tag}: {e}")
            
    if not data:
        print("No matching data found.")
        return
        
    df_results = pd.DataFrame(data)
    df_results = df_results.sort_values('nelems')
    
    print(f"Plotting results for {len(df_results)} simulation points...")
    
    plt.figure(figsize=(10, 8))
    
    # Plot each method
    plt.plot(df_results['nelems'], df_results['HyperDual'], 'o-', label='HyperDual')
    plt.plot(df_results['nelems'], df_results['FiniteDiff'], 'x-', label='FiniteDiff')
    plt.plot(df_results['nelems'], df_results['ComplexStep'], 's-', label='ComplexStep')
    plt.plot(df_results['nelems'], df_results['Otis'], 'd-', label='Otis')
    
    plt.xlabel('Number of Elements')
    plt.ylabel('Error (vs Exact UEL)')
    plt.title('Error vs Number of Elements (Force = 10000)')
    plt.yscale('log')
    plt.grid(True, which="both", ls="-")
    plt.legend()
    
    plt.savefig(output_plot)
    print(f"Plot saved to {output_plot}")

if __name__ == "__main__":
    main()
