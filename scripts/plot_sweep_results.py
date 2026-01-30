import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import glob

def main():
    sim_dir = os.path.join(os.path.dirname(__file__), '../sim')
    runs_dir = os.path.join(sim_dir, 'runs')
    plots_dir = os.path.join(sim_dir, 'plots')
    
    if not os.path.exists(plots_dir):_timings['force']==fv)]
        #        if subset.empty: continue
               
        #        plt.figure(figsize=(10,6))
               
        os.makedirs(plots_dir)
        
    print("Generating plots in " + plots_dir)
    
    # 1. Plot individual error histories
    # Find all tip_history.csv
    # history_files = glob.glob(os.path.join(runs_dir, '*', 'tip_history.csv'))
    
    # 2. Gather Summary Data (Timings and Max Errors)
    timings_file = os.path.join(sim_dir, 'timings.csv')
    if os.path.exists(timings_file):
        df_timings = pd.read_csv(timings_file)
        
        # CPU Time vs Elements (for different methods)
        # Assuming Force/Nincs constant or grouping
        
        # Filter for typical case (e.g. nincs=10, force=-100) if multiple exist, 
        # or just plot everything if distinct.
        
        # Check unique params
        unique_nincs = df_timings['nincs'].unique()
        unique_force = df_timings['force'].unique()
        
        # Plot Time vs Nelems for each method
        # for ni in unique_nincs:
        #    for fv in unique_force:
        #        subset = df_timings[(df_timings['nincs']==ni) & (df_timings['force']==fv)]
        #        if subset.empty: continue
               
        #        plt.figure(figsize=(10,6))
               
        #        for method in subset['method'].unique():
        #             data = subset[subset['method'] == method].sort_values('nelems')
        #             plt.plot(data['ne  python3 scripts/plot_sweep_results.pylems'], data['cpu_time_seconds'], marker='o', label=method)
                    
        #        plt.xlabel('Number of Elements')
        #        plt.ylabel('CPU Time (s)')
        #        plt.title(f'Runtime Scaling (Incs={ni}, Force={fv})')
        #        plt.legend()
        #        plt.grid(True)
        #        plt.xscale('log')
        #        plt.yscale('log')
        #        plt.savefig(os.path.join(plots_dir, f'time_vs_nelems_ni{ni}_f{fv}.png'))
        #        plt.close()
        
        # Simpler approach: Just plot all 'method' groups if params vary, maybe facet?
        # Let's pivot
        try:
             # Pivot to get mean time for each method/nelem combo (in case of dupes)
             pivot = df_timings.pivot_table(index='nelems', columns='method', values='cpu_time_seconds', aggfunc='mean')
             
             plt.figure(figsize=(10,6))
             for col in pivot.columns:
                 plt.plot(pivot.index, pivot[col], marker='o', label=col)
             
             plt.xlabel('Number of Elements')
             plt.ylabel('CPU Time (s)')
             plt.title('Execution Time Comparison')
             plt.legend()
             plt.grid(True)
             plt.xscale('log')
             plt.yscale('log')
             plt.savefig(os.path.join(plots_dir, 'summary_cpu_time.png'))
             plt.close()
        except:
            print("Could not create summary timing plot.")
    
    # 3. Max Error Summary
    # We need to read all history files to get max error for each run
    summary_data = []
    
    # Walk through runs
    run_folders = [f.path for f in os.scandir(runs_dir) if f.is_dir()]
    
    for rf in run_folders:
        tag = os.path.basename(rf)
        hist_file = os.path.join(rf, 'tip_history.csv')
        
        if not os.path.exists(hist_file):
            continue
            
        try:
            df = pd.read_csv(hist_file)
            
            # Extract parameters from tag or lookup?
            # Tag format: ne{ne}_f{abs(force)}_inc{ninc} or "default" or "test_otis"
            # We can try to parse tag, or just rely on plotting per tag.
            # But summary plots need X-axis (nelems). 
            # We can try to join with timings.csv by TAG if tags match?
            # Yes, run_comparison.sh uses TAG in timings.csv
            
            # Compute max errors (ignore NaNs)
            # Use small epsilon for log plot safety
            eps = 1e-20
            
            max_err_hd = df['err_hd'].max() if 'err_hd' in df.columns else None
            max_err_fd = df['err_fd'].max() if 'err_fd' in df.columns else None
            max_err_cs = df['err_cs'].max() if 'err_cs' in df.columns else None
            max_err_ot = df['err_otis'].max() if 'err_otis' in df.columns else None
            
            # Append 
            summary_data.append({
                'tag': tag,
                'err_hd': max(max_err_hd, eps) if pd.notnull(max_err_hd) else None,
                'err_fd': max(max_err_fd, eps) if pd.notnull(max_err_fd) else None,
                'err_cs': max(max_err_cs, eps) if pd.notnull(max_err_cs) else None,
                'err_ot': max(max_err_ot, eps) if pd.notnull(max_err_ot) else None
            })
            
            # Individual Plot
            plt.figure(figsize=(10,6))
            plt.plot(df['frame'], df['err_hd'], label='HyperDual', marker='.')
            plt.plot(df['frame'], df['err_fd'], label='FiniteDiff', marker='.')
            plt.plot(df['frame'], df['err_cs'], label='ComplexStep', marker='.')
            plt.plot(df['frame'], df['err_otis'], label='Otis', marker='.')
            plt.yscale('log')
            plt.xlabel('Increment')
            plt.ylabel('Tip Error (Log)')
            plt.title(f'Error History: {tag}')
            plt.legend()
            plt.grid(True)
            plt.savefig(os.path.join(plots_dir, f'error_hist_{tag}.png'))
            plt.close()
            
        except Exception as e:
            print(f"Failed to process {tag}: {e}")

    # Merge with timings to get parameters
    if os.path.exists(timings_file) and summary_data:
        df_sum = pd.DataFrame(summary_data)
        df_time = pd.read_csv(timings_file)
        
        # Timings csv has multiple rows per tag (one per method).
        # We need parameters (nelems, force, nincs) from it. 
        # Drop duplicates on tag.
        df_params = df_time[['tag', 'nelems', 'force', 'nincs']].drop_duplicates()
        
        df_merged = pd.merge(df_sum, df_params, on='tag', how='inner')
        
        if not df_merged.empty:
            df_merged = df_merged.sort_values('nelems')
            
            # Plot Error vs Nelems
            plt.figure(figsize=(10,6))
            plt.plot(df_merged['nelems'], df_merged['err_hd'], 'o-', label='HyperDual')
            plt.plot(df_merged['nelems'], df_merged['err_fd'], 'x--', label='FiniteDiff')
            plt.plot(df_merged['nelems'], df_merged['err_cs'], 's-.', label='ComplexStep')
            plt.plot(df_merged['nelems'], df_merged['err_ot'], 'd:', label='Otis')
            
            plt.xscale('log')
            plt.yscale('log')
            plt.xlabel('Number of Elements')
            plt.ylabel('Max Tip Error')
            plt.title('Convergence: Max Error vs Mesh Size')
            plt.legend()
            plt.grid(True)
            plt.savefig(os.path.join(plots_dir, 'convergence_summary.png'))
            plt.close()
            
    print("Plots saved.")

if __name__ == "__main__":
    main()
