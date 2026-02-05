import os
import pandas as pd
import matplotlib.pyplot as plt
import re

def main():
    sim_dir = os.path.join(os.path.dirname(__file__), '../sim')
    timings_file = os.path.join(sim_dir, 'timings.csv')
    plot_file = os.path.join(sim_dir, 'plots', 'time_vs_elements.png')
    
    if not os.path.exists(os.path.dirname(plot_file)):
        os.makedirs(os.path.dirname(plot_file))
        
    if not os.path.exists(timings_file):
        print("timings.csv not found")
        return

    # Load data
    df = pd.read_csv(timings_file)
    
    # Filter for the linear sweep we just ran
    # tag format: lin_ne{N}_f10000_inc10
    mask = df['tag'].str.match(r'lin_ne\d+_f10000_inc10')
    df_filtered = df[mask].copy()
    
    if df_filtered.empty:
        print("No matching data found in timings.csv")
        return

    # --- Validation Step ---
    # Check each run directory to see if the analysis completed (reached last frame)
    # We assume 'runs_dir' structure
    runs_dir = os.path.join(sim_dir, 'runs')
    
    valid_runs = set() # Store (tag, method) tuples that are valid
    
    # Method name mapping timings -> tip_history column suffix
    method_map = {
        'HD': 'hd',
        'FD': 'fd',
        'CS': 'cs',
        'Exact': 'exact',
        'Otis': 'otis'
    }

    unique_tags = df_filtered['tag'].unique()
    print(f"Validating {len(unique_tags)} runs...")

    for tag in unique_tags:
        hist_file = os.path.join(runs_dir, tag, 'tip_history.csv')
        if not os.path.exists(hist_file):
            print(f"Warning: History file missing for {tag}")
            continue
            
        try:
            df_hist = pd.read_csv(hist_file)
            if df_hist.empty:
                continue
                
            # Check last row
            last_row = df_hist.iloc[-1]
            
            # Strict Check: All methods must be valid
            all_valid = True
            
            for t_method, h_suffix in method_map.items():
                col = f'u2_{h_suffix}'
                if col not in df_hist.columns:
                    all_valid = False
                    break
                val = last_row[col]
                if pd.isna(val):
                    all_valid = False
                    break
            
            if all_valid:
                for t_method in method_map.keys():
                    valid_runs.add((tag, t_method))

        except Exception as e:
            print(f"Error reading {hist_file}: {e}")

    # Filter df_filtered based on valid_runs
    # We can create a tuple column for filtering
    df_filtered['valid_check'] = list(zip(df_filtered['tag'], df_filtered['method']))
    
    initial_count = len(df_filtered)
    df_filtered = df_filtered[df_filtered['valid_check'].isin(valid_runs)]
    final_count = len(df_filtered)
    
    print(f"Filtered out {initial_count - final_count} incomplete runs.")
    print(f"Remaining: {final_count} data points.")

    # Convert columns to numeric
    df_filtered['nelems'] = pd.to_numeric(df_filtered['nelems'])
    df_filtered['cpu_time_seconds'] = pd.to_numeric(df_filtered['cpu_time_seconds'], errors='coerce')
    
    # Pivot or group by method
    methods = df_filtered['method'].unique()
    
    plt.figure(figsize=(10, 6))
    
    # Map method names for consistent legend
    label_map = {
        'HD': 'HyperDual (UEL)',
        'FD': 'Finite Difference',
        'CS': 'Complex Step',
        'Otis': 'Otis (HyperDual)',
        'Exact': 'Exact UEL'
    }
    
    # Plot style map
    style_map = {
        'HD': 'o-',
        'FD': 's-',
        'CS': '^-',
        'Otis': 'd-',
        'Exact': 'x-'
    }

    # Sort by elements
    df_filtered = df_filtered.sort_values('nelems')

    for method in methods:
        subset = df_filtered[df_filtered['method'] == method]
        # Remove duplicates if any (e.g. if run multiple times) - take the last one
        subset = subset.drop_duplicates(subset=['nelems'], keep='last')
        
        lbl = label_map.get(method, method)
        st = style_map.get(method, '.-')
        
        plt.plot(subset['nelems'], subset['cpu_time_seconds'], st, label=lbl, alpha=0.8, markersize=4)

    plt.xlabel('Number of Elements', fontsize=12)
    plt.ylabel('CPU Time (seconds)', fontsize=12)
    plt.title('Execution Time vs Number of Elements (Force = 10000)', fontsize=14)
    
    plt.yscale('linear')
    plt.xscale('linear')
    
    # Customize grid and ticks
    plt.grid(True, which="major", ls="-", alpha=0.8)
    plt.grid(True, which="minor", ls=":", alpha=0.5)
    
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.legend(fontsize=10)
    
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {plot_file}")

if __name__ == "__main__":
    main()
