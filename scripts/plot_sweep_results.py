import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import os
import sys
import glob

def main():
    sim_dir = os.path.join(os.path.dirname(__file__), '../sim')
    runs_dir = os.path.join(sim_dir, 'runs')
    plots_dir = os.path.join(sim_dir, 'plots')
    timings_file = os.path.join(sim_dir, 'timings.csv')
    
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)
        
    print(f"Generating comprehensive plots in {plots_dir}...")
    
    # --- 1. Load and Verify Data ---
    
    # Load Timings
    if not os.path.exists(timings_file):
        print("Error: timings.csv not found.")
        return
        
    df_time = pd.read_csv(timings_file)
    # timings.csv cols: tag,nelems,force,nincs,method,cpu_time_seconds
    
    # Load Errors from individual runs
    summary_data = []
    # Find all run directories (subfolders in runs_dir)
    if os.path.exists(runs_dir):
        run_folders = [f.path for f in os.scandir(runs_dir) if f.is_dir()]
    else:
        run_folders = []
    
    print(f"Processing {len(run_folders)} run directories...")
    
    for rf in run_folders:
        tag = os.path.basename(rf)
        hist_file = os.path.join(rf, 'tip_history.csv')
        
        if not os.path.exists(hist_file):
            continue
            
        try:
            df = pd.read_csv(hist_file)
            eps = 1e-20
            
            # Get max errors
            err_hd = df['err_hd'].max() if 'err_hd' in df.columns else None
            err_fd = df['err_fd'].max() if 'err_fd' in df.columns else None
            err_cs = df['err_cs'].max() if 'err_cs' in df.columns else None
            err_ot = df['err_otis'].max() if 'err_otis' in df.columns else None
            
            summary_data.append({
                'tag': tag,
                'err_hd': max(err_hd, eps) if pd.notnull(err_hd) else None,
                'err_fd': max(err_fd, eps) if pd.notnull(err_fd) else None,
                'err_cs': max(err_cs, eps) if pd.notnull(err_cs) else None,
                'err_ot': max(err_ot, eps) if pd.notnull(err_ot) else None
            })
        except Exception as e:
            # print(f"Skipping {tag}: {e}")
            pass

    if not summary_data:
        print("No valid simulation data found.")
        return

    df_err = pd.DataFrame(summary_data)
    
    # Merge Errors with Parameters (from timings - unique tags only)
    # We need to drop duplicates in timings for the parameter lookup
    df_params = df_time[['tag', 'nelems', 'force', 'nincs']].drop_duplicates()
    
    # Master DataFrame: tag, nelems, force, nincs, err_hd, err_fd...
    df_master = pd.merge(df_err, df_params, on='tag', how='inner')
    
    # We also need a Master Time DataFrame.
    # The timings df has 'method' column (long format). Error df is wide.
    # Let's keep Time analysis separate or pivot it.
    
    # --- 2. Plotting Functions ---
    
    def plot_metric(data, x_col, y_cols, y_label, title, filename, x_log=False, y_log=False):
        plt.figure(figsize=(10, 6))
        
        markers = {'err_hd': 'o', 'err_fd': 'x', 'err_cs': 's', 'err_ot': 'd',
                   'HD': 'o', 'FD': 'x', 'CS': 's', 'Otis': 'd', 'Exact': '+'}
        labels = {'err_hd': 'HyperDual', 'err_fd': 'FiniteDiff', 'err_cs': 'ComplexStep', 'err_ot': 'Otis',
                  'HD': 'HyperDual', 'FD': 'FiniteDiff', 'CS': 'ComplexStep', 'Otis': 'Otis', 'Exact': 'Exact'}
        colors = {'err_hd': 'blue', 'err_fd': 'red', 'err_cs': 'green', 'err_ot': 'magenta',
                  'HD': 'blue', 'FD': 'red', 'CS': 'green', 'Otis': 'magenta', 'Exact': 'orange'}

        # Sort by x
        data = data.sort_values(x_col)
        
        # Ensure numeric
        data[x_col] = pd.to_numeric(data[x_col])
        if x_col == 'force': data[x_col] = data[x_col].abs() # Plot against magnitude if negative
        
        plotted_something = False
        for col in y_cols:
            if col in data.columns and not data[col].isna().all():
                 lbl = labels.get(col, col)
                 mk = markers.get(col, '.')
                 clr = colors.get(col, None)
                 plt.plot(data[x_col], data[col], marker=mk, label=lbl, color=clr, linestyle='-', linewidth=1.5)
                 plotted_something = True
        
        if not plotted_something:
            plt.close()
            return

        if x_log: plt.xscale('log')
        if y_log: plt.yscale('log')
        
        plt.xlabel(x_col.replace('_', ' ').title() + (" (Abs)" if x_col=='force' else ""))
        plt.ylabel(y_label)
        plt.title(title)
        plt.grid(True, which="both", ls="-", alpha=0.4)
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, filename))
        plt.close()

    def plot_metric_plotly(data, x_col, y_cols, y_label, title, filename, x_log=False, y_log=False):
        fig = go.Figure()
        
        # Mapping for symbols and names
        markers = {'err_hd': 'circle', 'err_fd': 'x', 'err_cs': 'square', 'err_ot': 'diamond',
                   'HD': 'circle', 'FD': 'x', 'CS': 'square', 'Otis': 'diamond', 'Exact': 'cross'}
        labels = {'err_hd': 'HyperDual', 'err_fd': 'FiniteDiff', 'err_cs': 'ComplexStep', 'err_ot': 'Otis',
                  'HD': 'HyperDual', 'FD': 'FiniteDiff', 'CS': 'ComplexStep', 'Otis': 'Otis', 'Exact': 'Exact'}
        colors = {'err_hd': 'blue', 'err_fd': 'red', 'err_cs': 'green', 'err_ot': 'magenta',
                  'HD': 'blue', 'FD': 'red', 'CS': 'green', 'Otis': 'magenta', 'Exact': 'orange'}
        
        # Sort by x
        data = data.sort_values(x_col)
        data[x_col] = pd.to_numeric(data[x_col])
        if x_col == 'force': data[x_col] = data[x_col].abs()

        for col in y_cols:
            if col in data.columns and not data[col].isna().all():
                lbl = labels.get(col, col)
                mk = markers.get(col, 'circle')
                clr = colors.get(col, None)
                
                fig.add_trace(go.Scatter(
                    x=data[x_col], 
                    y=data[col],
                    mode='lines+markers',
                    name=lbl,
                    marker=dict(symbol=mk, color=clr, size=8),
                    line=dict(color=clr, width=2)
                ))

        if x_log: fig.update_xaxes(type="log")
        if y_log: fig.update_yaxes(type="log")
        
        fig.update_layout(
            title=title,
            xaxis_title=x_col.replace('_', ' ').title() + (" (Abs)" if x_col=='force' else ""),
            yaxis_title=y_label,
            template="plotly_white",
            hovermode="x unified"
        )
        
        out_path = os.path.join(plots_dir, filename.replace('.png', '.html'))
        fig.write_html(out_path)

    # --- 3. Generate Dependency Plots ---

    # Define "Central" Parameters (Median or standard values) to slice through
    # Available values
    uniq_ne = sorted(df_master['nelems'].unique())
    uniq_f  = sorted(df_master['force'].unique())
    uniq_ni = sorted(df_master['nincs'].unique())
    
    if not uniq_ne or not uniq_f or not uniq_ni:
        print("Insufficient parameter data.")
        return
    
    # Pick a "middle" slice for fixed params
    # Prefer standardized values if available: 16, -1000, 50 ?
    # From sweep: [4, 16, 64], [-100, -1000, -5000], [10, 50, 100]
    # Median is perfect.
    fixed_ne = uniq_ne[len(uniq_ne)//2] 
    fixed_f  = uniq_f[len(uniq_f)//2]   
    fixed_ni = uniq_ni[len(uniq_ni)//2] 
    
    print(f"Base Slice Params: Ne={fixed_ne}, F={fixed_f}, Inc={fixed_ni}")

    # A. VARYING FORCE (Fixed Ne, Fixed Nincs)
    subset_f = df_master[(df_master['nelems'] == fixed_ne) & (df_master['nincs'] == fixed_ni)]
    if not subset_f.empty:
        # Error vs Force
        plot_metric(subset_f, 'force', ['err_hd', 'err_fd', 'err_cs', 'err_ot'], 
                    'Max Tip Error', f'Error vs Force (Ne={fixed_ne}, Inc={fixed_ni})', 'error_vs_force.png', x_log=True, y_log=True)
        plot_metric_plotly(subset_f, 'force', ['err_hd', 'err_fd', 'err_cs', 'err_ot'], 
                    'Max Tip Error', f'Error vs Force (Ne={fixed_ne}, Inc={fixed_ni})', 'error_vs_force.png', x_log=True, y_log=True)

        # Time vs Force
        subset_f_tags = subset_f['tag'].unique()
        subset_time = df_time[df_time['tag'].isin(subset_f_tags)]
        pivot_time = subset_time.pivot_table(index='force', columns='method', values='cpu_time_seconds').reset_index()
        plot_metric(pivot_time, 'force', ['HD', 'FD', 'CS', 'Otis', 'Exact'], 
                    'CPU Time (s)', f'Time vs Force (Ne={fixed_ne}, Inc={fixed_ni})', 'time_vs_force.png', x_log=True)
        plot_metric_plotly(pivot_time, 'force', ['HD', 'FD', 'CS', 'Otis', 'Exact'], 
                    'CPU Time (s)', f'Time vs Force (Ne={fixed_ne}, Inc={fixed_ni})', 'time_vs_force.png', x_log=True)

    # B. VARYING NELEMS (Fixed Force, Fixed Nincs)
    subset_n = df_master[(df_master['force'] == fixed_f) & (df_master['nincs'] == fixed_ni)]
    if not subset_n.empty:
        # Error vs Nelems
        plot_metric(subset_n, 'nelems', ['err_hd', 'err_fd', 'err_cs', 'err_ot'], 
                    'Max Tip Error', f'Error vs Mesh Size (F={fixed_f}, Inc={fixed_ni})', 'error_vs_nelems.png', x_log=True, y_log=True)
        plot_metric_plotly(subset_n, 'nelems', ['err_hd', 'err_fd', 'err_cs', 'err_ot'], 
                    'Max Tip Error', f'Error vs Mesh Size (F={fixed_f}, Inc={fixed_ni})', 'error_vs_nelems.png', x_log=True, y_log=True)
        
        subset_n_tags = subset_n['tag'].unique()
        subset_time = df_time[df_time['tag'].isin(subset_n_tags)]
        pivot_time = subset_time.pivot_table(index='nelems', columns='method', values='cpu_time_seconds').reset_index()
        plot_metric(pivot_time, 'nelems', ['HD', 'FD', 'CS', 'Otis', 'Exact'], 
                    'CPU Time (s)', f'Time vs Mesh Size (F={fixed_f}, Inc={fixed_ni})', 'time_vs_nelems.png', x_log=True)
        plot_metric_plotly(pivot_time, 'nelems', ['HD', 'FD', 'CS', 'Otis', 'Exact'], 
                    'CPU Time (s)', f'Time vs Mesh Size (F={fixed_f}, Inc={fixed_ni})', 'time_vs_nelems.png', x_log=True)

    # C. VARYING NINCS (Fixed Force, Fixed Nelems)
    subset_i = df_master[(df_master['force'] == fixed_f) & (df_master['nelems'] == fixed_ne)]
    if not subset_i.empty:
        # Error vs Nincs
        plot_metric(subset_i, 'nincs', ['err_hd', 'err_fd', 'err_cs', 'err_ot'], 
                    'Max Tip Error', f'Error vs Increments (F={fixed_f}, Ne={fixed_ne})', 'error_vs_nincs.png', y_log=True)
        plot_metric_plotly(subset_i, 'nincs', ['err_hd', 'err_fd', 'err_cs', 'err_ot'], 
                    'Max Tip Error', f'Error vs Increments (F={fixed_f}, Ne={fixed_ne})', 'error_vs_nincs.png', y_log=True)
        
        subset_i_tags = subset_i['tag'].unique()
        subset_time = df_time[df_time['tag'].isin(subset_i_tags)]
        pivot_time = subset_time.pivot_table(index='nincs', columns='method', values='cpu_time_seconds').reset_index()
        plot_metric(pivot_time, 'nincs', ['HD', 'FD', 'CS', 'Otis', 'Exact'], 
                    'CPU Time (s)', f'Time vs Increments (F={fixed_f}, Ne={fixed_ne})', 'time_vs_nincs.png')
        plot_metric_plotly(pivot_time, 'nincs', ['HD', 'FD', 'CS', 'Otis', 'Exact'], 
                    'CPU Time (s)', f'Time vs Increments (F={fixed_f}, Ne={fixed_ne})', 'time_vs_nincs.png')

    # --- 4. Average Statistics Plots ---
    
    print("Generating Averaged Statistics Plots...")
    
    # 4.1 Normalized Runtime (Time per Element per Increment)
    # Use full df_time dataset
    df_time['norm_time'] = df_time['cpu_time_seconds'] / (df_time['nelems'] * df_time['nincs'])
    
    # Group by method and average
    avg_time = df_time.groupby('method')['norm_time'].mean().reset_index()
    std_time = df_time.groupby('method')['norm_time'].std().reset_index()
    
    # Plot Bar Chart
    plt.figure(figsize=(10, 6))
    methods = avg_time['method']
    y_pos = range(len(methods))
    plt.bar(y_pos, avg_time['norm_time'] * 1000, yerr=std_time['norm_time'] * 1000, align='center', alpha=0.7, capsize=10)
    plt.xticks(y_pos, methods)
    plt.ylabel('Avg Time per Element-Increment (ms)')
    plt.title('Normalized Computational Cost')
    plt.grid(axis='y', alpha=0.5)
    plt.savefig(os.path.join(plots_dir, 'avg_normalized_runtime.png'))
    plt.close()
    
    # Plotly Bar Chart
    fig_bar = px.bar(
        avg_time, x='method', y=avg_time['norm_time'] * 1000,
        error_y=std_time['norm_time'] * 1000,
        title='Normalized Computational Cost',
        labels={'y': 'Avg Time per Element-Increment (ms)', 'method': 'Method'},
        template="plotly_white"
    )
    fig_bar.write_html(os.path.join(plots_dir, 'avg_normalized_runtime.html'))
    
    # 4.2 Average Error Distribution
    # Use df_err (Max error per run)
    # Melt to long format
    df_err_long = df_err.melt(id_vars=['tag'], value_vars=['err_hd', 'err_fd', 'err_cs', 'err_ot'], 
                              var_name='method_key', value_name='error')
    
    # Map keys to display names
    name_map = {'err_hd': 'HyperDual', 'err_fd': 'FiniteDiff', 'err_cs': 'ComplexStep', 'err_ot': 'Otis'}
    df_err_long['method'] = df_err_long['method_key'].map(name_map)
    
    # Calculate Mean Error (Geometric Mean might be better for log data, but user asked for average)
    # Let's do Box Plot to show distribution
    
    plt.figure(figsize=(10, 6))
    
    # Prepare data list for boxplot
    plot_data = []
    labels = []
    for m in ['HyperDual', 'FiniteDiff', 'ComplexStep', 'Otis']:
        subset = df_err_long[df_err_long['method'] == m]['error'].dropna()
        if len(subset) > 0:
            plot_data.append(subset.values)
            labels.append(m)
            
    plt.boxplot(plot_data, labels=labels)
    plt.yscale('log')
    plt.ylabel('Max Tip Error (Log Scale)')
    plt.title('Error Distribution across all Sweeps')
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.savefig(os.path.join(plots_dir, 'avg_error_distribution.png'))
    plt.close()

    # Plotly Box Plot
    fig_box = px.box(
        df_err_long, x='method', y='error',
        title='Error Distribution across all Sweeps',
        log_y=True,
        labels={'error': 'Max Tip Error (Log Scale)', 'method': 'Method'},
        template="plotly_white"
    )
    fig_box.write_html(os.path.join(plots_dir, 'avg_error_distribution.html'))

    print("Plots generated successfully.")

if __name__ == "__main__":
    main()
