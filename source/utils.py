from pylr2 import regress2 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
"""from lmfit.models import ExponentialModel
"""
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
import anndata

"""def plot_exponential(mean_df, axes, color, treatment_col, x_col):
    # fit negative exponential using lmfit
    model = ExponentialModel()
    params = model.guess(mean_df[treatment_col], x=mean_df[x_col])
    result = model.fit(mean_df[treatment_col], params, x=mean_df[x_col])
    # plot the fit at 1,000 values between min and max
    minn = mean_df[x_col].min()
    maxx = mean_df[x_col].max()
    x_vals = np.linspace(
        minn, maxx, 1000
        )
    y_vals = result.eval(x=x_vals)
    axes.plot(x_vals, y_vals, color=color)
    # plot confidence interval
    ci = result.eval_uncertainty(x=x_vals)
    axes.fill_between(x_vals, y_vals - ci, y_vals + ci, alpha=0.3, color=color)"""
    
        

def _fxn_to_apply(name:str, values:np.array):
    if name=='inverse':
        return 1/values
    elif name=='log':
        return np.log10(values)
    elif name==None:
        return values

def plot_major_axis_regr(
    to_use_df, 
    x_str, 
    y_str, 
    axes, 
    fxn_to_apply = None
    ):
    # scale x vals optionally
    x = _fxn_to_apply(fxn_to_apply, to_use_df[x_str].values)
    y = to_use_df[y_str].values
        
    results = regress2(
        _x = x , _y = y, _need_intercept = True
    )
    
    slope, intercept, r_value, std_err, std_intercept, predicted_vals = results['slope'], results['intercept'], results['r'], results['std_slope'], results['std_intercept'], results['predict']
    print(f"R: {r_value},  slope: {slope}, intercept: {intercept}, std_err: {std_err}, std_intercept: {std_intercept}")
    
    # plot the line
    x_vals = np.linspace(x.min(), x.max(), 1000)
    y_vals = slope * x_vals + intercept
    axes.plot(_fxn_to_apply(fxn_to_apply, x_vals), y_vals, color = 'black')
    
    # plot confidence intervals
    # check if 0 is in the range of x values
    """if 0 in x_vals:
        # until 0, use upper bound 1 and lowest bound2
        x_vals1 = np.linspace(x.min(), 0, 1000)
        
        lower_bound1 = (slope - std_err) * x_vals1 + (intercept - std_intercept)
        upper_bound1 = (slope - std_err) * x_vals1 + (intercept + std_intercept)
        lower_bound2 = (slope + std_err) * x_vals1 + (intercept - std_intercept)
        upper_bound2 = (slope + std_err) * x_vals1 + (intercept + std_intercept)
        axes.fill_between(_fxn_to_apply(fxn_to_apply, x_vals1), lower_bound2, upper_bound1, alpha=0.3, color='black',edgecolor = 'none')
        # after 0, use upper bound 2 and lower bound 1
        x_vals2 = np.linspace(0, x.max(), 1000)
        
        lower_bound1 = (slope - std_err) * x_vals2 + (intercept - std_intercept)
        upper_bound1 = (slope - std_err) * x_vals2 + (intercept + std_intercept)
        lower_bound2 = (slope + std_err) * x_vals2 + (intercept - std_intercept)
        upper_bound2 = (slope + std_err) * x_vals2 + (intercept + std_intercept)
        axes.fill_between(_fxn_to_apply(fxn_to_apply, x_vals2), lower_bound1, upper_bound2, alpha=0.3, color='black',  edgecolor = 'none')
    elif 0 < x_vals[0]:
        lower_bound = (slope - std_err) * x_vals + (intercept - std_intercept)
        upper_bound = (slope + std_err) * x_vals + (intercept + std_intercept)
        axes.fill_between(_fxn_to_apply(fxn_to_apply, x_vals), lower_bound, upper_bound, alpha=0.3, color='black', edgecolor = 'none')
    else:
        lower_bound = (slope - std_err) * x_vals + (intercept - std_intercept)
        upper_bound = (slope + std_err) * x_vals + (intercept + std_intercept)
        axes.fill_between(_fxn_to_apply(fxn_to_apply, x_vals), lower_bound, upper_bound, alpha=0.3, color='black', edgecolor = 'none')"""
    
def read_dream_files(
    use_new: bool = False
    ) -> None:
    """Read in DREAM files from /data/bujarrabal_dueso
    ### Parameters:
    use_new : bool
        Wether to use the new dream file
    ### Returns:
    dream_regulated_genes : pd.DataFrame
        DataFrame containing the DREAM regulated genes
    """
    if use_new:
        # index 'Gene name'
        # columns: 'ensembl_id', 'hgnc_id', 'uniprotkb
        dream_regulated_genes = pd.read_parquet(
            '/cellar/users/zkoch/dream/data/bujarrabal_dueso/570_DREAM_revisited_list_processed.parquet'
            )
        # replace spaces with underscores and make lowercase
        dream_regulated_genes.columns = [
            x.lower().replace(" ", "_") for x in dream_regulated_genes.columns
            ]
        dream_regulated_genes.set_index('gene_stable_id', inplace = True, drop = False)
    else:
        # index 'Gene name'
        # columns: 'gene_stable_id', 'dna_repair_genes', 'harmine_dna_repair_up', 'indy_dna_repair_up'
        dream_regulated_genes = pd.read_csv(
            "/cellar/users/zkoch/dream/data/bujarrabal_dueso/tableS12_dream_promoter_binding.csv", index_col=None
            )
        # replace spaces with underscores and make lowercase
        dream_regulated_genes.columns = [
            x.lower().replace(" ", "_") for x in dream_regulated_genes.columns
            ]
        dream_regulated_genes.set_index('gene_stable_id', inplace=True)
    return dream_regulated_genes


def read_subset_h5ad(
    file_path: str,
    start_row: int,
    end_row: int
    ):
    """
    Reads a subset of rows from an h5ad file without loading the entire file into memory.

    Parameters:
    - file_path (str): Path to the h5ad file.
    - start_row (int): The starting row index (inclusive).
    - end_row (int): The ending row index (exclusive).

    Returns:
    - anndata.AnnData: An AnnData object containing the subset of data.
    """
    # Open the h5ad file in backed mode to avoid loading the entire file into memory
    adata = anndata.read_h5ad(file_path, backed='r')
    # Slice the data to get the desired rows and create a copy in memory
    adata_subset = adata[start_row:end_row].to_memory()
    # Ensure the file is closed to free up resources
    adata.file.close()
    return adata_subset


def plot_dbs_signatures(file_path, figsize=(15, 8), percentage=False):
    """
    Plot DBS (Doublet Base Substitution) signature patterns similar to SigProfiler's plotDBS.
    
    Parameters:
    - file_path: Path to tab-delimited file with MutationType column and signature columns
    - figsize: Figure size tuple
    - percentage: Whether to plot as percentages
    
    Returns:
    - fig: matplotlib figure object
    - axes: matplotlib axes object(s)
    """
    # Read data
    df = pd.read_csv(file_path, sep='\t', index_col=0)
    
    # DBS78 color scheme (same as SigProfiler)
    colors = [
        [3/256, 189/256, 239/256],    # Light blue (TT>NN)
        [162/256, 207/256, 99/256],   # Light green (AC>NN)  
        [255/256, 153/256, 153/256],  # Light red (CT>NN)
        [3/256, 102/256, 204/256],    # Blue (AT>NN)
        [255/256, 178/256, 102/256],  # Light orange (TG>NN)
        [228/256, 41/256, 38/256],    # Red (CC>NN)
        [76/256, 1/256, 153/256],     # Purple (CG>NN)
        [255/256, 128/256, 1/256],    # Orange (TC>NN)
        [1/256, 102/256, 1/256],      # Dark green (GC>NN)
        [204/256, 153/256, 255/256],  # Light purple (TA>NN)
    ]
    
    # DBS context mapping
    context_colors = {}
    contexts = ['TT>', 'AC>', 'CT>', 'AT>', 'TG>', 'CC>', 'CG>', 'TC>', 'GC>', 'TA>']
    for i, context in enumerate(contexts):
        context_colors[context] = colors[i]
    
    # Get number of signatures
    signatures = df.columns.tolist()
    n_sigs = len(signatures)
    
    # Create subplot grid
    cols = min(4, n_sigs)
    rows = int(np.ceil(n_sigs / cols))
    
    fig, axes = plt.subplots(rows, cols, figsize=(figsize[0]*cols/2, figsize[1]*rows/2))
    if n_sigs == 1:
        axes = [axes]
    elif rows == 1:
        axes = axes.reshape(1, -1)
    axes = axes.flatten()
    
    # Plot each signature
    for sig_idx, sig in enumerate(signatures):
        ax = axes[sig_idx]
        
        # Get values
        values = df[sig].values
        if percentage:
            values = values / values.sum() * 100
            
        # Assign colors based on context
        bar_colors = []
        for mut_type in df.index:
            context = mut_type[:3]
            bar_colors.append(context_colors.get(context, [0.5, 0.5, 0.5]))
        
        # Create bar plot
        x_pos = range(len(values))
        bars = ax.bar(x_pos, values, color=bar_colors, width=0.8, edgecolor='white', linewidth=0.5)
        
        # Formatting
        ax.set_title(f'{sig}', fontsize=12, fontweight='bold')
        ax.set_xlim(-0.5, len(values)-0.5)
        ax.grid(True, axis='y', alpha=0.3)
        ax.set_axisbelow(True)
        
        # X-axis labels (rotated)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(df.index, rotation=90, ha='center', fontsize=8)
        
        # Y-axis
        if percentage:
            ax.set_ylabel('Percentage', fontsize=10)
        else:
            ax.set_ylabel('Mutations', fontsize=10)
    
    # Remove empty subplots
    for i in range(n_sigs, len(axes)):
        fig.delaxes(axes[i])
    
    plt.tight_layout()
    
    return fig, axes


def plot_sbs_signatures(file_path, figsize=(20, 8), percentage=False):
    """
    Plot SBS (Single Base Substitution) signature patterns similar to SigProfiler's plotSBS.
    
    Parameters:
    - file_path: Path to tab-delimited file with MutationType column and signature columns
    - figsize: Figure size tuple
    - percentage: Whether to plot as percentages
    
    Returns:
    - fig: matplotlib figure object
    - axes: matplotlib axes object(s)
    """
    # Read data
    df = pd.read_csv(file_path, sep='\t', index_col=0)
    
    # SBS96 color scheme (same as SigProfiler)
    colors = {
        'C>A': [3/256, 189/256, 239/256],     # Light blue
        'C>G': [162/256, 207/256, 99/256],    # Light green
        'C>T': [255/256, 153/256, 153/256],   # Light red
        'T>A': [228/256, 41/256, 38/256],     # Red  
        'T>C': [255/256, 178/256, 102/256],   # Light orange
        'T>G': [76/256, 1/256, 153/256],      # Purple
    }
    
    # Get number of signatures
    signatures = df.columns.tolist()
    n_sigs = len(signatures)
    
    # Create subplot grid
    cols = min(3, n_sigs)
    rows = int(np.ceil(n_sigs / cols))
    
    fig, axes = plt.subplots(rows, cols, figsize=(figsize[0]*cols/3, figsize[1]*rows/3))
    if n_sigs == 1:
        axes = [axes]
    elif rows == 1:
        axes = axes.reshape(1, -1)
    axes = axes.flatten()
    
    # Plot each signature
    for sig_idx, sig in enumerate(signatures):
        ax = axes[sig_idx]
        
        # Get values
        values = df[sig].values
        if percentage:
            values = values / values.sum() * 100
            
        # Extract mutation types from index (e.g., A[C>A]A -> C>A)
        mut_types = []
        bar_colors = []
        for mut_context in df.index:
            # Extract mutation type from context like A[C>A]A
            mut_type = mut_context.split('[')[1].split(']')[0]
            mut_types.append(mut_type)
            bar_colors.append(colors.get(mut_type, [0.5, 0.5, 0.5]))
        
        # Create bar plot
        x_pos = range(len(values))
        bars = ax.bar(x_pos, values, color=bar_colors, width=0.8, edgecolor='white', linewidth=0.3)
        
        # Add mutation type labels at top
        unique_mut_types = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
        type_positions = []
        
        # Calculate positions for mutation type labels
        current_pos = 0
        for mut_type in unique_mut_types:
            count = sum(1 for mt in mut_types if mt == mut_type)
            if count > 0:
                type_positions.append((current_pos + count/2 - 0.5, mut_type))
                current_pos += count
        
        # Add mutation type labels
        for pos, label in type_positions:
            ax.text(pos, ax.get_ylim()[1] * 1.05, label, ha='center', va='bottom', 
                   fontweight='bold', fontsize=10)
        
        # Formatting
        ax.set_title(f'{sig}', fontsize=14, fontweight='bold')
        ax.set_xlim(-0.5, len(values)-0.5)
        ax.grid(True, axis='y', alpha=0.3)
        ax.set_axisbelow(True)
        
        # X-axis labels (trinucleotide contexts)
        ax.set_xticks(x_pos)
        contexts = [idx.replace('[', '').replace(']', '') for idx in df.index]
        ax.set_xticklabels(contexts, rotation=90, ha='center', fontsize=6)
        
        # Y-axis
        if percentage:
            ax.set_ylabel('Percentage', fontsize=10)
        else:
            ax.set_ylabel('Mutations', fontsize=10)
    
    # Remove empty subplots
    for i in range(n_sigs, len(axes)):
        fig.delaxes(axes[i])
    
    plt.tight_layout()
    
    return fig, axes


def plot_id_signatures(file_path, figsize=(15, 8), percentage=False):
    """
    Plot ID (Insertion/Deletion) signature patterns similar to SigProfiler's plotID.
    
    Parameters:
    - file_path: Path to tab-delimited file with MutationType column and signature columns
    - figsize: Figure size tuple  
    - percentage: Whether to plot as percentages
    
    Returns:
    - fig: matplotlib figure object
    - axes: matplotlib axes object(s)
    """
    # Read data
    df = pd.read_csv(file_path, sep='\t', index_col=0)
    
    # ID color scheme (similar to SigProfiler)
    colors = {
        'C': [253/256, 192/256, 134/256],     # Light orange (C homopolymers)
        'T': [255/256, 127/256, 14/256],      # Orange (T homopolymers)
        '2': [227/256, 26/256, 28/256],       # Red (2bp repeats)
        '3': [166/256, 86/256, 40/256],       # Brown (3bp repeats)
        '4': [106/256, 61/256, 154/256],      # Purple (4bp repeats)
        '5': [31/256, 120/256, 180/256],      # Blue (5bp+ repeats)
        'M': [178/256, 223/256, 138/256],     # Light green (Microhomology)
    }
    
    # Get number of signatures
    signatures = df.columns.tolist()
    n_sigs = len(signatures)
    
    # Create subplot grid
    cols = min(3, n_sigs)
    rows = int(np.ceil(n_sigs / cols))
    
    fig, axes = plt.subplots(rows, cols, figsize=(figsize[0]*cols/3, figsize[1]*rows/3))
    if n_sigs == 1:
        axes = [axes]
    elif rows == 1:
        axes = axes.reshape(1, -1)
    axes = axes.flatten()
    
    # Plot each signature
    for sig_idx, sig in enumerate(signatures):
        ax = axes[sig_idx]
        
        # Get values
        values = df[sig].values
        if percentage:
            values = values / values.sum() * 100
            
        # Assign colors based on mutation type
        bar_colors = []
        for mut_type in df.index:
            # Simple heuristic for ID classification
            if 'C:' in mut_type:
                color_key = 'C'
            elif 'T:' in mut_type:
                color_key = 'T'
            elif '2:' in mut_type:
                color_key = '2'
            elif '3:' in mut_type:
                color_key = '3'
            elif '4:' in mut_type:
                color_key = '4'
            elif '5:' in mut_type or '6:' in mut_type:
                color_key = '5'
            elif 'M:' in mut_type:
                color_key = 'M'
            else:
                color_key = 'M'  # Default
            
            bar_colors.append(colors.get(color_key, [0.5, 0.5, 0.5]))
        
        # Create bar plot
        x_pos = range(len(values))
        bars = ax.bar(x_pos, values, color=bar_colors, width=0.8, edgecolor='white', linewidth=0.5)
        
        # Formatting
        ax.set_title(f'{sig}', fontsize=12, fontweight='bold')
        ax.set_xlim(-0.5, len(values)-0.5)
        ax.grid(True, axis='y', alpha=0.3)
        ax.set_axisbelow(True)
        
        # X-axis labels
        ax.set_xticks(x_pos)
        ax.set_xticklabels(df.index, rotation=90, ha='center', fontsize=8)
        
        # Y-axis
        if percentage:
            ax.set_ylabel('Percentage', fontsize=10)
        else:
            ax.set_ylabel('Mutations', fontsize=10)
    
    # Remove empty subplots
    for i in range(n_sigs, len(axes)):
        fig.delaxes(axes[i])
    
    plt.tight_layout()
    
    return fig, axes