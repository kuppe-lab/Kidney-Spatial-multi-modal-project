import warnings
warnings.filterwarnings("ignore")
import session_info
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from scipy.stats import gaussian_kde
from matplotlib.colors import LinearSegmentedColormap
from sklearn.preprocessing import MinMaxScaler
from scipy import sparse
from scipy.stats import norm
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import anndata as ad
import scanpy as sc
import plotly.io as pio
pio.renderers.default = 'iframe'

def select_slide(adata, 
                 s: str = None,
                 s_col: str = 'sample'):
    r""" This function selects the data for one slide from the spatial anndata object.

    :param adata: Anndata object with multiple spatial experiments
    :param s: name of selected experiment
    :param s_col: column in adata.obs listing experiment name for each location
    """
    if s is None:
        return "The sample id should be provided!"
    else:
        slide = adata[adata.obs[s_col].isin([s]), :]
        s_keys = list(slide.uns['spatial'].keys())
        s_spatial = np.array(s_keys)[[s in k for k in s_keys]][0]
    
        slide.uns['spatial'] = {s_spatial: slide.uns['spatial'][s_spatial]}
    
        return slide

def make_proportions_table(adata,
                           sample_id: str = 'Sample',
                           celltype_id: str = 'Clusters2',
                           label_id: str = 'label'):
    r""" This function creates cell types proportions for each sample.

    :param adata: Anndata object with multiple spatial experiments
    :param sample_id: name of sample id column used in anndata object obs
    :param celltype_id: name of celltype id column used in anndata object obs
    :param label_id: name of sample id column used in anndata object obs
    """
    
    df = adata.obs.groupby([sample_id, celltype_id]).count()[label_id].to_frame()
    df = df.reset_index() 
    df_pivot = df.pivot_table(index = sample_id, columns = celltype_id, values = label_id)
    df_pivot = df_pivot.div(adata.obs.groupby([sample_id]).count()[label_id].astype(float), axis=0)
    temp = adata.obs[[sample_id, label_id]].drop_duplicates()
    unique_map = dict(zip(temp[sample_id], temp[label_id]))
    df_pivot[label_id] = [unique_map[s] for s in df_pivot.index]
    return df_pivot

def plot_cell_types_distributions(adata,
                                  cell_types: list = None,
                                  figsize: tuple = (15, 7),
                                  my_pal: tuple = None,
                                  sample_id: str = 'Sample',
                                  celltype_id: str = 'Clusters2',
                                  label_id: str = 'label',
                                  save_path: str = None,
                                  font_size: int = 24
                                 ):
    r""" This function plots cell types proportions distribution.

    :param adata: Anndata object with multiple spatial experiments
    :param cell_types: list of cell types to show cell types proportions distribution
    :param figsize: the figure size as tuple
    :param my_pal: a tuple to determine the color of each label group
    :param sample_id: name of sample id column used in anndata object obs
    :param celltype_id: name of celltype id column used in anndata object obs
    :param label_id: name of sample id column used in anndata object obs
    :param save_path: string to specify the directory to save pdf plot 
    """
    
    proportions = make_proportions_table(adata, sample_id, celltype_id, label_id)

    groups = adata.obs[label_id].unique()
    group1 = groups[0]
    group2 = groups[1]
    group1_proportions = proportions[proportions[label_id] == group1]
    group2_proportions = proportions[proportions[label_id] == group2]
    
    stats_group12 = []
    scores_group12 = []
    # using Welch’s t-test to detect the statistically significant difference
    ## between group 1 and group 2
    for cell_type in cell_types:
        #### Cohen’s D Score
        def get_effect_size_and_p(x, y):
            # 1. Calculate Cohen's d (Your existing logic)
            mean_diff = np.mean(y) - np.mean(x)
            pooled_std = np.sqrt((np.std(x, ddof=1)**2 + np.std(y, ddof=1)**2) / 2)
            d = mean_diff / pooled_std
            
            # 2. Calculate p-value using an independent t-test
            # We use stats.ttest_ind which assumes equal variance by default
            t_stat, p_value = stats.ttest_ind(y, x, equal_var=True)
            
            return d, p_value

        d_stat, p_val = get_effect_size_and_p(group1_proportions[cell_type],
                          group2_proportions[cell_type])
        scores_group12.append(d_stat)
        stats_group12.append(p_val)
        
        
    # adjust the p-values
    adjstats_group12 = multipletests(stats_group12, method='fdr_bh')[1]

    stats_bar_group12 = pd.DataFrame()
    stats_bar_group12['cell_type'] = cell_types
    stats_bar_group12['adjPval'] = adjstats_group12
    stats_bar_group12['Pval'] = stats_group12
    stats_bar_group12['score'] = scores_group12
    
    # sort the data based on p-value and statistic score
    stats_bar_group12 = stats_bar_group12.sort_values(by = ['score'],
                                                      ascending = [False])
    
    col_names = list(stats_bar_group12['cell_type'])
    col_names.append(label_id)
    df_long = pd.melt(proportions.loc[:,col_names], label_id)
    fig = plt.figure(figsize = figsize)
    ax = sns.boxplot(x = celltype_id, hue = label_id, y = "value", data = df_long, palette = my_pal)
    plt.ylabel("Cell frequency", fontsize = font_size)
    plt.xlabel("")
    plt.xticks(fontsize = font_size, rotation = 45, ha = 'right', rotation_mode = 'anchor')
    plt.yticks(fontsize = font_size - 4)
    plt.legend(fontsize = font_size)
    fig.tight_layout()
    if save_path is not None:
        if Path(save_path).is_dir():
            plt.savefig(save_path + "/cell_proportions_distribution.pdf")
        else:
            print("save path does not exist!")
    plt.show()

def plot_together(adata = None,
                  cell_types: list = None,
                  file_path: str = None,
                  figsize: tuple = (7, 6),
                  my_pal: tuple = None,
                  font_size: int = 14,
                  group1: str = 'Young',
                  group2: str = 'Old',
                  pval_thr: float = 0.05,
                  sample_id: str = 'Sample',
                  celltype_id: str = 'Clusters2',
                  label_id: str = 'label'
                                          
                 ):
    r""" This function plots cell types proportions distribution along the ranking score showing significant differences between two groups.

    :param adata: Anndata object with multiple spatial experiments
    :param cell_types: list of cell types to show cell types proportions distribution
    :param figsize: the figure size as tuple
    :param my_pal: a tuple to determine the color of each label group
    :param sample_id: name of sample id column used in anndata object obs
    :param celltype_id: name of celltype id column used in anndata object obs
    :param label_id: name of sample id column used in anndata object obs
    :param save_path: string to specify the directory to save pdf plot 
    """
    
    proportions = make_proportions_table(adata, sample_id, celltype_id, label_id)
    
    col_names = list(cell_types)
    col_names.append(label_id)
    df_long = pd.melt(proportions.loc[:,col_names], label_id)
    
    group1_proportions = proportions[proportions[label_id] == group1]
    group2_proportions = proportions[proportions[label_id] == group2]
    
    stats_group12 = []
    scores_group12 = []
    # using Welch’s t-test to detect the statistically significant difference
    ## between group 1 and group 2
    for cell_type in cell_types:
        #### Cohen’s D Score
        def get_effect_size_and_p(x, y):
            # 1. Calculate Cohen's d (Your existing logic)
            mean_diff = np.mean(y) - np.mean(x)
            pooled_std = np.sqrt((np.std(x, ddof=1)**2 + np.std(y, ddof=1)**2) / 2)
            d = mean_diff / pooled_std
            
            # 2. Calculate p-value using an independent t-test
            # We use stats.ttest_ind which assumes equal variance by default
            t_stat, p_value = stats.ttest_ind(y, x, equal_var=True)
            
            return d, p_value

        d_stat, p_val = get_effect_size_and_p(group1_proportions[cell_type],
                          group2_proportions[cell_type])
        scores_group12.append(d_stat)
        stats_group12.append(p_val)
        
        
    # adjust the p-values
    adjstats_group12 = multipletests(stats_group12, method='fdr_bh')[1]

    stats_bar_group12 = pd.DataFrame()
    stats_bar_group12['cell_type'] = cell_types
    stats_bar_group12['adjPval'] = adjstats_group12
    stats_bar_group12['Pval'] = stats_group12
    stats_bar_group12['score'] = scores_group12
    
    # sort the data based on p-value and statistic score
    stats_bar_group12 = stats_bar_group12.sort_values(by = ['score'],
                                                      ascending = [False])

    # Logic to define significance stars
    def get_star(p):
        if p < 0.001: return '***'
        elif p < 0.01: return '**'
        elif p < 0.05: return '*'
        else: return ''

    # Apply star logic to the stats dataframe
    stats_bar_group12['stars'] = stats_bar_group12['adjPval'].apply(get_star)

    order = stats_bar_group12['cell_type'].tolist()
    
    data = stats_bar_group12
    x = 'score'
    y = 'cell_type'
    
    p_d = data[x][data[x] >=0]
    p_normalized_values = (p_d - p_d.min()) / (p_d.max() - p_d.min())
    n_d = data[x][data[x] < 0]
    n_normalized_values = (n_d - n_d.min()) / (n_d.max() - n_d.min())

    colors = []
    cmap = mpl.colors.LinearSegmentedColormap.from_list('my_list', ["#f6d6c2", "#d47264", "#ae282c"])
    colors.extend(cmap(p_normalized_values))
    cmap = mpl.colors.LinearSegmentedColormap.from_list('my_list', ["#2066a8", "#8ec1da", "#cde1ec"])
    colors.extend(cmap(n_normalized_values))
    

    fig, axes = plt.subplots(1, 2, sharey=True, figsize=figsize)

    sns.barplot(data=stats_bar_group12, x = 'score', y = 'cell_type', ax=axes[0], ci="sd", orient="h", 
                palette=colors, order=order)
    sns.boxplot(x = "value", hue = label_id, y = celltype_id, data = df_long, palette = my_pal, ax=axes[1],
                orient="h", order=order)

    # --- ADDING THE STARS ---
    for i, row in enumerate(stats_bar_group12.itertuples()):
        star_text = row.stars
        if star_text == '': continue # Skip if not significant
        
        x_pos_bar = row.score + (0.1 if row.score > 0 else -0.1)
        axes[0].text(x_pos_bar, i, star_text, 
                     ha='center', va='center', color='black', fontweight='bold')
        max_val = df_long[df_long[celltype_id] == row.cell_type]['value'].max()
        axes[1].text(max_val * 1.05, i, star_text, 
                     ha='left', va='center', color='black', fontweight='bold')
    
    axes[1].set_xlabel("Cell frequency", fontsize = font_size)
    axes[0].set_xlabel("Statistic score", fontsize = font_size)
    axes[0].set_ylabel("", fontsize = font_size)
    plt.xticks(fontsize = font_size)

    axes[0].tick_params(axis="x", labelsize=font_size)
    axes[1].tick_params(axis="x", labelsize=font_size)
    axes[0].tick_params(axis="y", labelsize=font_size)
    axes[1].tick_params(axis="y", labelsize=font_size)
    plt.legend(fontsize = font_size)
    plt.tight_layout()
    if file_path is not None:
        if Path(file_path).is_dir():
            plt.savefig(file_path + "/cell_composition.pdf")
        else:
            print("file path does not exist!")
    plt.show()

def add_trajectory_density(adata = None,
                           sample_id: str = 'Sample',
                           trajectory_id: str = 'Trajectory',
                           adjust: float = 0.3,
                           coord_x: str = 'pxl_col_in_fullres',
                           coord_y: str = 'pxl_row_in_fullres'
                           
                          ):
    names = adata.obs[sample_id].unique()
    all_coords = []
    # define global bandwidth first
    for name in names:
        slide = select_slide(adata, name, sample_id)
        coords = slide.obs[[coord_x, coord_y]].values
        all_coords.append(coords)
    
    all_coords = np.vstack(all_coords).T
    global_kde = gaussian_kde(all_coords)
    
    global_bw = global_kde.factor * adjust
    adata.obs[f"{trajectory_id}_density"] = np.nan
    
    for name in names:
        slide = select_slide(adata, name, sample_id)
        coords = slide.obs[[coord_x, coord_y]].values.T
        scores = slide.obs[trajectory_id].values
        mask = ~np.isnan(scores)
    
        clean_coords = coords[:, mask]
        clean_weights = scores[mask]
    
        kde = gaussian_kde(
            clean_coords,
            weights = clean_weights
        )
        kde.set_bandwidth(bw_method = global_bw)
        density = kde(clean_coords)
        density *= clean_coords.shape[1]

        scaler = MinMaxScaler(feature_range = (0, 1))
        density_scale = scaler.fit_transform(density.reshape(-1, 1))
        
        adata.obs.loc[slide.obs_names[mask], f"{trajectory_id}_density"] = density_scale * 100
    return adata

def get_blueYellow_color(reverse: bool = False):
    archr_blueYellow = [
        # '#f9fa60', '#f0b945', '#2bb4a4', '#0260dd','#4a4194'
        "#F8FA0D",   # Bright Yellow
        "#F6DA23",
        "#F8BA43",
        "#A5BE6A",
        "#2DB7A3",  # Teal transition
        "#1389D2",
        "#0262E0",
        "#343DAE",
        "#352A86"  # Deep Blue
    ]

    if reverse:
        archr_blueYellow = list(reversed(archr_blueYellow))
    
    # Create a Matplotlib colormap
    cmap_blueYellow = LinearSegmentedColormap.from_list("blueYellow", archr_blueYellow)
    cmap_blueYellow.set_bad('#f0f4f9')             
    cmap_blueYellow.set_under('#f0f4f9')
    return cmap_blueYellow

def get_blueRed_color():
    archr_blueRed = [
        "#3361a5", "#2580e1", "#15aafc", "#86cdee", "#c2d4d6", "#ebce86", "#f58d1c", "#cf2424", "#a41f1f"
    ]
    
    # Create a Matplotlib colormap
    cmap_blueRed = LinearSegmentedColormap.from_list("blueRed", archr_blueRed)
    cmap_blueRed.set_bad('#f0f4f9')             
    cmap_blueRed.set_under('#f0f4f9')
    return cmap_blueRed

def compute_TF_activity_kernel_density(adata = None,
                                       TF_name: str = None,
                                       trajectory_cells: bool = False,
                                       TFactivity_name: str = 'TFactivity',
                                       adjust: float = 0.1,
                                       sample_id: str = 'Sample',
                                       trajectory_id: str = 'Trajectory',
                                       coord_x: str = 'pxl_col_in_fullres',
                                       coord_y: str = 'pxl_row_in_fullres'
                                      ):
    
    names = adata.obs[sample_id].unique()
    if trajectory_cells:
        # cell_ids = adata.obs_names[adata.obs[celltype_id].isin(celltypes)]
        cell_ids = adata.obs_names[adata.obs[trajectory_id].notna()]
        adata.obs['TF_name'] = np.nan
        adata.obs.loc[cell_ids, 'TF_name'] = adata.obsm[TFactivity_name].loc[cell_ids, TF_name].values
    else:
        adata.obs['TF_name'] = adata.obsm[TFactivity_name][TF_name].values
    
        
    adata.obs["TF_name_density"] = np.nan
    for name in names:
        slide = select_slide(adata, name, sample_id)
        
        coords = slide.obs[[coord_x, coord_y]].values.T #slide.obsm['spatial'][:, :2].T
        scores = slide.obs['TF_name'].values 
        
        mask = ~np.isnan(scores)
        
        clean_coords = coords[:, mask]
        clean_weights = scores[mask]
            
        pos_weights = np.clip(clean_weights, a_min = 0, a_max = None)
        neg_weights = np.clip(clean_weights, a_min = None, a_max = 0) * -1 
        
        # Calculate densities
        kde_pos = gaussian_kde(clean_coords, weights = pos_weights)
        kde_neg = gaussian_kde(clean_coords, weights = neg_weights)
        
        kde_pos.set_bandwidth(bw_method = adjust)
        kde_neg.set_bandwidth(bw_method = adjust)
        
        # Combine them into a net activity score
        # (High positive = high activation density; Low negative = high repression density)
        net_density = kde_pos(clean_coords) - kde_neg(clean_coords)

        net_density *= clean_coords.shape[1]

        max_abs = np.max(np.abs(net_density))

        if max_abs > 0:
            net_density = (net_density / max_abs) * 2
        adata.obs.loc[slide.obs_names[mask], "TF_name_density"] = net_density

    return adata

def get_trajectory(
    adata,
    trajectory_name="Trajectory",
    activity_key="TFactivity",
    tf_names=None,
    group_every=1,
    log2_norm=False,
    scale_to=None,
    smooth_window=7,
):
    """
    Compute binned TF activity along a trajectory.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing trajectory metadata and TF activity.

    trajectory_name : str
        Column in adata.obs containing trajectory values (0–100).

    activity_key : str
        Key in adata.obsm containing the cell x TF activity matrix.

    tf_names : list, optional
        TF names corresponding to the columns of the activity matrix.
        If None, names are inferred from a DataFrame or generated.

    group_every : float
        Width of trajectory bins, in trajectory units.

    log2_norm : bool
        Whether to apply log2(x + 1) to nonnegative values.

    scale_to : float or None
        Column-wise normalization target. Set None for signed TF scores.

    smooth_window : int or None
        Window size for centered rolling mean. None disables smoothing.

    Returns
    -------
    dict
        {
            "mat": DataFrame (TFs x trajectory bins),
            "smoothMat": DataFrame (TFs x trajectory bins), if requested
        }
    """

    # 1. Validate trajectory metadata
    if trajectory_name not in adata.obs.columns:
        raise ValueError(
            f"Cannot find trajectory '{trajectory_name}' in adata.obs"
        )

    trajectory = pd.to_numeric(
        adata.obs[trajectory_name], errors="raise"
    )

    if ((trajectory.dropna() < 0) |
        (trajectory.dropna() > 100)).any():
        raise ValueError("Trajectory values must be between 0 and 100.")

    if group_every <= 0 or not np.isfinite(group_every):
        raise ValueError("group_every must be a positive finite number.")

    if not np.isclose(100 / group_every, round(100 / group_every)):
        raise ValueError("group_every must evenly divide 100.")

    # 2. Extract TF activity (cells x TFs)
    if activity_key not in adata.obsm:
        raise ValueError(f"Cannot find adata.obsm['{activity_key}']")

    activity = adata.obsm[activity_key]

    if isinstance(activity, pd.DataFrame):
        if tf_names is None:
            tf_names = list(activity.columns)
        activity = activity.to_numpy()

    if tf_names is None:
        tf_names = [f"TF_{i}" for i in range(activity.shape[1])]

    if len(tf_names) != activity.shape[1]:
        raise ValueError("Length of tf_names does not match TF columns.")

    if sparse.issparse(activity):
        activity = activity.toarray()
    else:
        activity = np.asarray(activity)

    # 3. Create trajectory bins
    breaks = np.linspace(
        0, 100, int(round(100 / group_every)) + 1
    )

    valid = trajectory.notna().to_numpy()
    traj_values = trajectory.to_numpy()[valid]
    activity = activity[valid, :]

    # (lower, upper] matching R; excludes trajectory == 0
    bin_idx = np.searchsorted(
        breaks, traj_values, side="left"
    ) - 1

    n_bins = len(breaks) - 1

    counts = np.bincount(
        bin_idx[bin_idx >= 0], minlength=n_bins
    )

    bin_names = [
        f"T.{lo:g}_{hi:g}"
        for lo, hi in zip(breaks[:-1], breaks[1:])
    ]

    if np.any(counts == 0):
        empty = [
            bin_names[i] for i in np.flatnonzero(counts == 0)
        ]
        raise ValueError(
            "Cannot compute means for empty trajectory bins: "
            + ", ".join(empty)
            + ". Try increasing group_every."
        )

    # 4. Mean TF activity for each trajectory bin
    print("Creating Trajectory Group Matrix...")

    group_mat = np.column_stack([
        activity[bin_idx == i].mean(axis=0)
        for i in range(n_bins)
    ])

    # 5. Optional depth normalization
    if scale_to is not None:
        if np.any(group_mat < 0):
            print(
                "Negative TF activity detected. "
                "Skipping depth normalization."
            )
        else:
            sums = group_mat.sum(axis=0)
            if np.any(sums == 0):
                raise ValueError(
                    "Cannot normalize bins with zero total activity."
                )
            group_mat = group_mat / sums * scale_to

    # 6. Optional log2 normalization
    if log2_norm:
        if np.any(group_mat < 0):
            print(
                "Negative TF activity detected. "
                "Skipping log2 normalization."
            )
        else:
            group_mat = np.log2(group_mat + 1)

    # 7. Store unsmoothed matrix (TFs x bins)
    mat = pd.DataFrame(
        group_mat,
        index=tf_names,
        columns=bin_names,
    )

    result = {"mat": mat}

    # 8. Centered rolling mean along trajectory
    if smooth_window is not None:
        if (
            not isinstance(smooth_window, (int, np.integer))
            or smooth_window < 1
            or smooth_window % 2 == 0
        ):
            raise ValueError(
                "smooth_window must be a positive odd integer."
            )

        print("Smoothing...")

        # Edge windows use the available neighboring bins
        smooth_mat = mat.T.rolling(
            window=smooth_window,
            center=True,
            min_periods=1,
        ).mean().T

        result["smoothMat"] = smooth_mat

    return result


def trajectory_heatmap(
    trajectory,
    assay = "smoothMat",
    var_cutoff = 0,
    max_features = None,
    limits = (-2, 2),
    scale_rows = True,
    row_order = None,
    label_rows = True,
    label_top = 50,
    cmap = "RdYlBu_r",
    name = "TF activity",
    figsize = (12, 10),
    return_mat = False,
    save_path: str = None
):
    """
    Python equivalent of ArchR TrajectoryHeatmap.

    Parameters
    ----------
    trajectory : dict or pandas.DataFrame
        Output of get_trajectory(), or a TF x trajectory-bin DataFrame.

    assay : str
        "smoothMat" or "mat" from get_trajectory() output.

    var_cutoff : float or None
        Fraction of lowest-variance TFs to exclude.
        0 keeps all TFs; 0.2 keeps the top 80%.

    max_features : int or None
        Maximum number of TFs retained.

    limits : tuple
        Lower and upper clipping limits after Z-score normalization.

    scale_rows : bool
        Whether to Z-score each TF across trajectory bins.

    row_order : list, optional
        Explicit ordered list of TF names.

    label_rows : bool
        Show TF labels.

    label_top : int or None
        Number of highest-variance TFs to label.

    cmap : str
        Matplotlib colormap.

    name : str
        Heatmap colorbar label.

    figsize : tuple
        Figure size.

    return_mat : bool
        If True, return the processed and ordered matrix.

    Returns
    -------
    pandas.DataFrame, optional
        Processed TF x trajectory-bin matrix when return_mat=True.
    """

    # 1. Extract activity matrix
    if isinstance(trajectory, dict):
        if assay not in trajectory:
            raise ValueError(
                f"'{assay}' is not available. "
                f"Available: {list(trajectory.keys())}"
            )
        mat = trajectory[assay].copy()

    elif isinstance(trajectory, pd.DataFrame):
        mat = trajectory.copy()

    else:
        raise TypeError(
            "trajectory must be a dict or pandas DataFrame."
        )

    # 2. Remove rows with NA
    na_rows = mat.isna().any(axis=1)

    if na_rows.any():
        print(f"Removing {na_rows.sum()} TFs containing NA values...")
        mat = mat.loc[~na_rows]

    # 3. Remove zero-variance rows
    row_sd = mat.std(axis=1, ddof=1)
    zero_var = row_sd.isna() | (row_sd == 0)

    if zero_var.any():
        print(f"Removing {zero_var.sum()} TFs without variation...")
        mat = mat.loc[~zero_var]

    if mat.empty:
        raise ValueError("No Features Remaining!")

    # 4. Calculate row variance
    row_var = mat.var(axis=1, ddof=1)

    # Percentile rank of variance (ArchR-style varQ)
    var_q = row_var.rank(method="average", pct=True)

    ordered_var = False

    # 5. Select top variable TFs
    if row_order is None:
        mat = mat.loc[
            var_q.sort_values(ascending=False).index
        ]

        ordered_var = True
        n = len(mat)
        if var_cutoff is not None:
            if not 0 <= var_cutoff <= 1:
                raise ValueError(
                    "var_cutoff must be between 0 and 1."
                )

            n = min(
                n,
                int(np.floor((1 - var_cutoff) * len(mat)))
            )

        if max_features is not None:
            n = min(n, int(max_features))

        mat = mat.iloc[:n]

    if mat.empty:
        raise ValueError("No Features Remaining!")

    # 6. Identify TFs to label (highest variance)
    if label_top is not None and label_top > 0:
        if ordered_var:
            idx_label = mat.index[:label_top].tolist()
        else:
            idx_label = (
                var_q.reindex(mat.index)
                .sort_values(ascending=False)
                .index[:label_top]
                .tolist()
            )
    else:
        idx_label = []

    # 7. Row-wise Z-score normalization
    if scale_rows:
        row_mean = mat.mean(axis=1)
        row_sd = mat.std(axis=1, ddof=1)

        mat = mat.sub(row_mean, axis=0).div(row_sd, axis=0)

        # Clip extreme scores
        mat = mat.clip(
            lower=min(limits),
            upper=max(limits)
        )

    # 8. Order TFs by position of maximum activity
    if row_order is not None:
        if all(isinstance(x, str) for x in row_order):
            missing = set(row_order) - set(mat.index)

            if missing:
                raise ValueError(
                    f"TFs not found in matrix: {missing}"
                )

            mat = mat.loc[row_order]

        else:
            # Integer positional ordering (Python 0-based)
            mat = mat.iloc[row_order]

    else:
        peak_position = np.argmax(mat.to_numpy(), axis=1)

        order_idx = np.argsort(
            peak_position, kind="stable"
        )

        mat = mat.iloc[order_idx]

    # 9. Plot heatmap
    fig, ax = plt.subplots(figsize=figsize)
    fig.subplots_adjust(right=0.80)

    sns.heatmap(
        mat,
        cmap=cmap,
        vmin=min(limits) if scale_rows else None,
        vmax=max(limits) if scale_rows else None,
        center=0 if scale_rows else None,
        xticklabels=False,
        yticklabels=False,
        cbar_kws={"location": "top", "orientation": "horizontal", "pad": 0.08, "label": name},
        ax=ax,
    )

    ax.set_xlabel("Pseudotime")
    ax.set_ylabel("")
    ax.set_title(name)


    # 10. Label top variable TFs without overlap
    # 10. Evenly spaced TF labels with connector lines
    if label_rows and idx_label:
    
        # Get original row positions
        label_info = [
            (tf, mat.index.get_loc(tf) + 0.5)
            for tf in idx_label
            if tf in mat.index
        ]
    
        # Sort labels by original heatmap position
        label_info.sort(key=lambda x: x[1])
    
        n_labels = len(label_info)
        n_rows = mat.shape[0]
        n_cols = mat.shape[1]
    
        if n_labels > 0:
    
            # Evenly distribute labels vertically
            label_positions = np.linspace(
                0.5,
                n_rows - 0.5,
                n_labels
            )
    
            # Keep heatmap width unchanged; put labels outside axes
            x_start = n_cols
            x_end = n_cols + 3
            x_text = n_cols + 4
    
            for (tf, original_y), label_y in zip(
                label_info, label_positions
            ):
    
                # Connector from true TF row to label position
                ax.plot(
                    [x_start, x_end],
                    [original_y, label_y],
                    color="gray",
                    linewidth=0.6,
                    alpha=0.7,
                    clip_on=False,
                )
    
                # TF label
                ax.text(
                    x_text,
                    label_y,
                    tf,
                    fontsize=8,
                    ha="left",
                    va="center",
                    clip_on=False,
                )

    # Reduce crowded x-axis labels
    n_cols = mat.shape[1]
    step = max(1, n_cols // 10)

    plt.tight_layout()
    if save_path is not None:
        if Path(save_path).is_dir():
            plt.savefig(save_path + "/trajectpry_heatmap.pdf")
        else:
            print("save path does not exist!")
    plt.show()

    if return_mat:
        return mat

def plot_trajectpry_activity_heatmap(adata = None,
                            assay: str ="smoothMat",
                            trajectory_name: str ="Trajectory",
                            activity_key: str ="TFactivity",
                            tf_names = None,
                            group_every = 1,
                            log2_norm = False,
                            scale_to = None,
                            smooth_window = 7,
                            var_cutoff = 0,
                            max_features = None,
                            limits = (-2, 2),
                            scale_rows = True,
                            row_order = None,
                            label_rows = True,
                            label_top = 50,
                            cmap = 'RdYlBu_r',
                            name = 'TF activity',
                            figsize = (4, 10),
                            return_mat = False,
                            save_path = None
                           ):
    result = get_trajectory(adata,
                            trajectory_name,
                            activity_key,
                            tf_names,
                            group_every,
                            log2_norm,
                            scale_to,
                            smooth_window
                           )
    heatmap_mat = trajectory_heatmap(result,
                                     assay,
                                     var_cutoff,
                                     max_features,
                                     limits,
                                     scale_rows,
                                     row_order,
                                     label_rows,
                                     label_top,
                                     cmap,
                                     name,
                                     figsize,
                                     return_mat,
                                     save_path
                                    )

def adjust_p_values(p_values):
    p = np.asarray(p_values)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]

def plot_top_trajectory_activity_terminal_fold_change(adata = None,
                                                      trajectory_name: str ="Trajectory",
                                                      activity_key: str ="TFactivity",
                                                      tf_names = None,
                                                      group_every = 1,
                                                      log2_norm = False,
                                                      scale_to = None,
                                                      smooth_window = 7,
                                                      num_top = 10,
                                                      adj_pval_thr = 0.01,
                                                      return_mat = False,
                                                      save_path = None,
                                                      figsize = (3,6),
                                                      cmap_color = 'viridis_r',
                                                      fontsize = 14,
                                                      point_size: int = 150
                                                     ):
    result = get_trajectory(adata,
                            trajectory_name,
                            activity_key,
                            tf_names,
                            group_every,
                            log2_norm,
                            scale_to,
                            smooth_window
                           )
    TFs_activities_heatmap = result['mat']
    
    zscoreTFs_activities_heatmap = stats.zscore(TFs_activities_heatmap.transpose().copy()).transpose()

    zscoreTFs_activities_heatmap = pd.DataFrame(zscoreTFs_activities_heatmap)
    zscoreTFs_activities_heatmap.index = TFs_activities_heatmap.index
    zscoreTFs_activities_heatmap.columns = TFs_activities_heatmap.columns

    tf_activity_score = np.subtract(zscoreTFs_activities_heatmap[zscoreTFs_activities_heatmap.columns[-1]],
                                zscoreTFs_activities_heatmap[zscoreTFs_activities_heatmap.columns[0]]) 
    z_score = tf_activity_score
    p_values = norm.sf(abs(z_score)) * 2

    detf = pd.DataFrame()
    detf['TF'] = tf_activity_score.index
    detf['values'] = z_score.values
    detf['p_values'] = p_values
    detf['adj_p_values'] = adjust_p_values(p_values)

    detf_sorted = detf[detf['adj_p_values'] < adj_pval_thr].sort_values('values', ascending = False)
    detf = pd.concat([detf_sorted.head(num_top), detf_sorted.tail(num_top)])
    detf = detf.sort_values('values')

    plt.style.use("default")

    fig, ax = plt.subplots(figsize = figsize, dpi = 100)
    
    bar_df = detf.sort_values("values")  # optional but looks nicer
    y_positions = np.arange(len(bar_df))
    
    ax.hlines(
        y = y_positions,
        xmin = 0,
        xmax = detf["values"],
        color = "black",
        linewidth = 1,
        zorder = 1
    )
    
    norm_colors = mcolors.Normalize(
        vmin = 0,
        vmax = -np.log10(detf['adj_p_values']).max() + 5
    )
    
    cmap = cm.get_cmap(cmap_color)
    
    sca = ax.scatter(
        detf["values"],
        y_positions,
        c = -np.log10(detf['adj_p_values']),
        cmap = cmap,
        norm = norm_colors,
        s = point_size
    )
    
    cbar = fig.colorbar(sca,
                        ax = ax,
                        fraction = 0.04,   # thickness (default ≈0.15)
                        pad = 0.02,        # space between plot and colorbar
                        aspect = 30,       # length/width ratio
                        shrink = 0.4
                       )
    cbar.set_label("-log adj p.value")
    
    ax.set_yticks(y_positions)
    ax.set_yticklabels(detf["TF"])
    
    ax.set_xlabel("Terminal Fold change", fontsize = fontsize)
    
    sns.despine(bottom=False, left=True)
    plt.yticks(fontsize = fontsize - 2)
    if save_path is not None:
        if Path(save_path).is_dir():
            plt.savefig(save_path + "/trajectpry_activity_terminal_fold_change.pdf")
        else:
            print("save path does not exist!")
    plt.show()

    if return_mat:
        return mat

def compute_gene_activity_kernel_density(adata = None,
                                         gene_name: str = None,
                                         sample_id: str = 'Sample',
                                         global_bw: float = 0.01,
                                         gene_activity_name: str = 'gene_activity',
                                         coord_x: str = 'pxl_col_in_fullres',
                                         coord_y: str = 'pxl_row_in_fullres'
                                        ):

    names = adata.obs[sample_id].unique()
    if gene_name in adata.obsm[gene_activity_name].columns:
        for name in names:
            slide = select_slide(adata, name, sample_id)

            slide.obs['gene_name'] = adata.obsm[gene_activity_name][gene_name]
            slide.obs["gene_name_density"] = np.nan
                
            coords = slide.obs[[coord_x, coord_y]].values.T
            scores = slide.obs['gene_name'].values 
            
            mask = ~np.isnan(scores)
            
            clean_coords = coords[:, mask]
            clean_weights = scores[mask]
                
            kde = gaussian_kde(
                clean_coords,
                weights=clean_weights
            )

            kde.set_bandwidth(bw_method=global_bw)
            density = kde(clean_coords)
            density *= clean_coords.shape[1]
            
            slide.obs.loc[slide.obs_names[mask], "gene_name_density"] = density
            
            slide.obs["gene_name_density"] = (
                        slide.obs
                        .groupby("Sample")["gene_name_density"]
                        .transform(lambda x: -2 + 4 * (x - x.min()) / (x.max() - x.min()))
                    )
            adata.obs.loc[slide.obs_names, "gene_name"] = slide.obs['gene_name']
            adata.obs.loc[slide.obs_names, "gene_name_density"] = slide.obs['gene_name_density']
    else:
        print("The gene does not exist!")
        return None
    return adata

def plot_gene_activity_marker_heatmap(adata = None,
                                      gene_set: list = None,
                                      celltypes: list = None,
                                      celltype_name: str = 'Clusters2',
                                      gene_activity_name: str = 'gene_activity',
                                      save_path: str = None,
                                      figsize = (4,8),
                                      return_mat = False
                                     ):
    
    gene_matrix = adata.obsm[gene_activity_name]
    sel_gene_matrix = gene_matrix[gene_set]
    sel_gene_matrix['celltype'] = adata.obs[celltype_name]
    data = sel_gene_matrix.groupby('celltype').mean().T
    data = data[celltypes]
        
    row_order = np.argsort(np.argmax(data.values, axis=1))
    
    data_sorted = data.iloc[row_order]
    g = sns.clustermap(data_sorted,
                       row_cluster = False,
                       col_cluster = False,
                       z_score = 0,  
                       figsize = figsize,
                       yticklabels = True,
                       vmin = -3,
                       vmax = 3,
                       cmap = get_blueYellow_color(reverse = True), 
                       center = 0,
                       cbar_pos = (0.05, .35, .02, .2),
                       cbar_kws = {"label": 'Z-score Gene activity'}
                  )
    g.cax.yaxis.set_label_coords(-2.2,0.5)
    g.ax_heatmap.set_xlabel('')
    if save_path is not None:
        if Path(save_path).is_dir():
            plt.savefig(save_path + "/gene_activity_marker_heatmap.pdf")
        else:
            print("save path does not exist!")
    plt.show()
    if return_mat:
        return data_sorted

def count_samples_clone(adata, 
                        sample_id: str = 'Sample',
                        celltype_id: str = 'Clusters2',
                        clone_id: str = 'clone'
                       ):
    count_clones = adata.obs.groupby([sample_id, celltype_id, clone_id]).count().iloc[:, 0]
    count_clones = count_clones.reset_index()

    count_clones = count_clones[
        ~count_clones.iloc[:, -1].eq(0).groupby([count_clones[sample_id], count_clones[celltype_id]]).transform('any')
    ]
    count_clones = count_clones.reset_index(drop = True)
    return count_clones

def plot_quality_control(adata = None,
                         features: list = ['TSSEnrichment', 'nFrags', 'FRIP'],
                         sample_id: str = 'Sample',
                         sort_feature: str = 'FRIP',
                         label_id: str = 'label',
                         save_path: str = None
                        ):

    colors = px.colors.qualitative.Plotly

    df = adata.obs[[sample_id, sort_feature, label_id]]
    samples_order = df.groupby(sample_id).agg(
        Feature = (sort_feature, 'mean'),
        sort_label = (label_id, 'first')  
    ).sort_values(['sort_label', 'Feature']).index
    
    sample_colors = {
        sample: colors[i % len(colors)]
        for i, sample in enumerate(samples_order)
    }
    
    fig = make_subplots(
        rows = 1,
        cols = len(features),
        subplot_titles = features
    )
    
    for col, feature in enumerate(features, start = 1):
        for sample in samples_order:
            mask = adata.obs[sample_id] == sample
            color = sample_colors[sample]
    
            fig.add_trace(
                go.Violin(
                    x = adata.obs.loc[mask, sample_id],
                    y = adata.obs.loc[mask, feature],
                    name = sample,
                    legendgroup = sample,
                    line_color = color,
                    fillcolor = f"rgba({int(color[1:3], 16)}, {int(color[3:5], 16)}, {int(color[5:7], 16)}, 0.2)",
                    opacity = 0.7,
                    width = 0.5,
                    showlegend = (col == 1),
                    box_visible = True,
                    meanline_visible = True,
                    points = False
                ),
                row = 1,
                col = col
            )
            fig.update_xaxes(tickangle = 45)
    
    fig.update_layout(
        width = 300 * len(features),
        height = 500,
        showlegend = True,
        violinmode = 'group',
        autosize = False,
        template = 'plotly_white',
        margin = dict(l=10, r=10, t=100, b=50)
    )
    if save_path is not None:
        if Path(save_path).is_dir():
            fig.write_image(save_path + "/quality_control_plot.pdf")
        else:
            print("save path does not exist!")
    fig.show()

def custom_deg_df(adata, group):
    results = adata.uns['rank_genes_groups']
    
    df = pd.DataFrame({
        'gene': results['names'][group],
        'logfoldchange': results['logfoldchanges'][group],
        'pvals_adj': results['pvals_adj'][group],
        'score': results['scores'][group]
    })
    
    return df

def Differential_TF_activity_clone_vs_nonClone(adata = None,
                                               sample_name: str = None,
                                               celltypes: list = None,
                                               clone_id: str = 'clone',
                                               sample_id: str = 'Sample',
                                               celltype_id: str = 'Clusters2',
                                               TFactivity_id: str = 'TFactivity',
                                               save_path: str = None,
                                               font_size: int = 12,
                                               number_top: int = 10
                                              ):
    
    cell_ids = adata.obs_names[ (adata.obs[sample_id] == sample_name) & (adata.obs[celltype_id].isin(celltypes))]
    TFadata = ad.AnnData(adata.obsm[TFactivity_id].loc[cell_ids])
    TFadata.obs = adata.obs.loc[cell_ids]
    TFadata.obs[clone_id] = TFadata.obs[clone_id].astype(str)
    print(TFadata.obs[clone_id].value_counts())
    if TFadata.obs[clone_id].value_counts()[1] == 1:
        print("Not enough cells!")
    else:
        sc.tl.rank_genes_groups(
            TFadata, 
            groupby = clone_id, 
            groups = ['1'], 
            reference = '0', 
            method = 'wilcoxon'  
        )
        
        rec_array = TFadata.uns['rank_genes_groups']['names']
        group_genes = rec_array['1'] 
        result = custom_deg_df(TFadata, '1')
        
        detf = pd.concat([result.sort_values('score', ascending = False).head(number_top),
                          result.sort_values('score', ascending = False).tail(number_top)])
        print(list(detf['gene'].values))
        
        plt.style.use("default")
    
        fig, ax = plt.subplots(figsize = (3,5), dpi = 100)
        
        bar_df = detf.sort_values("score")  
        y_positions = np.arange(len(bar_df))
        
        ax.hlines(
            y = y_positions,
            xmin = 0,
            xmax = bar_df["score"],
            color = "black",
            linewidth = 1,
            zorder = 1
        )
        
        norm = mcolors.Normalize(
            vmin = 0,
            vmax = 10
        )
        
        cmap = cm.get_cmap("viridis_r")
        
        sca = ax.scatter(
            bar_df["score"],
            y_positions,
            c = -np.log10(bar_df['pvals_adj']),
            cmap = cmap,
            norm = norm,
            s = 150
        )
        
        cbar = fig.colorbar(sca,
                            ax = ax,
                            fraction = 0.04,   
                            pad = 0.02,        
                            aspect = 30,       
                            shrink = 0.4
                           )
        cbar.set_label("-log adj p.value")

        ax.set_xlabel('Fold change\nnon-clone vs clone')
        ax.set_yticks(y_positions)
        ax.set_yticklabels(bar_df["gene"])
    
        sns.despine(bottom = False, left = True)
        plt.yticks(fontsize = font_size)
        if save_path is not None:
            if Path(save_path).is_dir():
                if len(celltypes) == 1:
                    plt.savefig(save_path + f"/Differential_TF_activity_clone_vs_nonClone_{sample_name}_{celltypes[0]}.pdf")
                else:
                    plt.savefig(save_path + f"/Differential_TF_activity_clone_vs_nonClone_{sample_name}_PTs.pdf")
            else:
                print("save path does not exist!")
        plt.show()

def plot_clone_spatial_image(adata = None,
                             sample_name: str = None,
                             celltypes: list = None,
                             clone_id: str = 'clone',
                             sample_id: str = 'Sample',
                             celltype_id: str = 'Clusters2',
                             spatial_id: str = 'spatial_rotated',
                             save_path: str = None,
                             spot_size: int = 50,
                             fig_size = [6, 6],
                             pal = None
                            ):

    slide = select_slide(adata, sample_name, sample_id)
    slide.obs['plot_PT_clone'] = 'Other'
    PT_cells = slide.obs_names[slide.obs[celltype_id].isin(celltypes)]
    PT_clone_cells = slide.obs_names[ (slide.obs[celltype_id].isin(celltypes)) & (slide.obs[clone_id] == 1)]
    if len(celltypes) == 1:
        slide.obs.loc[PT_cells, 'plot_PT_clone'] = f'Non-clonal {celltypes[0]}'
        slide.obs.loc[PT_clone_cells, 'plot_PT_clone'] = f'Clonal {celltypes[0]}'
    else:
        slide.obs.loc[PT_cells, 'plot_PT_clone'] = 'Non-clonal PTs'
        slide.obs.loc[PT_clone_cells, 'plot_PT_clone'] = 'Clonal PTs'

    if pal is None:
        if len(celltypes) == 1:
            pal = {f'Clonal {celltypes[0]}': 'tab:red',
                   f'Non-clonal {celltypes[0]}': 'tab:green',
                   'Other': '#f0f4f9'}
        else:
            pal = {'Clonal PTs': 'tab:red',
                   'Non-clonal PTs': 'tab:green',
                   'Other': '#f0f4f9'}

    with mpl.rc_context({'figure.figsize': fig_size,
                         'axes.facecolor': 'white'}):
        sc.pl.embedding(adata = slide, basis = spatial_id,
                        size = spot_size, color = ['plot_PT_clone'],
                        na_color = '#f0f4f9', palette = pal)
    if save_path is not None:
        if Path(save_path).is_dir():
            if len(celltypes) == 1:
                with mpl.rc_context({'figure.figsize': fig_size,
                                     'axes.facecolor': 'white'}):
                    sc.pl.embedding(adata = slide, basis = spatial_id,
                                    size = spot_size, color = ['plot_PT_clone'],
                                    na_color = '#f0f4f9', palette = pal,
                                    show = False).figure.savefig(save_path + f"clone_spatial_{sample_name}_{celltypes[0]}.pdf",
                                                               bbox_inches = 'tight', dpi = 300)
            else:
                with mpl.rc_context({'figure.figsize': fig_size,
                                     'axes.facecolor': 'white'}):
                    sc.pl.embedding(adata = slide, basis = spatial_id,
                                    size = spot_size, color = ['plot_PT_clone'],
                                    na_color = '#f0f4f9', palette = pal,
                                    show = False).figure.savefig(save_path + f"clone_spatial_{sample_name}_PTs.pdf",
                                                               bbox_inches = 'tight', dpi = 300)
        else:
            print("save path does not exist!")
