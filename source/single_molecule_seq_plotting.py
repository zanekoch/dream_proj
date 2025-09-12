import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm
import warnings
import colorcet as cc
from scipy.stats import ttest_ind
import re
DIPLOID_MOUSE_GENOME = 5.4e9

def plot_signature_activities(activity_df, signature_type='signatures', colors=cc.glasbey_cool, mut_cols=None, figsize=(6, 4), fig=None, axes=None, normalize=False):
    """
    create a stacked barplot of signature activities with subplots for each group

    Parameters:
    -----------
    activity_df : pd.DataFrame
        DataFrame with signature activities containing 'sample', 'group', and signature columns
    signature_type : str
        Type of signatures for labeling (e.g., 'SBS signatures', 'ID signatures')
    colors : list or colormap
        Colors for the signatures
    mut_cols : list, optional
        List of mutation columns to plot
    figsize : tuple
        Figure size for the plot
    fig : matplotlib figure object (optional)
        Existing figure to use
    axes : tuple of matplotlib axes objects (optional)
        Tuple of (ax1, ax2) existing axes to use
    normalize : bool, optional
        If True, show within-sample proportions instead of absolute counts (default False)

    Returns:
    --------
    fig, (ax1, ax2) : matplotlib figure and axes
    """
    # get signature columns (exclude sample, group, and age_days)
    if mut_cols is None:
        cols = [col for col in activity_df.columns if col not in ['sample', 'group', 'age_days', 'effective_coverage', 'number_of_effective_read_families']]
    else:
        cols = mut_cols

    # create the stacked barplot with subplots for each group
    if axes is None:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize, sharey=False, gridspec_kw={'width_ratios': [5, 4]})
    else:
        ax1, ax2 = axes
        fig = ax1.figure
        
    pivot_df = activity_df.set_index('sample')[cols]
    
    # normalize to proportions if requested
    if normalize:
        pivot_df = pivot_df.div(pivot_df.sum(axis=1), axis=0)
    
    # sort each group separately by their summed values
    wt_samples = activity_df.query("group == 'control'")['sample']
    ko_samples = activity_df.query("group == 'dream_ko'")['sample']
    wt_pivot = pivot_df.loc[wt_samples]
    ko_pivot = pivot_df.loc[ko_samples]

    wt_order = wt_pivot.sum(axis=1).sort_values(ascending=False).index
    ko_order = ko_pivot.sum(axis=1).sort_values(ascending=False).index

    # plot control group
    wt_pivot_ordered = wt_pivot.loc[wt_order]
    wt_pivot_ordered.plot(kind='bar', stacked=True, ax=ax1, color=colors, edgecolor='black', linewidth=1)
    ax1.set_xlabel('Control')
    # update y-axis label based on normalization
    ylabel = f'{signature_type} (proportion)' if normalize else f'{signature_type}'
    ax1.set_ylabel(ylabel)
    # move the legend above the plot with two columns
    ax1.legend(title='Signature', bbox_to_anchor=(0.5, 1.4), loc='upper center', ncol=2)
    # plot dream_ko group
    ko_pivot_ordered = ko_pivot.loc[ko_order]
    ko_pivot_ordered.plot(kind='bar', stacked=True, ax=ax2, color=colors, edgecolor='black', linewidth=1)
    ax2.set_xlabel('DREAM KO')
    ax2.set_ylabel(ylabel)
    # remove the legend from the second plot
    ax2.legend_.remove()
    
    

    # manually set y-axis range to be the same for both plots
    max_y = max(wt_pivot_ordered.sum(axis=1).max(), ko_pivot_ordered.sum(axis=1).max())
    # for proportions, max should be 1.0, but we add a small buffer
    if normalize:
        ax1.set_ylim(0, 1.05)
        ax2.set_ylim(0, 1.05)
    else:
        ax1.set_ylim(0, max_y*1.1)
        ax2.set_ylim(0, max_y*1.1)

    sns.despine()

    return fig, (ax1, ax2)

def fit_mutation_burden_model(
    df,
    mut_col,
    group_col='group',
    sample_col='sample',
    cov_col='effective_coverage',
    fam_col='number_of_effective_read_families',
    formula=None,
    offset_col=None,
    exclude_samples=None
    ):
    """
    fit negative binomial regression model for mutation burden

    Parameters:
    - df: pandas DataFrame with mutation burden and covariates
    - mut_col: str, column name for mutation burden to model
    - group_col: str, column name for group (default 'group')
    - sample_col: str, column name for sample (default 'sample')
    - cov_col: str, column name for effective coverage (default 'effective_coverage')
    - fam_col: str, column name for number of effective read families (default 'number_of_effective_read_families')
    - formula: str, formula for statsmodels GLM (if None, uses default)
    - offset_col: str or None, column name for offset variable
    - exclude_samples: list or None, samples to exclude from analysis

    Returns:
    - model: fitted statsmodels GLM model
    - plot_df: processed DataFrame used for modeling
    """
    import numpy as np
    import pandas as pd
    import statsmodels.formula.api as smf
    import statsmodels.api as sm
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        plot_df = df.copy()
        if exclude_samples is not None:
            plot_df = plot_df[~plot_df[sample_col].isin(exclude_samples)]
        else:
            plot_df = plot_df.copy()
        
        # negative binomial regression
        if formula is None:
            formula = f"{mut_col} ~ age_days*{group_col} + {cov_col} + {fam_col}"
            
        # scan across alpha values to find the one that minimizes deviance/df
        results = []
        for alpha in np.logspace(-4, 1, 1000):
            model = smf.glm(
                formula,
                data=plot_df, family=sm.families.NegativeBinomial(alpha=alpha),
                offset=np.log(plot_df[offset_col]) if offset_col is not None else None
            ).fit(cov_type='HC0')
            results.append({'alpha': alpha, 'deviance/df': model.deviance / model.df_resid})

        results_df = pd.DataFrame(results)
        chosen_alpha = results_df.iloc[np.argmin(np.abs(results_df['deviance/df'] - 1))]['alpha']
        model = smf.glm(
            formula,
            data=plot_df, family=sm.families.NegativeBinomial(alpha=chosen_alpha),
            offset=np.log(plot_df[offset_col]) if offset_col is not None else None
        ).fit(cov_type='HC0')
    return model, plot_df

def plot_mutation_burden(
    model,
    plot_df,
    mut_col,
    group_col='group',
    ylabel='SBS mutations per cell',
    group_order=['control', 'dream_ko'],
    plot_predicted=False,
    palette=None,
    ax=None
    ):
    """
    plot mutation burden per cell by group

    Parameters:
    - model: fitted statsmodels GLM model
    - plot_df: pandas DataFrame with mutation burden and covariates (with predicted values)
    - mut_col: str, column name for mutation burden to plot
    - group_col: str, column name for group (default 'group')
    - ylabel: str, label for y-axis
    - group_order: list, order of groups for plotting
    - plot_predicted: bool, whether to plot predicted values instead of actual
    - palette: list or None, colors for stripplot
    - ax: matplotlib axis or None

    Returns:
    - fig, ax: matplotlib figure and axis
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    from statannotations.Annotator import Annotator

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    else:
        fig = ax.figure
    
    # determine what to plot
    if plot_predicted:
        y_col = 'predicted_muts'
        plot_df['eff_coverage_tens_of_mbs'] = DIPLOID_MOUSE_GENOME / 1e7
        plot_df['n_effective_rf_millions'] = plot_df['n_effective_rf_millions'].mean()
        plot_df['predicted_muts'] = model.predict(plot_df)# + np.exp(model.params['Intercept'])
    else:
        y_col = mut_col
    
    # pointplot
    sns.pointplot(
        data=plot_df, x=group_col, y=y_col,
        order=group_order,
        linewidth=1, ax=ax,  err_kws={'color': 'black', 'linewidth': 1.3},
        capsize=0.1, errorbar='ci', dodge=True, color = 'black',  marker="D", markersize=5
    )
    
    # stripplot
    if palette is None:
        magma = plt.colormaps['magma']
        palette = [magma(3.33/10), magma(6.66/10)]
    # fix: assign hue to x variable and set legend=False to avoid deprecation warning
    sns.stripplot(
        data=plot_df, x=group_col, y=y_col,
        hue=group_col, palette=palette, s=5, order=group_order, ax=ax, legend=False,
        edgecolor='black', linewidth=1, alpha=0.75, jitter=0.1
    )
    
    sns.despine()
    ax.set_ylabel(ylabel)
    # set x labels
    # fix: set ticks before set_ticklabels to avoid warning
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Control\n(p107D/D;\np130fl/fl)', 'DREAM K.O.\n(p107D/D;\np130–/–)'])
    ax.set_xlabel('')
    
    # add statistical annotation from model
    # extract p-value for the group coefficient
    group_coef_names = ['group[T.dream_ko]', "C(group, Treatment(reference='control'))[T.dream_ko]"]
    group_pvalue = None
    
    for coef_name in group_coef_names:
        if coef_name in model.pvalues.index:
            group_pvalue = model.pvalues[coef_name]
            break
    
    if group_pvalue is not None:
        # set up annotator for statistical annotation
        pairs = [(group_order[0], group_order[1])]  # ('control', 'dream_ko')
        annotator = Annotator(ax, pairs, data=plot_df, x=group_col, y=y_col, order=group_order)
        
        # configure to use custom p-values (no statistical test)
        annotator.configure(text_format='simple', loc='outside', fontsize=6)
        
        # set the custom p-value and apply annotation without performing statistical test
        annotator.set_pvalues([group_pvalue])
        annotator.annotate()

    return fig, ax, plot_df


def plot_group_coefficients(activity_df, mut_cols, formula_template, color_map=None, figsize=(2,2), fig=None, ax=None):
    """
    Plot group coefficients from GLM models for multiple mutation types.
    
    Parameters:
    - activity_df: DataFrame with mutation data
    - mut_cols: list of mutation column names
    - formula_template: string template for GLM formula (should contain {prop_col} placeholder)
    - color_map: dict mapping mutation types to colors
    - figsize: tuple for figure size
    - fig: matplotlib figure object (optional)
    - ax: matplotlib axes object (optional)
    
    Returns:
    - fig, ax: matplotlib figure and axes objects
    """
    if color_map is None:
        color_map = {}
    activity_df = activity_df.copy(deep=True)
    if not 'total_mutations_incl_ffpe' in activity_df.columns:
        activity_df['total_mutations_incl_ffpe'] = activity_df[mut_cols].sum(axis=1)
    
    models = {}
    for mut_col in mut_cols:
        prop_col = f'{mut_col}_prop'
        activity_df[prop_col] = activity_df[mut_col] / activity_df['total_mutations_incl_ffpe']

        # use var_weights=total_mutations so the GLM knows the number of trials
        formula = formula_template.format(prop_col=prop_col)
        model = smf.glm(
            formula=formula,
            data=activity_df,
            family=sm.families.Binomial(),
            var_weights=activity_df['total_mutations_incl_ffpe']
        ).fit(cov_type='HC0')  # robust SEs to guard against overdispersion
        models[mut_col] = model

    # pick only the main-effect group coefficient (exclude interactions)
    group_term_regex = r"C\(group, Treatment\(reference='control'\)\)\[T\.dream_ko\]$"

    rows = []
    for mut_col, model in models.items():
        terms = [t for t in model.params.index if re.fullmatch(group_term_regex, t)]
        if not terms:
            continue
        t = terms[0]
        ci = model.conf_int().loc[t]
        rows.append({
            'mut_col': mut_col,
            'coef': model.params[t],
            'pvalue': model.pvalues[t],
            'lo': ci[0],
            'hi': ci[1],
        })

    group_df = pd.DataFrame(rows).sort_values('coef').reset_index(drop=True)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    y = np.arange(len(group_df))

    # draw vertical bars from 0 to coefficient value
    for i, row in group_df.iterrows():
        color = color_map.get(row['mut_col'], 'C0')
        
        # draw the main bar from 0 to coefficient value
        if row['coef'] >= 0:
            ax.bar(i, row['coef'], width=0.6, bottom=0, 
                   color=color, edgecolor='black', linewidth=0.5)
        else:
            ax.bar(i, row['coef'], width=0.6, bottom=0,
                   color=color, edgecolor='black', linewidth=0.5)
        
        # add confidence interval whiskers
        ax.plot([i, i], [row['lo'], row['hi']], 
               color='black', linewidth=1.2, solid_capstyle='round')
        # add caps
        ax.plot([i-0.1, i+0.1], [row['lo'], row['lo']], 
               color='black', linewidth=1.2)
        ax.plot([i-0.1, i+0.1], [row['hi'], row['hi']], 
               color='black', linewidth=1.2)

    ax.axhline(0, color='0.3', ls='--', lw=1)
    ax.set_xticks(y)
    ax.set_xticklabels(group_df['mut_col'])
    ax.set_ylabel('Log-odds\n(DREAM K.O. vs. control)')
    sns.despine()
    
    return fig, ax, group_df

def plot_model_coefficients(
    model,
    include_terms=None,
    exclude_terms=("Intercept",),
    exponentiate=False,
    ci_level=0.95,
    sort_by="value",
    ascending=True,
    ax=None,
    figsize=(5, 5),
    color="C0",
    capsize=3,
    elinewidth=1.2,
    marker_size=6,
    annotate=True,
    rename_map=None,
    ):
    """
    plot coefficients of a statsmodels GLM/GLMResults as a dot-and-whisker plot with p-value annotations.

    Parameters:
    - model: fitted statsmodels model results (e.g., from fit_mutation_burden_model)
    - include_terms: list[str] or None, terms to include (if None include all except excluded)
    - exclude_terms: tuple/list[str], terms to exclude (defaults to dropping Intercept)
    - exponentiate: bool, if True show exp(coef) and exp(CI) (e.g., IRR for log-link models)
    - ci_level: float, confidence level for intervals (default 0.95)
    - sort_by: str, one of {"value", "abs", None}; how to sort rows
    - ascending: bool, sort order
    - ax: matplotlib axes or None
    - figsize: tuple, used if ax is None
    - color: str or color, marker/line color
    - capsize: float, cap size for whiskers
    - elinewidth: float, whisker line width
    - marker_size: float, marker size
    - annotate: bool, whether to annotate each point with its p-value
    - rename_map: dict or None, mapping from term name to display name

    Returns:
    - fig, ax, coef_df: the figure, axis, and DataFrame of terms with coefficients, CI, and p-values
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns

    # default rename map for common terms
    if rename_map is None:
        rename_map = {
            'age_months': 'Age\n(months)',
            'group[T.dream_ko]': 'DREAM\nK.O. = True',
            'n_effective_rf_millions': 'Read\nfamilies (M)',
            'eff_coverage_tens_of_mbs': 'Coverage\n(10Mb)',
            'n_effective_rf_millions:eff_coverage_tens_of_mbs': 'Read families\n× Coverage'
        }

    # extract parameters, confidence intervals, and p-values
    params = model.params.copy()
    pvalues = model.pvalues.copy()
    alpha = 1.0 - ci_level
    conf_int_df = model.conf_int(alpha=alpha)

    # build dataframe of coefficients
    coef_df = pd.DataFrame({
        "term": params.index,
        "coef": params.values,
        "pvalue": pvalues.values,
        "ci_lo": conf_int_df[0].values,
        "ci_hi": conf_int_df[1].values,
    })

    # filter terms
    if include_terms is not None:
        coef_df = coef_df[coef_df["term"].isin(include_terms)]
    if exclude_terms is not None and len(exclude_terms) > 0:
        coef_df = coef_df[~coef_df["term"].isin(exclude_terms)]

    # handle empty result
    if coef_df.empty:
        raise ValueError("no terms to plot after applying include/exclude filters")

    # apply exponentiation if requested (e.g., IRR for log-link)
    if exponentiate:
        coef_df = coef_df.assign(
            coef=np.exp(coef_df["coef"].astype(float)),
            ci_lo=np.exp(coef_df["ci_lo"].astype(float)),
            ci_hi=np.exp(coef_df["ci_hi"].astype(float)),
        )

    # optional renaming for readability
    if rename_map is not None:
        coef_df["display_term"] = coef_df["term"].map(lambda t: rename_map.get(t, t))
    else:
        coef_df["display_term"] = coef_df["term"]

    # sorting logic
    if sort_by == "value":
        sort_key = coef_df["coef"].values
    elif sort_by == "abs":
        sort_key = np.abs(coef_df["coef"].values)
    else:
        sort_key = np.arange(len(coef_df))
    coef_df = coef_df.iloc[np.argsort(sort_key)]
    if not ascending:
        coef_df = coef_df.iloc[::-1]
    coef_df = coef_df.reset_index(drop=True)

    # prepare axis
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # compute x-range for annotation padding
    xmin = float(coef_df["ci_lo"].min())
    xmax = float(coef_df["ci_hi"].max())
    xpad = 0.02 * (xmax - xmin) if xmax > xmin else 0.1

    # helper to format p-values nicely
    def _format_p_value(p: float) -> str:
        # prefer compact yet readable formatting
        if p < 1e-3:
            return f"p={p:.2e}"
        if p < 0.01:
            return f"p={p:.3f}"
        if p < 0.1:
            return f"p={p:.3f}"
        return f"p={p:.2f}"

    # plot a vertical reference line at 0 (or 1 if exponentiated)
    ref_value = 1.0 if exponentiate else 0.0
    ax.axvline(ref_value, color="0.3", ls="--", lw=1)

    # y positions
    y_positions = np.arange(len(coef_df))

    # draw error bars and points
    for i, row in coef_df.iterrows():
        center = float(row["coef"]) if np.isfinite(row["coef"]) else np.nan
        lo = float(row["ci_lo"]) if np.isfinite(row["ci_lo"]) else np.nan
        hi = float(row["ci_hi"]) if np.isfinite(row["ci_hi"]) else np.nan
        # xerr requires [[lower],[upper]] distances from center
        lower_dist = center - lo
        upper_dist = hi - center
        ax.errorbar(
            x=center,
            y=i,
            xerr=np.array([[lower_dist], [upper_dist]]),
            fmt="o",
            color=color,
            markerfacecolor=color,
            markeredgecolor=color,
            ecolor=color,
            elinewidth=elinewidth,
            capsize=capsize,
            markersize=5,
        )
        if annotate:
            label = _format_p_value(float(row["pvalue"]))
            y_text = i + 0.1
            ax.text(
                center,
                y_text,
                label,
                va="bottom",
                ha="center",
                fontsize=6,
                color="0.25",
            )

    # aesthetics
    ax.set_yticks(y_positions)
    ax.set_yticklabels(coef_df["display_term"].tolist())
    # irr = incidence rate ratio, which is the exponentiated coefficient from a log-linear model
    x_label = "IRR" if exponentiate else "Coefficient (log scale)"
    ax.set_xlabel(x_label)
    sns.despine()

    return fig, ax, coef_df

def detect_model_outliers(
    model,
    data=None,
    sample_col: str = 'sample',
    group_col: str = 'group',
    total_mut_col: str = 'total_mutations',
    resid_threshold: float = 3.0,
    cooks_factor: float = 4.0,
    leverage_factor: float = 2.0,
    zscore_threshold: float = 3.0,
):
    """
    identify potential outliers for a fitted GLM model using multiple diagnostics:
    1) absolute deviance residuals > resid_threshold (default ~3)
    2) cook's distance > cooks_factor / n (default 4/n)
    3) high leverage: leverage > leverage_factor * mean_leverage (default 2 * p/n)
    4) z-score of total mutations relative to all samples > zscore_threshold
    5) z-score of total mutations relative to within-group samples > zscore_threshold

    parameters
    -----------
    model : statsmodels GLMResults
        fitted model returned by fit_mutation_burden_model
    data : pd.DataFrame or None
        optional dataframe used to fit the model; if provided and it contains
        a column named by sample_col, that column will be attached for readability
    sample_col : str
        column name in data with sample identifiers (default 'sample')
    group_col : str
        column name in data with group identifiers (default 'group')
    total_mut_col : str
        column name in data with total mutation counts (default 'total_mutations')
    resid_threshold : float
        threshold for absolute deviance residuals
    cooks_factor : float
        factor for cook's distance threshold (factor / n)
    leverage_factor : float
        factor relative to mean leverage used as the high leverage threshold
    zscore_threshold : float
        threshold for absolute z-scores of total mutations

    returns
    -------
    outlier_df : pd.DataFrame
        dataframe indexed by observation with columns:
        ['deviance_resid', 'abs_deviance_resid', 'cooks_d', 'leverage',
         'n', 'p_estimated', 'mean_leverage', 'threshold_abs_deviance',
         'threshold_cooks_d', 'threshold_leverage', 'flag_deviance',
         'flag_cooks_d', 'flag_leverage', 'total_mutations', 'zscore_all',
         'abs_zscore_all', 'zscore_within_group', 'abs_zscore_within_group',
         'flag_zscore_all', 'flag_zscore_within_group', 'is_outlier'] 
         and optionally 'sample' and 'group'
    """
    import numpy as np
    import pandas as pd

    # compute influence diagnostics from the fitted model
    influence = model.get_influence(observed=True)

    # deviance residuals
    try:
        deviance_residuals = np.asarray(model.resid_deviance)
    except Exception:
        # fallback if property not present
        deviance_residuals = np.asarray(getattr(influence, 'resid_deviance'))

    # cook's distance may be returned as (values, pvals)
    cooks_raw = influence.cooks_distance
    cooks_d = np.asarray(cooks_raw[0] if isinstance(cooks_raw, tuple) else cooks_raw)

    # leverage (hat) values
    leverage = np.asarray(influence.hat_matrix_diag)

    # sizes and thresholds
    n_obs = int(getattr(model, 'nobs', len(leverage)))
    # p can be estimated as the trace of the hat matrix
    p_estimated = float(np.sum(leverage))
    mean_leverage = float(np.mean(leverage)) if n_obs > 0 else float('nan')
    cooks_threshold = float(cooks_factor) / float(n_obs) if n_obs > 0 else float('nan')
    leverage_threshold = float(leverage_factor) * mean_leverage if np.isfinite(mean_leverage) else float('nan')

    # choose an index that mirrors the original data index when possible
    if hasattr(model, 'fittedvalues') and hasattr(model.fittedvalues, 'index'):
        index = model.fittedvalues.index
    elif hasattr(model, 'model') and hasattr(model.model, 'data') and hasattr(model.model.data, 'row_labels'):
        index = model.model.data.row_labels
    else:
        index = pd.Index(range(n_obs), name='row')

    outlier_df = pd.DataFrame(
        {
            'deviance_resid': deviance_residuals,
            'abs_deviance_resid': np.abs(deviance_residuals),
            'cooks_d': cooks_d,
            'leverage': leverage,
        },
        index=index,
    )

    # attach thresholds and flags
    outlier_df = outlier_df.assign(
        n=n_obs,
        p_estimated=p_estimated,
        mean_leverage=mean_leverage,
        threshold_abs_deviance=float(resid_threshold),
        threshold_cooks_d=cooks_threshold,
        threshold_leverage=leverage_threshold,
    )

    outlier_df['flag_deviance'] = outlier_df['abs_deviance_resid'] > outlier_df['threshold_abs_deviance']
    outlier_df['flag_cooks_d'] = outlier_df['cooks_d'] > outlier_df['threshold_cooks_d']
    outlier_df['flag_leverage'] = outlier_df['leverage'] > outlier_df['threshold_leverage']
    
    # initialize z-score columns with NaN
    outlier_df['total_mutations'] = np.nan
    outlier_df['zscore_all'] = np.nan
    outlier_df['abs_zscore_all'] = np.nan
    outlier_df['zscore_within_group'] = np.nan
    outlier_df['abs_zscore_within_group'] = np.nan
    outlier_df['flag_zscore_all'] = False
    outlier_df['flag_zscore_within_group'] = False

    # optionally add sample identifiers and compute z-scores
    if data is not None:
        try:
            # align data with outlier_df
            aligned_data = None
            if data.index.equals(outlier_df.index):
                aligned_data = data
            elif len(data) == len(outlier_df):
                # position-wise align if lengths match but indexes differ
                aligned_data = data.reset_index(drop=True)
                aligned_data.index = outlier_df.index
            
            if aligned_data is not None:
                # add sample and group identifiers if available
                if sample_col in aligned_data.columns:
                    outlier_df['sample'] = aligned_data[sample_col]
                if group_col in aligned_data.columns:
                    outlier_df['group'] = aligned_data[group_col]
                
                # compute z-scores for total mutations if available
                if total_mut_col in aligned_data.columns:
                    outlier_df['total_mutations'] = aligned_data[total_mut_col]
                    total_mut_values = aligned_data[total_mut_col].values
                    
                    # z-score relative to all other samples (leave-one-out)
                    for i, (idx, row) in enumerate(outlier_df.iterrows()):
                        # exclude current sample from calculation
                        other_values = np.concatenate([total_mut_values[:i], total_mut_values[i+1:]])
                        
                        if len(other_values) > 1:  # need at least 2 other samples
                            other_mean = float(np.mean(other_values))
                            other_std = float(np.std(other_values, ddof=1))  # sample std
                            
                            if other_std > 0:
                                zscore = (row['total_mutations'] - other_mean) / other_std
                                outlier_df.loc[idx, 'zscore_all'] = zscore
                                outlier_df.loc[idx, 'abs_zscore_all'] = abs(zscore)
                                outlier_df.loc[idx, 'flag_zscore_all'] = abs(zscore) > zscore_threshold
                    
                    # z-score relative to within-group samples (leave-one-out)
                    if group_col in aligned_data.columns:
                        for i, (idx, row) in enumerate(outlier_df.iterrows()):
                            current_group = row.get('group')
                            if pd.isna(current_group):
                                continue
                                
                            # get indices of samples in the same group, excluding current sample
                            same_group_mask = (aligned_data[group_col] == current_group).values
                            same_group_indices = np.where(same_group_mask)[0]
                            other_group_indices = same_group_indices[same_group_indices != i]
                            
                            if len(other_group_indices) > 1:  # need at least 2 other samples in group
                                other_group_values = total_mut_values[other_group_indices]
                                other_group_mean = float(np.mean(other_group_values))
                                other_group_std = float(np.std(other_group_values, ddof=1))  # sample std
                                
                                if other_group_std > 0:
                                    group_zscore = (row['total_mutations'] - other_group_mean) / other_group_std
                                    outlier_df.loc[idx, 'zscore_within_group'] = group_zscore
                                    outlier_df.loc[idx, 'abs_zscore_within_group'] = abs(group_zscore)
                                    outlier_df.loc[idx, 'flag_zscore_within_group'] = abs(group_zscore) > zscore_threshold
                        
        except Exception:
            # stay robust if alignment fails
            pass

    # update overall outlier flag to include z-score flags
    flag_cols = ['flag_deviance', 'flag_cooks_d', 'flag_leverage', 'flag_zscore_all', 'flag_zscore_within_group']
    outlier_df['is_outlier'] = outlier_df[flag_cols].any(axis=1)

    return outlier_df

