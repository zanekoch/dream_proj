import pandas as pd
import seaborn as sns
import statsmodels.formula.api as smf
import statsmodels.api as sm
import numpy as np

def make_regr_df(
    dream_expression,
    cell_types_to_drop = [],
    ):
    regression_df = dream_expression.copy(deep=True)
    
    
    # counts of expression
    regression_df['log_n_counts'] = np.log2(regression_df['n_counts'])
    regression_df['z_log_n_counts'] = (regression_df['log_n_counts'] - regression_df['log_n_counts'].mean()) / regression_df['log_n_counts'].std()
    
    # neurotypical reference
    regression_df['neurotypical_reference'] = regression_df['Neurotypical reference']
    regression_df['neurotypical_reference'] = regression_df['Neurotypical reference'].map({
        'True': 1,
        'False': 0
    })
    
    # ADNC, map from 'High', 'Intermediate', 'Low', 'Not AD', 'Reference' to 4,3,2,1,0
    regression_df['adnc'] = regression_df['Overall AD neuropathological Change'].map({
        'High': 4,
        'Intermediate': 3,
        'Low': 2,
        'Not AD': 1,
        'Reference': 0
    }).astype(int)
    regression_df['adnc234'] = regression_df['adnc'] >= 2 # true is bad
    
    # Thal score, map from 'Thal 4', 'Thal 5', 'Thal 3', 'Thal 0', 'Thal 2', 'Reference', 'Thal 1' to 4, 5, 3, 0, 2, 0, 1
    regression_df['thal'] = regression_df['Thal'].map({
        'Thal 4': 4,
        'Thal 5': 5,
        'Thal 3': 3,
        'Thal 0': 0,
        'Thal 2': 2,
        'Reference': 0,
        'Thal 1': 1
    }).astype(int)
    regression_df['thal345'] = regression_df['thal'] >= 3 # true is bad
    
    # NIA reagan score
    # Map 'CERAD score' from 'Moderate', 'Frequent', 'Absent', 'Sparse', 'Reference' to 1,0,3,2,3
    regression_df['ceradsc'] = regression_df['CERAD score'].map({
        'Moderate': 1,
        'Frequent': 0,
        'Absent': 3,
        'Sparse': 2,
        'Reference': 3
    })
    regression_df['cerad12'] = regression_df['ceradsc'] <= 1 # true is bad
    
    # braak stage 
    # Map 'Braak stage' from 'Reference', 'Braak I', 'Braak II', 'Braak III', 'Braak IV', 'Braak V', 'Braak VI' to 0,1,2,3,4,5,6
    regression_df['braaksc'] = regression_df['Braak'].map({
        'Reference': 0,
        'Braak I': 1,
        'Braak II': 2,
        'Braak III': 3,
        'Braak IV': 4,
        'Braak V': 5,
        'Braak VI': 6
    })
    regression_df['braak456'] = regression_df['braaksc'] >= 4 # true is bad
    
    # nia reagan score
    #regression_df['nia_reagan12'] = regression_df['nia_reagan'] <= 2 # true is bad
    
    # congitive impairment
    #regression_df.query("cogdx != 6", inplace=True)
    #regression_df['cogdx45'] = regression_df['cogdx'] >= 4 # true is bad
    
    # apoe4 genotype
    regression_df['apoe4'] = regression_df['APOE Genotype'].astype(str).str.contains('4') == True # true is bad
    regression_df['apoe44'] = regression_df['APOE Genotype'] == '4/4' # true is bad
    # post mortem interval
    # regression_df['pmi_rescaled'] = regression_df['pmi'] / 10
    
    # age of death 
    regression_df['age_death'] = regression_df['age']
    regression_df['age_death_rescaled'] = regression_df['age'] / 10
    regression_df['z_age_death_rescaled'] = (regression_df['age_death_rescaled'] - regression_df['age_death_rescaled'].mean()) / regression_df['age_death_rescaled'].std()

    # cell level variables
    # number of captured genes
    regression_df['log_n_genes'] = np.log2(regression_df['n_genes'])
    regression_df['z_log_n_genes'] = (regression_df['log_n_genes'] - regression_df['log_n_genes'].mean()) / regression_df['log_n_genes'].std()
    
    # reads per gene
    regression_df['reads_per_gene'] = regression_df['log_n_counts'] / regression_df['log_n_genes']
    regression_df['z_reads_per_gene'] = (regression_df['reads_per_gene'] - regression_df['reads_per_gene'].mean()) / regression_df['reads_per_gene'].std()
    
    # z score dream
    regression_df['z_DREAM_normalized_enrichment_score_resid'] = (regression_df['DREAM_normalized_enrichment_score_resid'] - regression_df['DREAM_normalized_enrichment_score_resid'].mean()) / regression_df['DREAM_normalized_enrichment_score_resid'].std()

    # drop cell types
    regression_df.query("cell_type not in @cell_types_to_drop", inplace=True)
    
    # z score mutation burden
    regression_df['z_sc_mutation_burden_in_bulk_per_genome'] = (regression_df['sc_mutation_burden_in_bulk_per_genome'] - regression_df['sc_mutation_burden_in_bulk_per_genome'].mean()) / regression_df['sc_mutation_burden_in_bulk_per_genome'].std()
    
    regression_df['log_DREAM_normalized_enrichment_score_resid'] = np.log2(regression_df['DREAM_normalized_enrichment_score_resid'])
    regression_df['z_log_DREAM_normalized_enrichment_score_resid'] = (regression_df['log_DREAM_normalized_enrichment_score_resid'] - regression_df['log_DREAM_normalized_enrichment_score_resid'].mean()) / regression_df['log_DREAM_normalized_enrichment_score_resid'].std()
    # transform mutation burden
    regression_df['log_sc_mutation_burden_in_bulk_per_kb'] = np.log2(regression_df['sc_mutation_burden_in_bulk_per_kb'])
    regression_df['z_log_sc_mutation_burden_in_bulk_per_kb'] = (regression_df['log_sc_mutation_burden_in_bulk_per_kb'] - regression_df['log_sc_mutation_burden_in_bulk_per_kb'].mean()) / regression_df['log_sc_mutation_burden_in_bulk_per_kb'].std()
    # transform mutation burden per genome
    regression_df['log_sc_mutation_burden_in_bulk_per_genome'] = np.log2(regression_df['sc_mutation_burden_in_bulk_per_genome'])
    regression_df['z_log_sc_mutation_burden_in_bulk_per_genome'] = (regression_df['log_sc_mutation_burden_in_bulk_per_genome'] - regression_df['log_sc_mutation_burden_in_bulk_per_genome'].mean()) / regression_df['log_sc_mutation_burden_in_bulk_per_genome'].std()
        
    regression_df['log_sc_mutation_burden_in_bulk_per_genome_nanfilled'] = np.log2(1+regression_df['sc_mutation_burden_in_bulk_per_genome_nanfilled'])
    regression_df['z_log_sc_mutation_burden_in_bulk_per_genome_nanfilled'] = (regression_df['log_sc_mutation_burden_in_bulk_per_genome_nanfilled'] - regression_df['log_sc_mutation_burden_in_bulk_per_genome_nanfilled'].mean()) / regression_df['log_sc_mutation_burden_in_bulk_per_genome_nanfilled'].std()
        
    # mask map 'nan' to np.nan and convert to float
    regression_df['MMSE'] = regression_df['Last MMSE Score'].map(lambda x: np.nan if x == 'nan' or x == 'Reference' or x == '' else float(x))
    regression_df['MOCA'] = regression_df['Last MOCA Score'].map(lambda x: np.nan if x == 'nan' or x == 'Reference' or x == '' else float(x))
    regression_df['CASI'] = regression_df['Last CASI Score'].map(lambda x: np.nan if x == 'nan' or x == 'Reference' or x == '' else float(x))
    
    return regression_df

def plot_coefficients(
    ax, 
    model = None,
    coef_df=pd.DataFrame(),
    xlim=None,
    colors = {True: 'black', False: 'grey'},
    use_std_dev = True
    ):
    # check if coef_Df is empty
    if coef_df.shape[0] == 0:
        # Step 1: Extract the coefficients and standard errors
        coefficients = model.params[1:]
        std_errors = model.bse[1:]

        # Calculate the confidence intervals (95%)
        conf_int = model.conf_int()
        conf_int.columns = ['Lower CI', 'Upper CI']

        # Combine coefficients, standard errors, and confidence intervals into a DataFrame
        coef_df = pd.DataFrame({
            'Coefficient': coefficients,
            'Std Error': std_errors,
            'Lower CI': conf_int['Lower CI'],
            'Upper CI': conf_int['Upper CI']
        })
        
        coef_df['Variable'] = coef_df.index
        # remove last row
        coef_df = coef_df.loc[coef_df.index[coef_df.index != 'Intercept']]

    # set as sig if CI does not include 0
    coef_df['sig'] = (coef_df['Lower CI'] > 0) | (coef_df['Upper CI'] < 0)
    # Step 2: Sort by coefficient value for better visualization
    coef_df = coef_df.sort_values(by='Coefficient', ascending=True)

    # remove any coefficients that have "tissue", "celltype", or "Subclass" in the name
    coef_df = coef_df.loc[~coef_df['Variable'].str.contains('tissue|cell_type|Subclass')]


    # Create the plot
    sns.pointplot(
        x="Coefficient", y="Variable", data=coef_df,
        linestyle='none',  hue = 'sig', palette = colors,
        ax = ax, legend = False
        )

    # Add error bars manually, with caps
    if use_std_dev:
        multiplier = 2
    else:
        multiplier = 1
    for index, row in coef_df.iterrows():
        if row['sig']:
            ax.errorbar(x=row['Coefficient'], y=row['Variable'], xerr=multiplier*row['Std Error'], fmt='none', ecolor=colors[True], capsize=2.5)
        else :
            ax.errorbar(x=row['Coefficient'], y=row['Variable'], xerr=multiplier*row['Std Error'], fmt='none', ecolor=colors[False], capsize=2.5)

    # Add a vertical line at 0 to denote no effect
    ax.axvline(0, color='black', linestyle='--', linewidth = 0.8)

    # Titles and labels
    ax.set_xlabel('Coefficient Value')
    ax.set_ylabel('Variable')
    sns.despine()
    if xlim:
        ax.set_xlim(xlim)
    
def do_regr_one_by_one(
    regression_df, 
    target: str,
    variables: list,
    count_vars = [],
    ext_covars = [],
    offset = True, 
    offset_str = 'log_n_counts',
    correct_hierarchy = False,
    family = sm.families.NegativeBinomial(alpha=1),
    ):
    coefs, ses, pvals, ci_highs, ci_lows = [], [], [], [], []
    for var in variables:
        count_str = ' + '.join(count_vars)
        ext_covar_str = ' + '.join(ext_covars)
        formula = f'{target} ~ {count_str} + {ext_covar_str}'
        formula += f' + {var}'
        # fit the model
        if offset:
            if correct_hierarchy:
                # for all cols, drop nans
                cols = count_vars + ext_covars + [var] + [target] + ['donor_id'] + [offset_str]
                regression_df = regression_df.dropna(subset=cols)
                model = smf.glm(formula=formula, data=regression_df, 
                                family=family, 
                                offset=regression_df[offset_str]
                                ).fit(cov_type='cluster', cov_kwds={'groups': regression_df['donor_id'].astype(str)})
            else:
                model = smf.glm(formula=formula, data=regression_df, 
                                family=family, 
                                offset=regression_df[offset_str]
                                ).fit()
        else:
            if correct_hierarchy:
                # for all cols, drop nans
                cols = count_vars + ext_covars + [var] + [target] + ['donor_id'] + [offset_str]
                regression_df = regression_df.dropna(subset=cols)
                model = smf.glm(formula=formula, data=regression_df, 
                                family=family, 
                                ).fit(cov_type='cluster', cov_kwds={'groups': regression_df['donor_id'].astype(str)})
            else:
                model = smf.glm(formula=formula, data=regression_df, 
                                family=family, 
                                ).fit()
        # check if the var string contains a number, if so its an indicator variable and we need to add [T.True]
        if any(char.isdigit() for char in var):
            var = var + '[T.True]'
        elif var == 'neurotypical_reference':
            var = var + '[T.1]'
        # get the coefficients, standard errors, p-value, and confidence interval for the var
        coef = model.params[var]
        se = model.bse[var]
        pval = model.pvalues[var]
        ci_high = model.conf_int().loc[var,1]
        ci_low = model.conf_int().loc[var,0]
        # append to the lists
        coefs.append(coef)
        ses.append(se)
        pvals.append(pval)
        ci_highs.append(ci_high)
        ci_lows.append(ci_low)
    # create df, matching format of plot_coefficients()
    coef_df = pd.DataFrame({
        'Coefficient': coefs,
        'Std Error': ses,
        'p-value': pvals,
        'Lower CI': ci_lows,
        'Upper CI': ci_highs,
        'Variable': variables
    })
    return coef_df

def do_regr_all_together(
    regression_df, 
    target: str,
    variables: list,
    count_vars: list,
    ext_covars: list,
    offset=True,
    offset_str='log_n_counts',
    correct_hierarchy=False,
    family = sm.families.NegativeBinomial(alpha=.5),
    alpha = 0.01
    ):
    var_str = ' + '.join(variables)
    count_str = ' + '.join(count_vars)
    ext_covar_str = ' + '.join(ext_covars)
    formula = f'{target} ~ {count_str} + {ext_covar_str} + {var_str}'
    
    # for each variable, count_Var, and ext_covar, check if there are any nans in the col + target
    num_rows = regression_df.shape[0]
    regression_df = regression_df.dropna(subset=[target] + variables + count_vars + ext_covars + ['donor_id'])
    num_rows_after = regression_df.shape[0]
    print(f'Dropped {num_rows - num_rows_after} rows with nans')
    
    if offset:
        # for all cols, drop nans
        cols = count_vars + ext_covars + [target] + variables + ['donor_id'] + [offset_str]
        regression_df = regression_df.dropna(subset=cols)
        print("Number of cells after dropping nans: ", regression_df.shape[0])
        if correct_hierarchy:
            model = smf.glm(formula=formula, data=regression_df, 
                            family=family, 
                            offset=regression_df[offset_str], 
                            ).fit(cov_type='cluster', cov_kwds={'groups': regression_df['donor_id'].astype(str)})
        else:
            model = smf.glm(formula=formula, data=regression_df, 
                            family=family, 
                            offset=regression_df[offset_str], 
                            ).fit()
    else:
        # for all cols, drop nans
        cols = count_vars + ext_covars + [target] + ['donor_id']
        regression_df = regression_df.dropna(subset=cols)
        print("Number of cells after dropping nans: ", regression_df.shape[0])
        if correct_hierarchy:
            model = smf.glm(formula=formula, data=regression_df, 
                            family=family, 
                            ).fit(cov_type='cluster', cov_kwds={'groups': regression_df['donor_id'].astype(str)})
        else:
            model = smf.glm(formula=formula, data=regression_df, 
                            family=family, 
                            ).fit()
    
    return model

def create_nia_reagan_score(dream_expression):
    dream_expression['nia_reagan'] = 4
    dream_expression['nia_reagan'] = np.where(
        (dream_expression['ceradsc'] == 4),
        4,
        dream_expression['nia_reagan']
    )
    dream_expression['nia_reagan'] = np.where(
        (dream_expression['braaksc'] <= 2),
        4,
        dream_expression['nia_reagan']
    )
    dream_expression['nia_reagan'] = np.where(
        (dream_expression['ceradsc'] == 3) & (dream_expression['braaksc'] >= 3),
        3,
        dream_expression['nia_reagan']
    )
    dream_expression['nia_reagan'] = np.where(
        (dream_expression['ceradsc'] == 2) & (dream_expression['braaksc'] >= 3) & (dream_expression['braaksc'] <= 4),
        2,
        dream_expression['nia_reagan']
    )
    dream_expression['nia_reagan'] = np.where(
        (dream_expression['ceradsc'] == 2) & (dream_expression['braaksc'] >= 5),
        1,
        dream_expression['nia_reagan']
    )
    dream_expression['nia_reagan'] = np.where(
        (dream_expression['ceradsc'] == 1) & (dream_expression['braaksc'] <= 2),
        2,
        dream_expression['nia_reagan']
    )
    dream_expression['nia_reagan'] = np.where(
        (dream_expression['ceradsc'] == 1) & (dream_expression['braaksc'] >= 3),
        1,
        dream_expression['nia_reagan']
    )

