import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import openpyxl

PBRM1_data = pd.read_excel('~/icr/simko/data/pbrm1_experiment_data/all_pbrm1_results.xlsx')

PBRM1_data 

#selecting the CEN proteins

CEN_proteins = ['MIS18a','MIS18BP1','RSF1','KAT7','SMC2','SMC4','NCAPD3','NCAPG2','NCAPH2',
                'CENPC','CENPT','CENPS','CENPX','CENPO','CENPU','CENPL','CENPN','CENPH','CENPI','CENPK',
                'KNL1','ZWINT','NDC80','NUF2','SPC24','SPC25','DSN1','MIS12','NSL1','PMF1',
                'INCENP','BIRC5','CDCA8','AURKB',
                'CENPE','CENPF','CENPB','CENPV','CENPJ','CDCA2',
                'CBX1','SIRT6','UHRF2','SUV39H1','SETDB1','DAXX','SSRP1','SMARCAD1',
                'ATRX','BAZ1A','BAZ1B','CBX3','CBX5','CHRAC1','DNMT1','EZH2','HELLS','KDM4A','KMT5C','LRWD1','MACROH2A1','POLE3','SMARCA5']

##Volcano plot - for experimental PBRM1 knockout data set
#changing colnames so its easier to work with
PBRM1_data.columns = PBRM1_data.columns.str.replace(' ', '_').str.replace('-', '').str.replace('/', '_')

PBRM1_data.columns

# specifying col name sto make it easier to read
logfc = 'mean_log2_KO_WT'
pval = 'Log_Ttest_pvalue_'
proteins = 'Gene_Names_(primary)'

PBRM1_data['color'] = 'black'  # Default color for all proteins
PBRM1_data.loc[PBRM1_data['Gene_Names_(primary)'] == 'PBRM1', 'color'] = 'green'  # PBRM1 in green
PBRM1_data.loc[PBRM1_data['Gene_Names_(primary)'].isin(CEN_proteins), 'color'] = 'purple'  # CEN proteins in purple

#basic plot without colour highlights
plt.figure(figsize=(10, 8))
# Plot black points (default)
sns.scatterplot(
    data=PBRM1_data[PBRM1_data['color'] == 'black'],  # Only plot black points
    x=logfc,
    y=pval,
    color='black',
    alpha=1,
    edgecolor=None,
    s=50  # Adjust size of the points
)

# Plot green points (PBRM1)
sns.scatterplot(
    data=PBRM1_data[PBRM1_data['color'] == 'green'],  # Only plot green points
    x=logfc,
    y=pval,
    color='green',
    alpha=1,
    edgecolor=None,
    s=50  # Slightly larger size to emphasize
)

# Plot purple points (CEN_proteins)
sns.scatterplot(
    data=PBRM1_data[PBRM1_data['color'] == 'purple'],  # Only plot purple points
    x=logfc,
    y=pval,
    color='purple',
    alpha=1,  # Make purple points slightly more opaque
    edgecolor=None,
    s=50  # Adjust size
)

# Annotate PBRM1
pbrm1_row = PBRM1_data[PBRM1_data['Gene_Names_(primary)'] == 'PBRM1']
if not pbrm1_row.empty:
    plt.text(
        pbrm1_row['mean_log2_KO_WT'].values[0],
        pbrm1_row['Log_Ttest_pvalue_'].values[0],
        'PBRM1',
        color='green',
        fontsize=12,
        fontweight='bold'
    )

# Add plot labels and lines
plt.title('Volcano Plot: Highlighted Proteins', fontsize=16)
plt.xlabel('Log2 Fold Change (KO/WT)', fontsize=14)
plt.ylabel('-Log10 p-value', fontsize=14)
plt.axhline(y=1.3, color='gray', linestyle='--', linewidth=0.8)  # Threshold for p-value significance
plt.axvline(x=0, color='black', linestyle='--', linewidth=0.8)   # LFC = 0 line

# Customize legend
handles = [
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=8, label='Others'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='green', markersize=8, label='PBRM1'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='purple', markersize=8, label='CEN Proteins')
]
plt.legend(handles=handles, loc='upper right', fontsize=12)

plt.tight_layout()
plt.show()








#now look at the proteins up and down regulated from PBRM1 simulated knockout
abundance = pd.read_csv('~/icr/simko/data/simko2_data/passport_prots.csv', index_col=0)
abundance.index = abundance.index.astype(str)


class_df = get_classes_by_mean_abundance(ko_proteins=['PBRM1'], abundance=abundance, n=30)
control_diffs = get_control_differentials(abundance, class_df, k=30)
diffs = get_ko_differentials(abundance=abundance, class_df=class_df)
diffs
PBRM1_row = diffs[diffs["protein"] == "PBRM1"]
PBRM1_row
diff_stats = get_significance(diffs, control_diffs, n=30)
diff_stats
#-np.log10() on pvals for volcano plot
diff_stats['log10_pval']= -np.log10(diff_stats['adjusted_p'])
diff_stats = diff_stats.sort_values(by='log10_pval', ascending = False)

diff_stats[:20]


# volcano plot for sim data
#specifying col names
logfcsim = 'diff'
pvalsim = 'log10_pval'
protein = 'protein'

#colours for different points
diff_stats['color'] = 'black'  # Default color for all proteins
diff_stats.loc[diff_stats['protein'] == 'PBRM1', 'color'] = 'green'  # PBRM1 in green
diff_stats.loc[diff_stats['protein'].isin(CEN_proteins), 'color'] = 'purple'  # CEN proteins in purple


#basic plot without colour highlights
plt.figure(figsize=(10, 8))
# Plot black points (default)
sns.scatterplot(
    data=diff_stats[diff_stats['color'] == 'black'],  # Only plot black points
    x=logfcsim,
    y=pvalsim,
    color='black',
    alpha=1,
    edgecolor=None,
    s=50  
)

# Plot green points (PBRM1)
sns.scatterplot(
    data=diff_stats[diff_stats['color'] == 'green'],  # Only plot green points
    x=logfcsim,
    y=pvalsim,
    color='green',
    alpha=1,
    edgecolor=None,
    s=50  
)

# Plot purple points (CEN_proteins)
sns.scatterplot(
    data=diff_stats[diff_stats['color'] == 'purple'],  # Only plot purple points
    x=logfcsim,
    y=pvalsim,
    color='purple',
    alpha=1,  
    edgecolor=None,
    s=50  
)

# Add plot labels and lines
plt.title('Volcano Plot: Highlighted Proteins', fontsize=16)
plt.xlabel('Log2 Fold Change (KO/WT)', fontsize=14)
plt.ylabel('-Log10 p-value', fontsize=14)
plt.axhline(y=1.3, color='gray', linestyle='--', linewidth=0.8)  # Threshold for p-value significance
plt.axvline(x=0, color='black', linestyle='--', linewidth=0.8)   # LFC = 0 line

# Customize legend
handles = [
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=8, label='Others'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='green', markersize=8, label='PBRM1'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='purple', markersize=8, label='CEN Proteins')
]
plt.legend(handles=handles, loc='upper right', fontsize=12)

plt.tight_layout()
plt.show()




#list of CEN proteins that are upregulated in the simulated data
sim_upreg = diff_stats.loc[diff_stats['protein'].isin(CEN_proteins) 
                           & (diff_stats['diff'] >= 0), 'protein']


#list of cen proteins that are upregulated in PBRM1 data
exp_upreg = PBRM1_data.loc[PBRM1_data['Gene_Names_(primary)'].isin(CEN_proteins)
                            & (PBRM1_data['mean_log2_KO_WT'] >=0), 'Gene_Names_(primary)' ]

len(sim_upreg)
len(exp_upreg)
common_upreg_proteins = set(sim_upreg).intersection(set(exp_upreg))
common_upreg_proteins


#list of CEN proteins thats are downregulated in both
#list of CEN proteins that are upregulated in the simulated data
sim_downreg = diff_stats.loc[diff_stats['protein'].isin(CEN_proteins) 
                           & (diff_stats['diff'] <= 0), 'protein']

#list of cen proteins that are upregulated in PBRM1 data
exp_downreg = PBRM1_data.loc[PBRM1_data['Gene_Names_(primary)'].isin(CEN_proteins)
                            & (PBRM1_data['mean_log2_KO_WT'] <=0), 'Gene_Names_(primary)' ]
common_downreg_proteins = set(sim_downreg).intersection(set(exp_downreg))
common_downreg_proteins

len(common_downreg_proteins)
len(common_upreg_proteins)

#comparing how many up or down regulate proteins there are in sim and exp data
#at least in both, there are less upregulated cen proteins
sim_downreg[20:34]
sim_upreg
exp_downreg[1:25]
exp_upreg


## comparing the differences between the experimental dn simulated data
#create a data frame contaning a column for proteins, difs from sim and diffs from experimental
#need to make sure they both contain the same protein

significant_proteins = diff_stats.loc[diff_stats['adjusted_p'] < 0.05, 'protein']
significant_proteins

common_proteins = set(significant_proteins).intersection(set(PBRM1_data['Gene_Names_(primary)'])) #7106 proteins

# Extract the "diff" values for the common proteins in both datasets
simulated_diff = diff_stats.loc[diff_stats['protein'].isin(common_proteins), 
                                ['protein', 'diff']].set_index('protein')
pbrm1_diff = PBRM1_data.loc[PBRM1_data['Gene_Names_(primary)'].isin(common_proteins),
                             ['Gene_Names_(primary)', 'mean_log2_KO_WT']]
pbrm1_diff = pbrm1_diff.rename(columns={'Gene_Names_(primary)': 'protein', 'mean_log2_KO_WT': 'diff'}).set_index('protein')

# Combine into a single DataFrame for comparison
comparison_df = simulated_diff.join(pbrm1_diff, lsuffix='_simulated', rsuffix='_pbrm1')

comparison_df

plt.figure(figsize=(8, 6))
sns.regplot(
    x=comparison_df['diff_simulated'], 
    y=comparison_df['diff_pbrm1'], 
    ci=None,  # Disable confidence intervals for clarity
    line_kws={'color': 'red', 'linewidth': 1.5}  # Styling for the line of best fit
)
# Add plot labels
plt.title('Comparison of Diff Values Between Simulated and PBRM1 Data')
plt.xlabel('Diff (Simulated Data)')
plt.ylabel('Diff (PBRM1 Data)')
plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
plt.axvline(0, color='gray', linestyle='--', linewidth=0.8)
plt.grid(True)

# Show the plot
plt.show()

from scipy.stats import pearsonr
correlation, p_value = pearsonr(comparison_df['diff_simulated'], comparison_df['diff_pbrm1'])
correlation
p_value

#scatter plot of only the cen proteins
cen_comparison_df = comparison_df.loc[comparison_df.index.isin(CEN_proteins)]

# Scatter plot
plt.figure(figsize=(8, 6))
sns.regplot(
    x=cen_comparison_df['diff_simulated'], 
    y=cen_comparison_df['diff_pbrm1'],
    ci=None,
    line_kws={'color': 'red', 'linewidth': 1.5}  # Style for the line of best fit
)

# Add plot labels
plt.title('Comparison of Diff Values for CEN Proteins')
plt.xlabel('Diff (Simulated Data)')
plt.ylabel('Diff (PBRM1 Data)')
plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
plt.axvline(0, color='gray', linestyle='--', linewidth=0.8)
plt.grid(True)

# Show the plot
plt.show()

correlation_cen, p_value_cen = pearsonr(cen_comparison_df['diff_simulated'], cen_comparison_df['diff_pbrm1'])
correlation_cen
p_value_cen








#now looking at ccle data
abundance_ccle = pd.read_csv('~/icr/simko/simko_git/simko2_analysis-main/data/simko2_data/CCLE_all.csv', index_col=0)
abundance_ccle.index = abundance_ccle.index.astype(str)


class_df = get_classes_by_mean_abundance(ko_proteins=['PBRM1'], abundance=abundance_ccle, n=30)
control_diffs = get_control_differentials(abundance_ccle, class_df, k=30)
diffs = get_ko_differentials(abundance=abundance_ccle, class_df=class_df)
diffs
PBRM1_row = diffs[diffs["protein"] == "PBRM1"]
PBRM1_row
diff_stats = get_significance(diffs, control_diffs, n=30)
diff_stats
#-np.log10() on pvals for volcano plot
diff_stats['log10_pval']= -np.log10(diff_stats['adjusted_p'])
diff_stats = diff_stats.sort_values(by='log10_pval', ascending = False)

diff_stats[:20]

significant_proteins = diff_stats.loc[diff_stats['adjusted_p'] < 0.05, 'protein']
significant_proteins

common_proteins = set(significant_proteins).intersection(set(PBRM1_data['Gene_Names_(primary)'])) #7106 proteins

# Extract the "diff" values for the common proteins in both datasets
simulated_diff = diff_stats.loc[diff_stats['protein'].isin(common_proteins), 
                                ['protein', 'diff']].set_index('protein')
pbrm1_diff = PBRM1_data.loc[PBRM1_data['Gene_Names_(primary)'].isin(common_proteins),
                             ['Gene_Names_(primary)', 'mean_log2_KO_WT']]
pbrm1_diff = pbrm1_diff.rename(columns={'Gene_Names_(primary)': 'protein', 'mean_log2_KO_WT': 'diff'}).set_index('protein')

# Combine into a single DataFrame for comparison
comparison_df = simulated_diff.join(pbrm1_diff, lsuffix='_simulated', rsuffix='_pbrm1')

comparison_df

plt.figure(figsize=(8, 6))
sns.regplot(
    x=comparison_df['diff_simulated'], 
    y=comparison_df['diff_pbrm1'], 
    ci=None,  # Disable confidence intervals for clarity
    line_kws={'color': 'red', 'linewidth': 1.5}  # Styling for the line of best fit
)
# Add plot labels
plt.title('Comparison of Diff Values Between Simulated and PBRM1 Data')
plt.xlabel('Diff (Simulated Data)')
plt.ylabel('Diff (PBRM1 Data)')
plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
plt.axvline(0, color='gray', linestyle='--', linewidth=0.8)
plt.grid(True)

# Show the plot
plt.show()

correlation, p_value = pearsonr(comparison_df['diff_simulated'], comparison_df['diff_pbrm1'])
correlation
p_value

#only cen proteins from ccle data
#scatter plot of only the cen proteins
cen_comparison_df = comparison_df.loc[comparison_df.index.isin(CEN_proteins)]

# Scatter plot
plt.figure(figsize=(8, 6))
sns.regplot(
    x=cen_comparison_df['diff_simulated'], 
    y=cen_comparison_df['diff_pbrm1'],
    ci=None,
    line_kws={'color': 'red', 'linewidth': 1.5}  # Style for the line of best fit
)

plt.title('Comparison of Diff Values for CEN Proteins')
plt.xlabel('Diff (Simulated Data)')
plt.ylabel('Diff (PBRM1 Data)')
plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
plt.axvline(0, color='gray', linestyle='--', linewidth=0.8)
plt.grid(True)
plt.show()

correlation_cen, p_value_cen = pearsonr(cen_comparison_df['diff_simulated'], cen_comparison_df['diff_pbrm1'])
correlation_cen
p_value_cen 

