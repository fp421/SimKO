import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import openpyxl
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from simko_func import simko


abundance1 = pd.read_csv('~/icr/simko/data/simko2_data/passport_prots.csv', index_col=0)
abundance1.index = abundance1.index.astype(str)


#finding correlation of all proteins to pbrm1
pbrm1 = abundance1.loc["PBRM1"]
pbrm1

correlations = abundance1.corrwith(pbrm1, axis=1)
correlations

#now putting these correlations onto the scatter plot i previously made
#doing simko
class_df = simko.get_classes_by_mean_abundance(ko_proteins=['PBRM1'], abundance=abundance1, n=30)
control_diffs = simko.get_control_differentials(abundance1, class_df, k=30)
diffs = simko.get_ko_differentials(abundance=abundance1, class_df=class_df)
diffs
PBRM1_row = diffs[diffs["protein"] == "PBRM1"]
PBRM1_row
diff_stats = simko.get_significance(diffs, control_diffs, n=30)
diff_stats
#-np.log10() on pvals for volcano plot
diff_stats['log10_pval']= -np.log10(diff_stats['adjusted_p'])
diff_stats = diff_stats.sort_values(by='log10_pval', ascending = False)
diff_stats[:20]

#downloading pbrm_1 experiment
PBRM1_data1 = pd.read_excel('~/icr/simko/data/pbrm1_experiment_data/all_pbrm1_results.xlsx')
PBRM1_data1.columns = PBRM1_data1.columns.str.replace(' ', '_').str.replace('-', '').str.replace('/', '_')
PBRM1_data1

#creating scatter plot (for all proteins not jsut significant)
proteins = diff_stats['protein']
common_proteins = set(proteins).intersection(set(PBRM1_data1['Gene_Names_(primary)']))
simulated_diff1 = diff_stats.loc[diff_stats['protein'].isin(common_proteins), 
                                ['protein', 'diff']].set_index('protein')
pbrm1_diff = PBRM1_data1.loc[PBRM1_data1['Gene_Names_(primary)'].isin(common_proteins),
                             ['Gene_Names_(primary)', 'mean_log2_KO_WT']]
pbrm1_diff = pbrm1_diff.rename(columns={'Gene_Names_(primary)': 'protein', 'mean_log2_KO_WT': 'diff'}).set_index('protein')

comparison_df = simulated_diff1.join(pbrm1_diff, lsuffix='_simulated', rsuffix='_pbrm1')
comparison_df['correlation'] = comparison_df.index.map(correlations)

#now doing the scatter plot
plt.figure(figsize=(18, 13))
sns.scatterplot(
    x = comparison_df['diff_simulated'],
    y = comparison_df['diff_pbrm1'],
    color = 'steelblue',
    edgecolor='None',
    alpha = 0.8,
    s = 100
)
#keeping line of best fit
sns.regplot(
    x=comparison_df['diff_simulated'],
    y=comparison_df['diff_pbrm1'],
    scatter=False,  # Prevents regplot from drawing its own scatter points
    color='red',
    ci=None,
    line_kws={'linewidth': 2}  # Customization of the line
)
plt.title('Comparison of LogFC Between Simulated and Experimental Data')
plt.xlabel('LogFC (Simulated Data)')
plt.ylabel('LogFC (PBRM1 Data)')
plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
plt.axvline(0, color='gray', linestyle='--', linewidth=0.8)
plt.grid(True)
plt.show()

#finding correlations
correlation, p_value = pearsonr(comparison_df['diff_simulated'], comparison_df['diff_pbrm1'])
correlation
p_value



#scatter plot only using significant proteins from simko
significant_protein = diff_stats.loc[diff_stats["adjusted_p"] < 0.05, 'protein']
common_proteins1 = set(significant_protein).intersection(set(PBRM1_data1['Gene_Names_(primary)']))
simulated_diff1 = diff_stats.loc[diff_stats['protein'].isin(common_proteins1), 
                                ['protein', 'diff']].set_index('protein')
pbrm1_diff1 = PBRM1_data1.loc[PBRM1_data1['Gene_Names_(primary)'].isin(common_proteins1),
                             ['Gene_Names_(primary)', 'mean_log2_KO_WT']]
pbrm1_diff1 = pbrm1_diff1.rename(columns={'Gene_Names_(primary)': 'protein', 'mean_log2_KO_WT': 'diff'}).set_index('protein')

# Combine into a single DataFrame for comparison
comparison_df1 = simulated_diff1.join(pbrm1_diff1, lsuffix='_simulated', rsuffix='_pbrm1')
comparison_df1['correlation'] = comparison_df1.index.map(correlations)
comparison_df1

#first scatter plot wont conatain correlation information or labels 
plt.figure(figsize=(18, 13))
sns.scatterplot(
    x = comparison_df1['diff_simulated'],
    y = comparison_df1['diff_pbrm1'],
    color = 'steelblue',
    edgecolor='None',
    alpha = 0.8,
    s = 100
)
#keeping line of best fit
sns.regplot(
    x=comparison_df1['diff_simulated'],
    y=comparison_df1['diff_pbrm1'],
    scatter=False,  # Prevents regplot from drawing its own scatter points
    color='red',
    ci=None,
    line_kws={'linewidth': 2}  # Customization of the line
)
plt.title('Comparison of LogFC Between Simulated and Experimental Data (significant proteins)')
plt.xlabel('LogFC (Simulated Data)')
plt.ylabel('LogFC (PBRM1 Data)')
plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
plt.axvline(0, color='gray', linestyle='--', linewidth=0.8)
plt.grid(True)
plt.show()

#finding correlations
correlation1, p_value1 = pearsonr(comparison_df1['diff_simulated'], comparison_df1['diff_pbrm1'])
correlation1
p_value1


#plotting the scatter plot - but now with correlation information
plt.figure(figsize=(18, 13))
sns.scatterplot(
    x = comparison_df1['diff_simulated'],
    y = comparison_df1['diff_pbrm1'],
    hue=comparison_df1['correlation'],
    palette='Spectral',
    edgecolor='k'
)
#keeping line of best fit
sns.regplot(
    x=comparison_df1['diff_simulated'],
    y=comparison_df1['diff_pbrm1'],
    scatter=False,  # Prevents regplot from drawing its own scatter points
    color='steelblue',
    ci=None,
    line_kws={'linewidth': 2}  # Customization of the line
)
plt.legend(
    title="Individual Protein Correlation to PBRM1", title_fontsize=12, fontsize=10, loc='best')    
plt.title('Comparison of logFC Values Between Simulated and PBRM1 Data with correlation information')
plt.xlabel('LogFC (Simulated Data)')
plt.ylabel('LogFC (PBRM1 Data)')
plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
plt.axvline(0, color='gray', linestyle='--', linewidth=0.8)
plt.grid(True)
plt.show()


#labelling specific proteins on the correlation plots
label_proteins = comparison_df1[
    ((comparison_df1['correlation']<-0.1)&
    (comparison_df1['diff_simulated']<-0.87)&
    (comparison_df1['diff_pbrm1']<-0))|
    ((comparison_df1['correlation']==1))|
    ((comparison_df1['correlation']>0.08)&
     (comparison_df1['diff_simulated']>0.3)&
     (comparison_df1['diff_pbrm1']>0.68))
]
plt.figure(figsize=(18, 13))
sns.scatterplot(
    x = comparison_df1['diff_simulated'],
    y = comparison_df1['diff_pbrm1'],
    hue=comparison_df1['correlation'],
    palette='Spectral',
    edgecolor='k'
)
#keeping line of best fit
sns.regplot(
    x=comparison_df1['diff_simulated'],
    y=comparison_df1['diff_pbrm1'],
    scatter=False,  # Prevents regplot from drawing its own scatter points
    color='steelblue',
    ci=None,
    line_kws={'linewidth': 2}  # Customization of the line
)
#annotating specific proteins
for index, row in label_proteins.iterrows():
    plt.text(
        row['diff_simulated'], 
        row['diff_pbrm1'], 
        row.name,  # Assuming the index contains the protein name
        fontsize=14, 
        ha='right', 
        color='black'
    )
plt.legend(
    title="Individual Protein Correlation to PBRM1", title_fontsize=12, fontsize=10, loc='best')    
plt.title('Comparison of logFC Values Between Simulated and PBRM1 Data with correlation information')
plt.xlabel('LogFC (Simulated Data)')
plt.ylabel('LogFC (PBRM1 Data)')
plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
plt.axvline(0, color='gray', linestyle='--', linewidth=0.8)
plt.grid(True)
plt.show()


#creating a table of just the labelle proteins
label_proteins_df = label_proteins.sort_values(by='correlation', ascending=True)
label_proteins_df.to_csv('labelled_proteins_table.csv')

#labelle dproteins but only their significance (from diff_stats)
label_protein_names = ['LIMD2', 'NOS3', 'GTSF1', 'SCNM1', 'GFRA1', 'SOD3', 'MARVELD2', 'AUTS2']
label_proteins_stats = diff_stats[diff_stats['protein'].isin(label_protein_names)]
selected_columns = ['protein', 'diff', 'p']
label_proteins_stats = label_proteins_stats.loc[:, selected_columns]
label_proteins_stats








#same scatter plot and gradient but onyl for the cen proteins - all cen proteins - not just sig ones

#cen proteins
CEN_proteins = ['MIS18a','MIS18BP1','RSF1','KAT7','SMC2','SMC4','NCAPD3','NCAPG2','NCAPH2',
                'CENPC','CENPT','CENPS','CENPX','CENPO','CENPU','CENPL','CENPN','CENPH','CENPI','CENPK',
                'KNL1','ZWINT','NDC80','NUF2','SPC24','SPC25','DSN1','MIS12','NSL1','PMF1',
                'INCENP','BIRC5','CDCA8','AURKB',
                'CENPE','CENPF','CENPB','CENPV','CENPJ','CDCA2',
                'CBX1','SIRT6','UHRF2','SUV39H1','SETDB1','DAXX','SSRP1','SMARCAD1',
                'ATRX','BAZ1A','BAZ1B','CBX3','CBX5','CHRAC1','DNMT1','EZH2','HELLS','KDM4A','KMT5C','LRWD1','MACROH2A1','POLE3','SMARCA5']

#want all cen proteins - not just the significant ones
#just quickly re did the code so that it isn't just significant
common_proteins_notsig = set(diff_stats['protein']).intersection(set(PBRM1_data1['Gene_Names_(primary)']))
simulated_diff1 = diff_stats.loc[diff_stats['protein'].isin(common_proteins_notsig), 
                                ['protein', 'diff']].set_index('protein')
pbrm1_diff1 = PBRM1_data1.loc[PBRM1_data1['Gene_Names_(primary)'].isin(common_proteins_notsig),
                             ['Gene_Names_(primary)', 'mean_log2_KO_WT']]
pbrm1_diff1 = pbrm1_diff1.rename(columns={'Gene_Names_(primary)': 'protein', 'mean_log2_KO_WT': 'diff'}).set_index('protein')

# Combine into a single DataFrame for comparison
comparison_df1 = simulated_diff1.join(pbrm1_diff1, lsuffix='_simulated', rsuffix='_pbrm1')
comparison_df1['correlation'] = comparison_df1.index.map(correlations)


cen_comparison_df1 = comparison_df1.loc[comparison_df1.index.isin(CEN_proteins)]
cen_comparison_df1['correlation'] = cen_comparison_df1.index.map(correlations)
cen_comparison_df1



#now scatter plot again - points coloured based on correaltion to pbrm1
plt.figure(figsize=(8, 6))
sns.scatterplot(
    x = cen_comparison_df1['diff_simulated'],
    y = cen_comparison_df1['diff_pbrm1'],
    hue=cen_comparison_df1['correlation'],
    palette = 'Spectral',
    edgecolor = 'k'
)
sns.regplot(
    x = cen_comparison_df1['diff_simulated'],
    y = cen_comparison_df1['diff_pbrm1'],
    scatter = False,
    ci = None,
    line_kws={'color': 'black', 'linewidth': 1}
)
# Add plot labels
plt.title('Comparison of LogFC Values for CEN proteins Between Simulated and PBRM1 Data')
plt.xlabel('Diff (Simulated Data)')
plt.ylabel('Diff (PBRM1 Data)')
plt.axhline(0, color='gray', linestyle='--', linewidth=0.6)
plt.axvline(0, color='gray', linestyle='--', linewidth=0.6)
plt.grid(True)
plt.show()


#scatter plot for cen proteins with line of best fit - not hue gradient - just redoing it
plt.figure(figsize=(8, 6))
sns.regplot(
    x = cen_comparison_df1['diff_simulated'],
    y = cen_comparison_df1['diff_pbrm1'],
    ci = None,
    line_kws={'color': 'red', 'linewidth': 1.5}
)
# Add plot labels
plt.title('Comparison of LogFC Values for CEN proteins Between Simulated and PBRM1 Data')
plt.xlabel('LogFC (Simulated Data)')
plt.ylabel('LogFC (PBRM1 Data)')
plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
plt.axvline(0, color='gray', linestyle='--', linewidth=0.8)
plt.grid(True)

plt.show()

correlation2, p_value2 = pearsonr(cen_comparison_df1['diff_simulated'], cen_comparison_df1['diff_pbrm1'])
correlation2
p_value2


#cen proteins but only the significant ones
significant_protein = diff_stats.loc[diff_stats["adjusted_p"] < 0.05, 'protein']
common_proteins1 = set(significant_protein).intersection(set(PBRM1_data1['Gene_Names_(primary)']))
simulated_diff1 = diff_stats.loc[diff_stats['protein'].isin(common_proteins1), 
                                ['protein', 'diff']].set_index('protein')
pbrm1_diff1 = PBRM1_data1.loc[PBRM1_data1['Gene_Names_(primary)'].isin(common_proteins1),
                             ['Gene_Names_(primary)', 'mean_log2_KO_WT']]
pbrm1_diff1 = pbrm1_diff1.rename(columns={'Gene_Names_(primary)': 'protein', 'mean_log2_KO_WT': 'diff'}).set_index('protein')

# Combine into a single DataFrame for comparison
comparison_df1 = simulated_diff1.join(pbrm1_diff1, lsuffix='_simulated', rsuffix='_pbrm1')
comparison_df1['correlation'] = comparison_df1.index.map(correlations)


cen_comparison_sig = comparison_df1.loc[comparison_df1.index.isin(CEN_proteins)]
cen_comparison_sig['correlation'] = cen_comparison_sig.index.map(correlations)



#scatter plot to compare significance with correlation?
merged_df = cen_comparison_sig.merge(diff_stats[['protein', 'adjusted_p']], on='protein', how='left')
#plotting correlation vs significance
plt.figure(figsize=(8, 6))
sns.scatterplot(
    x = merged_df['correlation'],
    y = merged_df['adjusted_p']
)

# Add plot labels
plt.title('Comparison of Correlation with Significance')
plt.xlabel('Correlation')
plt.ylabel('P Value')
plt.axhline(0, color='gray', linestyle='--', linewidth=0.6)
plt.axvline(0, color='gray', linestyle='--', linewidth=0.6)
plt.grid(True)
plt.show()

#correaltion between signifciance and correlation 
correlation3, p_value3 = pearsonr(merged_df['correlation'], merged_df['adjusted_p'])
correlation3
p_value3



#scatter plot for SIGNIFICANT cen proteins with line of best fit - not hue gradient - just redoing it
plt.figure(figsize=(8, 6))
sns.regplot(
    x = cen_comparison_sig['diff_simulated'],
    y = cen_comparison_sig['diff_pbrm1'],
    ci = None,
    line_kws={'color': 'red', 'linewidth': 1.5}
)
# Add plot labels
plt.title('Comparison of LogFC Values for CEN proteins Between Simulated and PBRM1 Data')
plt.xlabel('LogFC (Simulated Data)')
plt.ylabel('LogFC (PBRM1 Data)')
plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
plt.axvline(0, color='gray', linestyle='--', linewidth=0.8)
plt.grid(True)

plt.show()

correlation3, p_value3 = pearsonr(cen_comparison_sig['diff_simulated'], cen_comparison_sig['diff_pbrm1'])
correlation3
p_value2


#significant centromere proteins but with correlation labelled
#now scatter plot again - points coloured based on correaltion to pbrm1
plt.figure(figsize=(8, 6))
sns.scatterplot(
    x = cen_comparison_sig['diff_simulated'],
    y = cen_comparison_sig['diff_pbrm1'],
    hue=cen_comparison_sig['correlation'],
    palette = 'Spectral',
    edgecolor = 'k'
)
sns.regplot(
    x = cen_comparison_sig['diff_simulated'],
    y = cen_comparison_sig['diff_pbrm1'],
    scatter = False,
    ci = None,
    line_kws={'color': 'black', 'linewidth': 1}
)
# Add plot labels
plt.title('Comparison of LogFC Values for CEN proteins Between Simulated and PBRM1 Data')
plt.xlabel('Diff (Simulated Data)')
plt.ylabel('Diff (PBRM1 Data)')
plt.axhline(0, color='gray', linestyle='--', linewidth=0.6)
plt.axvline(0, color='gray', linestyle='--', linewidth=0.6)
plt.grid(True)
plt.show()

cen_comparison_sig


#looking at proteins associated with 
pbrm1_assoc_proteins = simko.get_associated_proteins(ko_protein='PBRM1')

correlations[pbrm1_assoc_proteins]

diff_stats_assoc = diff_stats[diff_stats['protein'].isin(pbrm1_assoc_proteins)]
diff_stats_assoc

##correlation graph for the proteins known to be associated with pbrm1
#will just only select the associate dportiens from comparison_df1
assoc_p_comparison_df1 = comparison_df1.loc[comparison_df1.index.isin(pbrm1_assoc_proteins)]
assoc_p_comparison_df1['correlation'] = assoc_p_comparison_df1.index.map(correlations)
assoc_p_comparison_df1

#plot with line of best fit
plt.figure(figsize=(8, 6))
sns.scatterplot(
    x = assoc_p_comparison_df1['diff_simulated'],
    y = assoc_p_comparison_df1['diff_pbrm1'],
    hue=assoc_p_comparison_df1['correlation'],
    palette = 'Spectral',
    edgecolor = 'k',
    size=70
)
sns.regplot(
    x=assoc_p_comparison_df1['diff_simulated'],
    y=assoc_p_comparison_df1['diff_pbrm1'],
    scatter=False,  # Prevents regplot from drawing its own scatter points
    color='steelblue',
    ci=None,
    line_kws={'linewidth': 2})

for i, row in assoc_p_comparison_df1.iterrows():
    plt.annotate(
        i,  # Use the index value (protein name) as the label
        (row['diff_simulated'], row['diff_pbrm1']),
        textcoords="offset points",
        xytext=(5, 5),  # Offset the text slightly from the point
        ha='right',
        fontsize=8
    )
# Add plot labels
plt.title('Comparison of Diff Values for PBRM1 assocuated proteins Between Simulated and PBRM1 Data')
plt.xlabel('Diff (Simulated Data)')
plt.ylabel('Diff (PBRM1 Data)')
plt.axhline(0, color='black', linestyle='--', linewidth=0.8)
plt.axvline(0, color='black', linestyle='--', linewidth=0.8)
plt.grid(True)
plt.show()

#correlation between experiment and simko for pbrm1 associated proteins
correlation, p_value = pearsonr(assoc_p_comparison_df1['diff_simulated'], assoc_p_comparison_df1['diff_pbrm1'])
correlation #correlation = 0.573
p_value #0.083 - not significant but not a lot of data



abundance1