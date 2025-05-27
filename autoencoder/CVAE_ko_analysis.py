import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import openpyxl
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy import stats
from simko_func import simko

cvae_pbrm1_results = pd.read_csv('protein_shift_summary.csv')
cvae_pbrm1_results['diff_copy'] = cvae_pbrm1_results['diff']
cvae_pbrm1_results

pathways = simko.get_pathway_assigned(cvae_pbrm1_results)

pathway_means = simko.get_pathway_stats(pathways, operation='median')
pathway_means

simko.boxplots_for_pathways(pathways, x='pathways', y='diff_copy',
                      title='Skin (n=55)',
                      figsize=(10,8),
                      xlabel=' ',
                      ylabel=' ')


#centromere protein correlations 
CEN_proteins = ['MIS18a','MIS18BP1','RSF1','KAT7','SMC2','SMC4','NCAPD3','NCAPG2','NCAPH2',
                'CENPC','CENPT','CENPS','CENPX','CENPO','CENPU','CENPL','CENPN','CENPH','CENPI','CENPK',
                'KNL1','ZWINT','NDC80','NUF2','SPC24','SPC25','DSN1','MIS12','NSL1','PMF1',
                'INCENP','BIRC5','CDCA8','AURKB',
                'CENPE','CENPF','CENPB','CENPV','CENPJ','CDCA2',
                'CBX1','SIRT6','UHRF2','SUV39H1','SETDB1','DAXX','SSRP1','SMARCAD1',
                'ATRX','BAZ1A','BAZ1B','CBX3','CBX5','CHRAC1','DNMT1','EZH2','HELLS','KDM4A','KMT5C','LRWD1','MACROH2A1','POLE3','SMARCA5']

cvae_cen_prots = cvae_pbrm1_results.loc[cvae_pbrm1_results['protein'].isin(CEN_proteins)]
cvae_cen_prots

#correlation to experimental results
pbrm1_exp = pd.read_excel('~/icr/simko/data/pbrm1_experiment_data/all_pbrm1_results.xlsx')
pbrm1_exp.columns = pbrm1_exp.columns.str.replace(' ', '_').str.replace('-', '').str.replace('/', '_')
pbrm1_exp

common_proteins = set(cvae_cen_prots['protein']).intersection(set(pbrm1_exp['Gene_Names_(primary)']))

cvae_df = cvae_cen_prots.loc[cvae_cen_prots['protein'].isin(common_proteins), 
                                ['protein', 'diff']].set_index('protein')
pbrm1_df = pbrm1_exp.loc[pbrm1_exp['Gene_Names_(primary)'].isin(common_proteins),
                             ['Gene_Names_(primary)', 'mean_log2_KO_WT']]
pbrm1_df = pbrm1_df.rename(columns={'Gene_Names_(primary)': 'protein', 'mean_log2_KO_WT': 'diff'}).set_index('protein')

# Combine into a single DataFrame for comparison
comparison_df = cvae_df.join(pbrm1_df, lsuffix='_simulated', rsuffix='_pbrm1')
comparison_df.sort_values(by='diff_pbrm1', ascending=False)
#correlation values
correlation, p_value = pearsonr(comparison_df['diff_simulated'], comparison_df['diff_pbrm1'])
correlation

#plot
plt.figure(figsize=(8, 8))
sns.scatterplot(
    x = comparison_df['diff_simulated'],
    y = comparison_df['diff_pbrm1'],
    edgecolor = 'k'
)
sns.regplot(
    x = comparison_df['diff_simulated'],
    y = comparison_df['diff_pbrm1'],
    scatter = False,
    ci = None,
    line_kws={'color': 'black', 'linewidth': 1}
)
plt.title('Comparison of LogFC Values for CEN proteins Between CVAE Simulated and PBRM1 Data')
plt.xlabel('Diff (CVAE Simulated Data)')
plt.ylabel('Diff (PBRM1 Data)')
plt.axhline(0, color='gray', linestyle='--', linewidth=0.6)
plt.axvline(0, color='gray', linestyle='--', linewidth=0.6)
plt.grid(True)
plt.show()

#checking cvae model hasnt just done it by correlations
prots_data = pd.read_csv('~/icr/simko/data/simko2_data/passport_prots.csv', index_col=0)
pbrm1 = prots_data.loc['PBRM1']
correlations = prots_data.corrwith(pbrm1, axis=1)
correlations.sort_values(ascending=False).head(50)
correlations.name = 'corr_to_pbrm1'



#comparing correlation and abundance change from cvae results side by side
cvae_diffs = cvae_pbrm1_results[['protein', 'diff']]
cvae_diffs['corr_to_pbrm1'] = cvae_diffs['protein'].map(correlations)
cvae_diffs.sort_values(by='corr_to_pbrm1', ascending=False)

plt.figure(figsize=(8, 8))
sns.scatterplot(
    x = cvae_diffs['diff'],
    y = cvae_diffs['corr_to_pbrm1'],
    edgecolor = 'k'
)
sns.regplot(
    x = cvae_diffs['diff'],
    y = cvae_diffs['corr_to_pbrm1'],
    scatter = False,
    ci = None,
    line_kws={'color': 'black', 'linewidth': 1}
)
plt.title('Comparison of LogFC changes (CVAE simulated) to protein correlations to PBRM1')
plt.xlabel('Diff (CVAE Simulated Data)')
plt.ylabel('Correlation to PBRM1)')
plt.axhline(0, color='gray', linestyle='--', linewidth=0.6)
plt.axvline(0, color='gray', linestyle='--', linewidth=0.6)
plt.grid(True)
plt.show()


what


