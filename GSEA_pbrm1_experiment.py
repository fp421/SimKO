import pandas as pd
import numpy as np
from scipy.stats import t
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
from gseapy.plot import barplot, dotplot
import csv
import os
from simko_func import simko

#GSEA on CVAE simulated PBRM1 knokcout results
#gsea on CVAE results 
cvae_pbrm1_results = pd.read_csv('protein_shift_summary.csv')
cvae_pbrm1_results = cvae_pbrm1_results.sort_values(by='diff', ascending=False)
cvae_pbrm1_results = cvae_pbrm1_results[['protein', 'diff']]

#saved the ranked file
cvae_pbrm1_results.to_csv('PBRM1_CVAE_ranked.rnk', sep='\t', header=False, index=False)

#need the gene id file
results = gp.prerank(
    rnk="PBRM1_CVAE_ranked.rnk",  # Path to your ranked file
    gene_sets="h.all.v2024.1.Hs.symbols.gmt",
    outdir="CVAE_PBRM1_GSEA_results",  # Output directory for results
    permutation_num=100,  # Number of permutations
    seed=123,  # For reproducibility
)

pbrm1_experiment = pd.read_excel('~/icr/simko/data/pbrm1_experiment_data/all_pbrm1_results.xlsx')
pbrm1_experiment.columns = pbrm1_experiment.columns.str.replace(' ', '_').str.replace('-', '').str.replace('/', '_')


#ranking proteins based in diff and creating a data frame of only this
pbrm1_exp_for_gsea = pbrm1_experiment[['Gene_Names_(primary)', 'mean_log2_KO_WT']]
pbrm1_exp_for_gsea = pbrm1_exp_for_gsea.sort_values(by='mean_log2_KO_WT', ascending=False)
pbrm1_exp_for_gsea

#saved the ranked file
pbrm1_exp_for_gsea.to_csv('PBRM1_avg_exp_results_ranked.rnk', sep='\t', header=False, index=False)

#need the gene id file
results = gp.prerank(
    rnk="PBRM1_avg_exp_results_ranked.rnk",  # Path to your ranked file
    gene_sets="h.all.v2024.1.Hs.symbols.gmt",
    outdir="PBRM1_exp_avg_GSEA_results",  # Output directory for results
    permutation_num=100,  # Number of permutations
    seed=123,  # For reproducibility
)

#same again but just for kidney cell line
pbrm1_experiment
#ranking proteins based in diff and creating a data frame of only this
pbrm1_exp_for_gsea_kid = pbrm1_experiment[['Gene_Names_(primary)', 'median_786O_log2_PBRM1_KO_WT_norm']]
pbrm1_exp_for_gsea_kid = pbrm1_exp_for_gsea_kid.sort_values(by='median_786O_log2_PBRM1_KO_WT_norm', ascending=False)
pbrm1_exp_for_gsea_kid

#saved the ranked file
pbrm1_exp_for_gsea_kid.to_csv('PBRM1_kidney_exp_results_ranked.rnk', sep='\t', header=False, index=False)

#need the gene id file
results = gp.prerank(
    rnk="PBRM1_kidney_exp_results_ranked.rnk",  # Path to your ranked file
    gene_sets="h.all.v2024.1.Hs.symbols.gmt",
    outdir="PBRM1_exp_kidney_GSEA_results",  # Output directory for results
    permutation_num=100,  # Number of permutations
    seed=123,  # For reproducibility
)







#now doing the same for Arid1a knockout
abundance_arid1a = pd.read_csv('data/simko2_data/passport_prots.csv', index_col=0)
abundance_arid1a.index = abundance_arid1a.index.astype(str)
abundance_arid1a

#doing only lung tissues
model_lists = pd.read_csv('~/icr/simko/data/simko2_data/model_list_20240110.csv', index_col=0)

#need to make sure that the cell lines in model_list matched us with cell lines in abundance data
abundance_cell_lines = set(abundance_arid1a.columns)
model_list_filt = model_lists[model_lists['model_name'].isin(abundance_cell_lines)]
unique_tissues = model_list_filt['tissue'].unique()
unique_tissues
#column names to sue: model_name and tissues
# making a dictionary
grouped = model_list_filt.groupby('tissue')['model_name'].apply(list)
tissue_to_cell_lines = grouped.to_dict()
print(tissue_to_cell_lines)
#finding out how many cell lines for each tissue
for key, value in tissue_to_cell_lines.items():
    if isinstance(value, list):  # Check if the value is a list
        print(f"The length of the list for {key} is: {len(value)}")


specific_CL = tissue_to_cell_lines["Lung"]
abundance_specific_CL = abundance_arid1a[specific_CL]
abundance_specific_CL

class_df = simko.get_classes_by_mean_abundance(ko_proteins=['ARID1A'], abundance=abundance_specific_CL, n=30)
control_diffs = simko.get_control_differentials(abundance_specific_CL, class_df, k=30)
diffs = simko.get_ko_differentials(abundance=abundance_specific_CL, class_df=class_df)
PBRM1_row = diffs[diffs["protein"] == "PBRM1"]
PBRM1_row
diff_stats = simko.get_significance(diffs, control_diffs, n=30)
diff_stats
#-np.log10() on pvals for volcano plot
diff_stats['log10_pval']= -np.log10(diff_stats['adjusted_p'])
diff_stats = diff_stats.sort_values(by='log10_pval', ascending = False)
diff_stats[:20]

#onyl usng significant proteins forthe gsea
diff_stats_sf = diff_stats.loc[diff_stats['adjusted_p'] < 0.05]
diff_stats_sf

#ranking proteins based in diff and creating a data frame of only this
diff_stats_for_gsea_arid1a = diff_stats_sf[['protein', 'diff']]
diff_stats_for_gsea_arid1a = diff_stats_for_gsea_arid1a.sort_values(by='diff', ascending=False)
diff_stats_for_gsea_arid1a.value_counts()

#saved the ranked file
diff_stats_for_gsea_arid1a.to_csv('ARID1A_LUNG_results_ranked.rnk', sep='\t', header=False, index=False)

#need the gene id file
results = gp.prerank(
    rnk="ARID1A_LUNG_results_ranked.rnk",  # Path to your ranked file
    gene_sets="h.all.v2024.1.Hs.symbols.gmt",
    outdir="ARID1A_LUNG_GSEA_results",  # Output directory for results
    permutation_num=100,  # Number of permutations
    seed=123, 
    processes=1 # For reproducibility
)


#now doing the same with TP53
class_df1 = simko.get_classes_by_mean_abundance(ko_proteins=['TP53'], abundance=abundance_arid1a, n=30)
control_diffs1 = simko.get_control_differentials(abundance_arid1a, class_df1, k=30)
diffs1 = simko.get_ko_differentials(abundance=abundance_arid1a, class_df=class_df1)
diff_stats1 = simko.get_significance(diffs1, control_diffs1, n=30)
diff_stats1
#-np.log10() on pvals for volcano plot
diff_stats1['log10_pval']= -np.log10(diff_stats1['adjusted_p'])
diff_stats1 = diff_stats1.sort_values(by='log10_pval', ascending = False)
diff_stats1[:20]

#onyl usng significant proteins forthe gsea
diff_stats_sf1 = diff_stats1.loc[diff_stats1['adjusted_p'] < 0.05]
diff_stats_sf1

#ranking proteins based in diff and creating a data frame of only this
diff_stats_for_gsea_tp53 = diff_stats_sf1[['protein', 'diff']]
diff_stats_for_gsea_tp53 = diff_stats_for_gsea_tp53.sort_values(by='diff', ascending=False)
diff_stats_for_gsea_tp53

#saved the ranked file
diff_stats_for_gsea_tp53.to_csv('TP53_results_ranked.rnk', sep='\t', header=False, index=False)

#need the gene id file
results = gp.prerank(
    rnk="TP53_results_ranked.rnk",  # Path to your ranked file
    gene_sets="h.all.v2024.1.Hs.symbols.gmt",
    outdir="TP53_GSEA_results",  # Output directory for results
    permutation_num=100,  # Number of permutations
    seed=123,  # For reproducibility
)