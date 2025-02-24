import pandas as pd
import numpy as np
from scipy.stats import t
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
from gseapy.plot import barplot, dotplot
import csv
import os

os.chdir('..')
from simko_func.simko import get_classes_by_mean_abundance, get_control_differentials, get_ko_differentials, get_significance, boxplots_for_pathways


abundance = pd.read_csv('data/simko2_data/passport_prots.csv', index_col=0)
abundance.index = abundance.index.astype(str)
abundance

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

#onyl usng significant proteins forthe gsea
diff_stats_sf = diff_stats.loc[diff_stats['adjusted_p'] < 0.05]
diff_stats_sf


#ranking proteins based in diff and creating a data frame of only this
diff_stats_for_gsea = diff_stats_sf[['protein', 'diff']]
diff_stats_for_gsea = diff_stats_for_gsea.sort_values(by='diff', ascending=False)
diff_stats_for_gsea

#saved the ranked file
diff_stats_for_gsea.to_csv('PBRM1_results_ranked.rnk', sep='\t', header=False, index=False)

#need the gene id file
results = gp.prerank(
    rnk="PBRM1_results_ranked.rnk",  # Path to your ranked file
    gene_sets="h.all.v2024.1.Hs.symbols.gmt",
    outdir="GSEA_results",  # Output directory for results
    permutation_num=100,  # Number of permutations
    seed=123,  # For reproducibility
)


#enrichment plots
gsea_hallmark_results = pd.read_csv('GSEA_results/gseapy.prerank.gene_sets.report.csv')
gsea_hallmark_results

#plots dna repair
gsea_hallmark_results["Regulation"] = gsea_hallmark_results["nes"].apply(lambda x: "Up-regulated" if x > 0 else "Down-regulated")
gsea_hallmark_results

gsea_hallmark_results = gsea_hallmark_results.sort_values("nes", ascending=False)

# Create the box plot
plt.figure(figsize=(8, 8))
sns.barplot(data=gsea_hallmark_results, x="nes", y="Term", hue="Regulation", palette={"Up-regulated": "blue", "Down-regulated": "red"})
plt.axvline(0, color="gray", linestyle="--", linewidth=0.8)  # Vertical line at NES = 0
plt.title("Pathway Regulation (Up vs Down)")
plt.xlabel("Normalized Enrichment Score (NES)")
plt.ylabel("Pathways")
plt.tight_layout()
plt.legend(title="Regulation", loc="lower right")
plt.show()
#gene set enrichment og biological processes

