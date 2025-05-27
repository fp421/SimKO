import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import openpyxl
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy import stats
from simko_func import simko
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.metrics import confusion_matrix

#loading abundance data
abundance = pd.read_csv('~/icr/simko/data/simko2_data/passport_prots.csv', index_col=0)
abundance.index = abundance.index.astype(str)

#doing simko on the abundance data
#note: the p values will slightly change as the control shuffles are always different 
class_df = simko.get_classes_by_mean_abundance(ko_proteins=['PBRM1'], abundance=abundance, n=30)
control_diffs = simko.get_control_differentials(abundance, class_df, k=30)
diffs = simko.get_ko_differentials(abundance=abundance, class_df=class_df)
diffs
PBRM1_row = diffs[diffs["protein"] == "PBRM1"]
PBRM1_row
diff_stats = simko.get_significance(diffs, control_diffs, n=30)
diff_stats
#-np.log10() on pvals for volcano plot
diff_stats['log10_pval']= -np.log10(diff_stats['adjusted_p'])
diff_stats = diff_stats.sort_values(by='diff', ascending = True)
diff_stats_sf = diff_stats.loc[diff_stats['adjusted_p'] < 0.05]
diff_stats_sf

#in this current data frame, the diff coloumn acts more like a function, rathwr than just a simple value
#for ease im making another coloumn (diff_copy) --> so its just a value not a 'fucntion'
diff_stats_sf['LogFC'] = diff_stats_sf['diff']

#assigning pathways to simko results (to gey selected group of proteins)
diff_stats_pathway_prots = simko.get_pathway_assigned(diff_stats_sf)
diff_stats_pathway_prots
#only keeping pathways involve din cell cycle regulation
#keep1 = ["e2f_targets"]
#keep = ["e2f_targets", "myc_targets", "dna_repair", "g2m_checkpoint"]
#diff_stats_pathway_prots = diff_stats_pathway_prots[diff_stats_pathway_prots["pathways"].isin(keep1)].copy()
diff_stats_pathway_prots = diff_stats_pathway_prots[['protein', 'LogFC']]
diff_stats_pathway_prots['pred_label'] = (diff_stats_pathway_prots['LogFC'] > 0).astype(int) #one is upregulated
diff_stats_pathway_prots



#now doing the same labelling with experimental data
pbrm1_exp = pd.read_excel('~/icr/simko/data/pbrm1_experiment_data/all_pbrm1_results.xlsx')
pbrm1_exp.columns = pbrm1_exp.columns.str.replace(' ', '_').str.replace('-', '').str.replace('/', '_')
pbrm1_exp[['protein', 'LogFC']] = pbrm1_exp[['Gene_Names_(primary)', 'mean_log2_KO_WT']]
pbrm1_exp = pbrm1_exp[['protein', 'LogFC']]
pbrm1_exp['true_label'] = (pbrm1_exp['LogFC'] > 0).astype(int)

merged_df = pd.merge(diff_stats_pathway_prots[['protein', 'LogFC', 'pred_label']],
                     pbrm1_exp[['protein', 'LogFC', 'true_label']],
                     on='protein',
                     suffixes=('_pred', '_true'))

merged_df

fpr, tpr, thresholds = roc_curve(merged_df['true_label'], merged_df['LogFC_pred'])
auc_score = roc_auc_score(merged_df['true_label'], merged_df['LogFC_pred'])

print(f'AUC = {auc_score:.3f}')

# --- 5. Plot the ROC curve ---
plt.figure(figsize=(8,8))
plt.plot(fpr, tpr, label=f'ROC (AUC = {auc_score:.3f})')
plt.plot([0,1], [0,1], linestyle='--', color='gray', label='Random chance')
plt.xlabel('False Positive Rate', fontsize=16)
plt.ylabel('True Positive Rate', fontsize=16)
plt.title('ROC Curve: SimKO vs. Experimental', fontsize=18)
plt.legend(loc='lower right', fontsize=14)
plt.grid(True)
plt.savefig("roc_for_gsea_poster.pdf", format='pdf')
plt.show()

# assuming df has columns ['logFC_pred','logFC_true']
pearson_r, pearson_p = pearsonr(merged_df['LogFC_pred'], merged_df['LogFC_true'])
spearman_rho, spearman_p = spearmanr(merged_df['LogFC_pred'], merged_df['LogFC_true'])

print(f"Pearson r = {pearson_r:.3f} (p = {pearson_p:.1e})")
print(f"Spearman ρ = {spearman_rho:.3f} (p = {spearman_p:.1e})")

# 1) Binarize predictions & truth (adjust threshold if you like)
y_pred = (merged_df['LogFC_pred'] > 0).astype(int)   # 1 means “predicted up”
y_true = (merged_df['LogFC_true'] > 0).astype(int)   # 1 means “actually up”

# 2) Compute confusion matrix
#    Returns array([[TN, FP],
#                    [FN, TP]])
tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()

# 3) Compute sensitivity & specificity
sensitivity = tp / (tp + fn)   # TP / (TP + FN)
specificity = tn / (tn + fp)   # TN / (TN + FP)

print(f"Sensitivity (recall for up-regulation) = {sensitivity:.3f}")
print(f"Specificity (recall for down-regulation) = {specificity:.3f}")


#


