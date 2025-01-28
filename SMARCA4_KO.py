import pandas as pd
import numpy as np
from scipy.stats import t
import matplotlib.pyplot as plt
import seaborn as sns
import random
import requests
from io import StringIO
import csv

from simko_func.simko import get_classes_by_mean_abundance, get_control_differentials, get_ko_differentials, get_significance, get_pathway_assigned

abundance = pd.read_csv('data/simko2_data/passport_prots.csv', index_col=0)
abundance.index = abundance.index.astype(str)

class_df = get_classes_by_mean_abundance(ko_proteins=['SMARCA4'], abundance=abundance, n=30)
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
diff_stats['diff_copy'] = diff_stats['diff']
diff_stats[:20]

#assignign pathways
diff_stats_pathways = get_pathway_assigned(diff_stats)

#boxplot of pathway changes
boxplots_for_pathways(diff_stats_pathways, x='pathways', y='diff_copy',
                      title='SMARCA4 SimKO effects',
                      figsize=(10,8),
                      xlabel='Pathways ',
                      ylabel='Changes in pathway regulation',
                      save_path='SMARCA4 SiKO')


