import pandas as pd
import numpy as np
from scipy.stats import t
import matplotlib.pyplot as plt
import seaborn as sns
import random

#from importlib import reload
#import simko.simko as simko
#reload(simko)


model_annotation = pd.read_csv('data/simko2_data/model_list_20240110.csv')
protein_set = ['BRCA1', 'BRCA2', 'BARD1', 'SMC6', 'SMC5', 'ARID1A', 'PBRM1', 'SMARCA1', 'SMARCC1', 'STAG1', 'STAG2', 'SMARCE1', 'SMARCD1', 'PLK1']


abundance = pd.read_csv('data/simko2_data/passport_prots.csv', index_col=0)
abundance.index = abundance.index.astype(str)



## Single SimKO
#created data frame of top 'median' and bottom 'low' cell lines based on abundance of the "ko_protein"
class_df = get_classes_by_mean_abundance(ko_proteins=['BRCA1'], abundance=abundance, n=30)

# random abundance difference of proteisn between randomly selected CLs (30) - iterated through 100 combos
control_diffs = get_control_differentials(abundance, class_df, k=30)



#computes difference between mean proten abundance 
diffs = get_ko_differentials(abundance=abundance, class_df=class_df)
diffs

#is there significant difference between observed diffs (from KO sim) and differences due to chance?
diff_stats = get_significance(diffs, control_diffs, n=30)
diff_stats
diff_stats.sort_values(by = 'adjusted_p', ascending=True)
papr1 = diff_stats[diff_stats['protein'] == 'PARP1']
papr1

class_df