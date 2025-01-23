import pandas as pd
import numpy as np
from scipy.stats import t
import matplotlib.pyplot as plt
import seaborn as sns
import random
import requests
from io import StringIO
import csv

#from importlib import reload
#import simko.simko as simko
#reload(simko)

data_dir = '~/icr/simko/simko_git/simko2_analysis-main/data/simko2_data/'


model_annotation = pd.read_csv('data/simko2_data/model_list_20240110.csv')
protein_set = ['BRCA1', 'BRCA2', 'BARD1', 'SMC6', 'SMC5', 'ARID1A', 'PBRM1', 'SMARCA1', 'SMARCC1', 'STAG1', 'STAG2', 'SMARCE1', 'SMARCD1', 'PLK1']


abundance = pd.read_csv('data/simko2_data/passport_prots.csv', index_col=0)
abundance.index = abundance.index.astype(str)



## Single SimKO
#created data frame of top 'median' and bottom 'low' cell lines based on abundance of the "ko_protein"
class_df = get_classes_by_mean_abundance(ko_proteins=['ARID1A'], abundance=abundance, n=30)
# random abundance difference of proteisn between randomly selected CLs (30) - iterated through 100 combos
control_diffs = get_control_differentials(abundance, class_df, k=30)

#box plot for proteins in 'protein set'  - plots abundance difference per protein fro each of the trials
#sns.boxplot(control_diffs.loc[control_diffs['protein'].isin(protein_set)], y='protein', x='diff')

#computes difference between mean proten abundance 
diffs = get_ko_differentials(abundance=abundance, class_df=class_df)
diffs

#is there significant difference between observed diffs (from KO sim) and differences due to chance?
diff_stats = get_significance(diffs, control_diffs, n=30)
diff_stats.sort_values(by = 'adjusted_p', ascending=True).head(20)


#getting associated proteins to the 'knockout protein 
prot_assoc = get_associated_proteins(ko_protein='ARID1A')
prot_assoc

brca1_results = diff_stats[diff_stats['protein'].isin(prot_assoc)]
brca1_results


### All SimKOs
proteins=sorted(protein_set) 
protein_count = len(proteins)
i = 1
n=30
diffs_list = []
for ko_prot in proteins:
    print("%s / %s - %s                       " % (i, protein_count, ko_prot),sep='',end="\r",flush=True)
    try:
        class_df = get_classes_by_mean_abundance(ko_proteins=[ko_prot], abundance=abundance, n=n)
        control_diffs = get_control_differentials(abundance, class_df, trials=30)
        ko_diffs = get_ko_differentials(abundance=abundance, class_df=class_df)
        ko_diffs = get_significance(ko_diffs , control_diffs, n=n)
        ko_diffs['ko'] = ko_prot
        diffs_list.append(ko_diffs)
    except ValueError:
        print(ko_prot, 'missing!')
        pass
    i += 1
simkos = pd.concat(diffs_list)
simkos = simkos[['ko', 'protein', 'diff', 'control_mean', 'control_std', 'p']]
simkos = simkos.round({'control_mean':2, 'control_std':2, 'p':5})

label = 'passport'
simkos.to_csv('%s/%s_simkos.csv' % (data_dir, label))





### Box plot
simkos = pd.read_csv('simko2_analysis-main/data/simko2_data/%s_simkos.csv' % label)
simkos
simko_set = simkos[['protein', 'ko', 'diff']]
simko_set = simko_set.loc[simko_set['protein'].isin(protein_set)]
simko_set['set'] = 'KO'

control_set = control_diffs.loc[control_diffs['protein'].isin(protein_set)]
control_set = control_set.rename(columns={'control_diff':'diff'})
control_set['set'] = 'Control'

box_set = pd.concat([simko_set, control_set])
box_set = box_set.sort_values('protein').reset_index()

palette = {
    'KO': '#F8766D',
    'Control': '#00BFC4',
}

box = sns.boxplot(box_set, y='protein', x='diff', hue='set', palette=palette, hue_order=palette.keys()).set(title='%s | %s control' % (label, control_method))
box.get_figure().savefig('figs/%s_cntl2%s_sample_boxplot.pdf' % (label, control_method))

### Heatmap

sk_hm = simkos.pivot_table(columns='protein', index='ko', values='diff')
sk_hm = sk_hm[proteins].fillna(0)
hm = sns.clustermap(sk_hm, row_cluster=False, col_cluster=False, center = 0, cmap=sns.diverging_palette(220, 20, as_cmap=True))
hm.savefig('figs/%s_cntl2%s_hm.pdf' % (label, control_method))





### COREAD ARID1A knockout predictions
abundance_filt = simko.filter_cls(abundance, model_annotation, ['Colorectal Carcinoma'], data_source=label)
abundance_filt

class_df = simko.get_classes_by_mean_abundance(ko_proteins=['ARID1A'], abundance=abundance_filt, n=15)
control_diffs = simko.get_control_differentials(abundance_filt, trials=100)
diffs = simko.get_ko_differentials(abundance=abundance_filt, class_df=class_df)
diff_stats = simko.get_significance(diffs , control_diffs, n=15)
diff_stats.loc[diff_stats['p'] < 0.01]
diff_stats.to_csv('%s/coread_arid1a_preds_%s_%s.csv' % (data_dir, label, control_method))





### Compare predictions

diff_stats_ccle = pd.read_csv('%s/coread_arid1a_preds_ccle_%s.csv' % (data_dir,  control_method))
diff_stats_passport = pd.read_csv('%s/coread_arid1a_preds_passport_%s.csv' % (data_dir,  control_method))
diff_merge = diff_stats_ccle[['protein', 'diff', 'p']].merge(diff_stats_passport[['protein', 'diff', 'p']], on='protein', suffixes=['_ccle', '_passport']).dropna()
diff_merge = diff_merge.loc[(diff_merge['p_ccle'] < 0.01) | (diff_merge['p_passport'] < 0.01)]
sns.heatmap(diff_merge[['p_ccle', 'p_passport']].sort_values('p_ccle'))