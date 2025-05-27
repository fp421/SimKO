import pandas as pd
import numpy as np
from scipy.stats import t
import matplotlib.pyplot as plt
import seaborn as sns
import random
import requests
from io import StringIO
import csv
from simko_func import simko
import matplotlib.patches as mpatches
from matplotlib.lines   import Line2D
from matplotlib.patches import Patch
from scipy.stats import ttest_rel, wilcoxon
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

abundance = pd.read_csv('data/simko2_data/passport_prots.csv', index_col=0)
abundance.index = abundance.index.astype(str)

## Single SimKO
#created data frame of top 'median' and bottom 'low' cell lines based on abundance of the "ko_protein"
# 2) get the top-30 / bottom-30 cell-lines for each KO

proteins = ('SMARCC1', 'SMARCC2', 'SMARCE1', 'SMARCB1',
        'SMARCD1', 'SMARCD2', 'SMARCD3', 'SMARCA2', 'SMARCA4', 'BCL7A',
        'BCL7B', 'BCL7C', 'ACTL6A')

class_pbrm1  = simko.get_classes_by_mean_abundance(ko_proteins=['PBRM1'], abundance=abundance, n=30)
class_arid1a = simko.get_classes_by_mean_abundance(ko_proteins=['ARID1A'], abundance=abundance, n=30)




# if you’re not sure what the labels are, inspect:
print(class_pbrm1['class'].unique())  
# → e.g. ['high','low']

# pull out the cell-line names in each group
high_p = class_pbrm1.index[class_pbrm1['class']=='median']
low_p  = class_pbrm1.index[class_pbrm1['class']=='low']

high_a = class_arid1a.index[class_arid1a['class']=='median']
low_a  = class_arid1a.index[class_arid1a['class']=='low']


# 3) compute Δ = mean(high) – mean(low)  for each protein
abundance = abundance.T
mean_p_high = abundance.loc[high_p, proteins].mean()
mean_p_low  = abundance.loc[low_p,  proteins].mean()
delta_pbrm1 = mean_p_high - mean_p_low

mean_a_high = abundance.loc[high_a, proteins].mean()
mean_a_low  = abundance.loc[low_a,  proteins].mean()
delta_arid1a = mean_a_high - mean_a_low


# 4) pack into a DataFrame
df_delta = pd.DataFrame({
    'PBRM1_Δ':  delta_pbrm1,
    'ARID1A_Δ': delta_arid1a
})

print(df_delta)


# 5) test whether those two Δ-vectors differ
#    (a) paired t-test across proteins
t_stat, p_val = ttest_rel(df_delta['PBRM1_Δ'], df_delta['ARID1A_Δ'])
print(f"Paired t-test: t={t_stat:.3f}, p={p_val:.3e}")

#    (b) non-parametric Wilcoxon signed-rank
w_stat, p_w = wilcoxon(df_delta['PBRM1_Δ'], df_delta['ARID1A_Δ'])
print(f"Wilcoxon test: W={w_stat}, p={p_w:.3e}")

plt.hist(delta_pbrm1 - delta_arid1a, bins=15)
plt.title("Histogram of ΔPBRM1 – ΔARID1A")
plt.show()






abundance1 = abundance.T
# 3) build a “long” DataFrame for statsmodels
def build_long(class_df, ko_name):
    # 1) pull out the slice you want
    df = abundance1.loc[class_df.index, proteins].copy()

    # 2) annotate it
    df['Class'] = class_df['class'].values   # “high” vs “low”
    df['KO']    = ko_name                   # “PBRM1” or “ARID1A”

    # 3) give the index a real name, then reset it
    df.index.name = 'CellLine'
    df = df.reset_index()   # now there's a column called 'CellLine'

    # 4) melt, keeping 'CellLine', 'KO', 'Class' as id_vars
    long = df.melt(
        id_vars=['CellLine','KO','Class'],
        value_vars=proteins,
        var_name='Protein',
        value_name='Abundance'
    )
    return long

long_p = build_long(class_pbrm1,  'PBRM1')
long_a = build_long(class_arid1a, 'ARID1A')
long   = pd.concat([long_p, long_a], ignore_index=True)


# 4) loop over proteins, fit Abundance ~ KO * Class, grab interaction p-value
p_int = {}
for prot in proteins:
    sub   = long[long['Protein'] == prot]
    model = smf.ols('Abundance ~ C(KO)*C(Class)', data=sub).fit()

    # inspect what terms you got (uncomment if you need to debug)
    # print(model.pvalues.index.tolist())

    # find the single term that represents KO×Class
    interaction_terms = [
        name for name in model.pvalues.index
        if 'C(KO)' in name and 'C(Class)' in name
    ]
    if len(interaction_terms) != 1:
        raise RuntimeError(f"Could not uniquely identify interaction term for {prot}: {interaction_terms}")
    inter = interaction_terms[0]

    p_int[prot] = model.pvalues[inter]

p_raw = pd.Series(p_int, name='p_raw').sort_values()

# 5) correct for multiple testing (e.g. Benjamini–Hochberg FDR)
reject, p_fdr, _, _ = multipletests(p_raw, alpha=0.05, method='fdr_bh')
results = pd.DataFrame({
    'p_raw':      p_raw,
    'p_fdr':      pd.Series(p_fdr, index=p_raw.index),
    'significant': pd.Series(reject, index=p_raw.index)
})

print(results)