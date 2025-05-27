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


#from importlib import reload
#import simko.simko as simko
#reload(simko)

data_dir = '~/icr/simko/simko_git/simko2_analysis-main/data/simko2_data/'


model_annotation = pd.read_csv('data/simko2_data/model_list_20240110.csv')
protein_set = ['BRCA1', 'BRCA2', 'BARD1', 'SMC6', 'SMC5', 'ARID1A', 'PBRM1', 'SMARCA1', 'SMARCC1', 'STAG1', 'STAG2', 'SMARCE1', 'SMARCD1', 'PLK1']
abundance = pd.read_csv('data/simko2_data/passport_prots.csv', index_col=0)
abundance.index = abundance.index.astype(str)
abundance

#metric test - what % of cell lines is ko protein expressed in
total_cell_lines = abundance.shape[1]
prot_non_nans = abundance.loc['CBL'].count()
prot_non_nans
abundance_percentage = (prot_non_nans/total_cell_lines)*100
abundance_percentage


## Single SimKO
#created data frame of top 'median' and bottom 'low' cell lines based on abundance of the "ko_protein"
class_df = simko.get_classes_by_mean_abundance(ko_proteins=['PBRM1'], abundance=abundance, n=30)
# random abundance difference of proteisn between randomly selected CLs (30) - iterated through 100 combos
control_diffs = simko.get_control_differentials(abundance, class_df, k=30)
#box plot for proteins in 'protein set'  - plots abundance difference per protein fro each of the trials
sns.boxplot(control_diffs.loc[control_diffs['protein'].isin(protein_set)], y='protein', x='diff')
#computes difference between mean proten abundance 
diffs = simko.get_ko_differentials(abundance=abundance, class_df=class_df)
#is there significant difference between observed diffs (from KO sim) and differences due to chance?
diff_stats = simko.get_significance(diffs, control_diffs, n=30)
diff_stats.sort_values(by = 'adjusted_p', ascending=True).head(20)





#creating a plot of the distributiom of protein across cell lines
protein_of_interest = 'ARID1A'  # change this to the protein you want
subset_df = control_diffs[control_diffs['protein'] == protein_of_interest]
subset_df

#finding the observed difference in a knockout
subset_df_ko = diff_stats[diff_stats['protein'] == protein_of_interest]
subset_df_ko
diff = subset_df_ko['diff'].values[0]
diff

# Create the KDE plot
sns.kdeplot(data=subset_df, x='diff', fill=True, color='teal', label='Control Distribution')
plt.xlabel('Difference')
plt.ylabel('Density')
plt.axvline(diff, color='lightcoral', linestyle='--', label='Observed diff (KO)')
plt.legend()
plt.title('Control Distribution of Abundance Changes for ARID1A')
plt.show()


#doing the above code but for a group of proteins and plotting it all together
proteins = ['BRCA1','BRCA2','SMC6','SMC5',
            'ARID1A','SMARCA1','SMARCC1','STAG1',
            'STAG2','SMARCE1','SMARCD1','PLK1']

n = len(proteins)
# make the figure very “short” per protein
fig, axes = plt.subplots(n, 1,
                         figsize=(6, 0.45 * n),
                         sharex=True, sharey=True)

for ax, prot in zip(axes, proteins):
    # 1) control distribution
    df_ctrl = control_diffs[control_diffs['protein'] == prot]
    sns.kdeplot(data=df_ctrl, x='diff',
                fill=True, alpha=0.6,
                color='teal', ax=ax)
    # 2) observed KO diff (single value)
    diff_val = diff_stats.loc[
        diff_stats['protein'] == prot, 'diff'
    ].item()
    ax.axvline(diff_val,
               color='lightcoral',
               linestyle='--',
               linewidth=1.5)
    ax.axhline(0,
               color='black',
               linestyle='--',
               linewidth=1,
               zorder=0)
    ax.tick_params(axis='x', which='both', length=0)
    # 3) kill the spines / box
    for spine in ax.spines.values():
        spine.set_visible(False)
    # 4) no y‐ticks, no margins
    ax.set_ylabel('')
    ax.set_yticks([])
    ax.margins(y=0)
    # 5) protein label on left & right
    ax.text(-0.02, 0.5, prot,
            transform=ax.transAxes,
            ha='right', va='center',
            fontsize=12)
# global labels
axes[-1].set_xlabel('LogFC')
fig.supylabel('')  # density axis is implicit
# create handles
plt.tight_layout()
plt.subplots_adjust(left=0.15, right=0.85, top=0.98, bottom=0.12)
plt.savefig("null_dist_plot_poster.pdf", format='pdf')
plt.show()

diff_for_selected_prots = diff_stats.loc[diff_stats['protein'].isin(proteins)]
diff_for_selected_prots[['protein', 'adjusted_p']]

#invetsigating baf and pbaf proteins for arid1a/pbrm1

BAF = ('ARID1A', 'SMARCC1', 'SMARCC2', 'SMARCE1', 'SMARCB1', 'SMARCD1', 'SMARCD2',
       'SMARCD3', 'DPF1', 'DPF2', 'DPF3', 'SMARCA2', 'SMARCA4', 'SS18', 'ACTB',
       'ACTL6A', 'BCL7A', 'BCL7B', 'BCL7C')

PBAF = ('ARID2', 'PHF10', 'BRD7', 'PBRM1', 'SMARCC1', 'SMARCC2', 'SMARCE1', 'SMARCB1',
        'SMARCD1', 'SMARCD2', 'SMARCD3', 'SMARCA2', 'SMARCA4', 'BCL7A',
        'BCL7B', 'BCL7C', 'ACTB', 'ACTL6A')

baf_results = diff_stats[diff_stats['protein'].isin(BAF)]
baf_results = baf_results.sort_values(by='diff', ascending=True)
baf_results
threshold = 0.05
colors = ['tab:blue' if p < threshold else 'grey' for p in baf_results['adjusted_p']]
legend_elements = [
    Patch(facecolor='tab:blue', edgecolor='black', label='Significant (p < 0.05)'),
    Patch(facecolor='grey', edgecolor='black', label='Not significant')
]
#plot in a bar plot
plt.figure(figsize=(10, 6))
plt.bar(baf_results['protein'], baf_results['diff'], color=colors, edgecolor='black')
plt.xticks(rotation=90)
plt.xlabel(' ')
plt.ylabel('Abundance Change (LogFC)', fontsize=12)
#plt.ylim(-2.5, 0.5)
plt.axhline(y=0, color='black', linestyle='-')
plt.title('Change in BAF Protein Abundance after ARID1A Knockdown', fontsize=14)
plt.legend(handles=legend_elements, loc='upper left')
plt.tight_layout()
plt.show()

pbaf_results = diff_stats[diff_stats['protein'].isin(PBAF)]
pbaf_results = pbaf_results.sort_values(by='diff', ascending=True)
pbaf_results
colors1 = ['tab:blue' if p < threshold else 'grey' for p in pbaf_results['adjusted_p']]

# Create custom legend
legend_elements = [
    Patch(facecolor='tab:blue', edgecolor='black', label='Significant (p < 0.05)'),
    Patch(facecolor='grey', edgecolor='black', label='Not significant')
]

plt.figure(figsize=(10, 6))
plt.bar(pbaf_results['protein'], pbaf_results['diff'], color=colors1, edgecolor='black')
plt.xticks(rotation=90)
plt.xlabel(' ')
plt.ylabel('Abundance Change (LogFC)', fontsize=12)
plt.axhline(y=0, color='black', linestyle='-')
plt.title('Change in PBAF Protein Abundances after PBRM1 Knockdown', fontsize=14)
plt.legend(handles=legend_elements, loc='upper left')
plt.tight_layout()
plt.show()


#plot that compares common protein in the complexes but in different KOs (arid1a vs pbrm1)
common = baf_results.merge(pbaf_results, on='protein', suffixes=('_BAF', '_PBAF'))
# Sort by BAF diff or any order you prefer
common_sorted = common.sort_values(by='diff_BAF')
proteins = common_sorted['protein'].to_list()
x = np.arange(len(proteins))  # positions for each protein

width = 0.35  # width of each bar

plt.figure(figsize=(12, 6))
plt.bar(x - width/2, common_sorted['diff_BAF'], width, label='ARID1A Knockout', color='steelblue', edgecolor='black')
plt.bar(x + width/2, common_sorted['diff_PBAF'], width, label='PBRM1 Knockout', color='orange', edgecolor='black')
plt.xticks(x, proteins, rotation=45)
plt.xlabel('Common BAF/BAF Proteins', fontsize=14)
plt.ylabel('LogFC', fontsize=14)
plt.axhline(y=0, color='black', linestyle='--')
plt.title('Comparison of Predicted PBAF/BAF Protein Changes: ARID1A vs PBRM1 Knockout', fontsize=16)
plt.legend(fontsize=14)
plt.ylim(-1.5, 0.6)  # Adjust based on your data
plt.tight_layout()
plt.savefig("baf_pbaf_poster.pdf", format='pdf')
plt.show()




#getting a graph/diagram to visualise high and low cell lines 
class_df
sorted_df = class_df.sort_values(by='mean')
cell_lines = sorted_df.index
mean_protein = sorted_df['mean'].values

# Create a figure with a wide layout to represent a single row of squares
fig, ax = plt.subplots(figsize=(20, 2))

# Plot a single-row heatmap.
# By wrapping mean_protein in a list, we create a 2D array with one row.
# The 'bwr' colormap will map low values to blue and high values to red.
heatmap = ax.imshow([mean_protein], cmap='bwr', aspect='auto')

# Set x-axis ticks to correspond to each cell line and label them accordingly
ax.set_xticks(np.arange(len(cell_lines)))
ax.set_xticklabels(cell_lines, rotation=90, ha='center')
ax.set_yticks([])  # Remove the y-axis ticks since we only have one row

# Add the numeric mean values inside each square, rotated 90° and in black.
for i, value in enumerate(mean_protein):
    ax.text(i, 0, f"{value:.2f}", ha='center', va='center', 
            color='black', fontsize=8, rotation=90)

plt.tight_layout()
plt.savefig('cl_heatmap.pdf', format='pdf')
plt.show()







#getting associated proteins to the 'knockout protein 
#prot_assoc = get_associated_proteins(ko_protein='ARID1A')
#prot_assoc

#brca1_results = diff_stats[diff_stats['protein'].isin(prot_assoc)]
#brca1_results


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