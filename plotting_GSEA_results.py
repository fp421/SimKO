import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import openpyxl
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from simko_func import simko

simko_all_gsea = pd.read_csv('~/icr/simko/GSEA/PBRM1_exp_avg_GSEA_results/gseapy.prerank.gene_sets.report.csv')

simko_all_gsea['Term'] = simko_all_gsea['Term'].str.replace('_', ' ', regex=False)
simko_all_gsea['Term'] = simko_all_gsea['Term'].str.replace('hallmark', '', case=False, regex=False)
simko_all_gsea['Term'] = simko_all_gsea['Term'].str.strip()
simko_all_gsea

#boxplot for gsea hallmark set
colors2 = ["red" if x > 0 else "blue" for x in simko_all_gsea["nes"]]

plt.figure(figsize=(8, 10))
sns.barplot(data=simko_all_gsea, x="nes", y="Term", palette=colors2)
plt.xlabel("NES")
plt.title("TP53")
plt.xticks(rotation=45, ha="right")  # Rotate x-axis labels for better readability
plt.tight_layout()
plt.show()

#now trying only the pathways from the paper 
selected_pathways = ['ALLOGRAFT REJECTION', 'APOPTOSIS', 'BILE ACID METABOLISM',
                     'CHOLESTEROL HOMEOSTASIS', 'COAGULATION', 'EPITHELIAL MESENCHYMAL TRANSITION',
                     'FATTY ACID METABOLISM', 'HEME METABOLISM', 'IL2 STAT5 SIGNALING',
                     'IL6 JAK STAT3 SIGNALING', 'INTERFERON GAMMA RESPONSE', 'MYOGENESIS',
                     'TNFA SIGNALING VIA NFKB', 'XENOBIOTIC METABOLISM', 'DNA REPAIR', 
                     'E2F TARGETS', 'G2M CHECKPOINT', 'MYC TARGETS V1']

selected_pathways2 = ['OXIDATIVE PHOSPHORYLATION', 'MYOGENESIS', 'ESTROGEN RESPONSE EARLY', 'XENIBIOTIC METABLOSIMS',
                      'EPITHELIAL MESENCHYMAL TRANSITION', 'ESTROGEN RESPONSE LATE', 'E2F TARGETS', 'APOPTOSIS',
                      'MYC TARGETS V1', 'INTERFERON GAMMA RESPONSE', 'COAGULATION', 'KRAS SIGNALING DN',
                      'ADIPOGENESIS']

# Create a pattern that matches any of the selected pathways
pattern = '|'.join(selected_pathways)

# Filter the DataFrame: keep only rows where 'Term' contains one of the selected pathways (case-insensitive)
filtered_df = simko_all_gsea[simko_all_gsea['Term'].str.contains(pattern, case=False, na=False)]
filtered_df

#box of only these pathways
colors3 = ["red" if x > 0 else "blue" for x in filtered_df["nes"]]

highlight = {"DNA REPAIR", "E2F TARGETS", "G2M CHECKPOINT", "MYC TARGETS V1"}

plt.figure(figsize=(10, 8))
ax = sns.barplot(data=filtered_df, x="nes", y="Term", palette=colors3)
# Add black edges to each bar
for bar in ax.patches:
    bar.set_edgecolor('black')
    bar.set_linewidth(1)  # Adjust the linewidth as needed
# bold specific y-labels
for label in ax.get_yticklabels():
    if label.get_text() in highlight:
        label.set_fontweight('bold')
ax.tick_params(axis='y', labelsize=12) 
ax.tick_params(axis='x', labelsize=12) 
plt.axvline(0, color='black', linestyle='-', linewidth=2)
plt.xlabel("NES", fontsize=14)
plt.title("Experimental GSEA Results", fontsize=16)
plt.xticks()  # Rotate x-axis labels for better readability
plt.tight_layout()
plt.savefig("exp_gsea_poster.pdf", format='pdf')
plt.show()


#plotting pbrm1 simko gsea results - but in same order as experimental
simko_gsea = pd.read_csv('~/icr/simko/GSEA/PBRM1_simko_GSEA/gseapy.prerank.gene_sets.report.csv')

simko_gsea['Term'] = simko_gsea['Term'].str.replace('_', ' ', regex=False)
simko_gsea['Term'] = simko_gsea['Term'].str.replace('hallmark', '', case=False, regex=False)
simko_gsea['Term'] = simko_gsea['Term'].str.strip()
simko_gsea

pattern = '|'.join(selected_pathways)

# Filter the DataFrame: keep only rows where 'Term' contains one of the selected pathways (case-insensitive)
filtered_df2 = simko_gsea[simko_gsea['Term'].str.contains(pattern, case=False, na=False)]

# turn Term into an *ordered* categorical, with categories in the ref order
filtered_df2['Term'] = pd.Categorical(
    filtered_df2['Term'],
    categories=filtered_df['Term'].tolist(),
    ordered=True
)
# now a simple sort_values will respect that order
filtered_df2_ordered = (
    filtered_df2
    .sort_values('Term')
    .reset_index(drop=True)
)

#box of only these pathways
colors4 = ["red" if x > 0 else "blue" for x in filtered_df2_ordered["nes"]]
highlight = {"DNA REPAIR", "E2F TARGETS", "G2M CHECKPOINT", "MYC TARGETS V1"}

plt.figure(figsize=(10, 8))
ax = sns.barplot(data=filtered_df2_ordered, x="nes", y="Term", palette=colors4)
# Add black edges to each bar
for bar in ax.patches:
    bar.set_edgecolor('black')
    bar.set_linewidth(1)  # Adjust the linewidth as needed
# bold specific y-labels
for label in ax.get_yticklabels():
    if label.get_text() in highlight:
        label.set_fontweight('bold')
ax.tick_params(axis='y', labelsize=12) 
ax.tick_params(axis='x', labelsize=12) 
plt.axvline(0, color='black', linestyle='-', linewidth=2)
plt.xlabel("NES", fontsize=14)
plt.title("SimKO GSEA Results", fontsize=16)
plt.xticks()  # Rotate x-axis labels for better readability
plt.tight_layout()
plt.savefig("simko_gsea_poster.pdf", format='pdf')
plt.show()





#plotting the downregulaetd pathways and the log 1o pvalues
downreg_pathways = filtered_df[filtered_df['nes'] < 0]
# Calculate -log10(p-value)
downreg_pathways['log10_pval'] = -np.log10(downreg_pathways['pval'])
# Sort for better visual
downreg_pathways = downreg_pathways.sort_values('nes')
y_pos = np.arange(len(downreg_pathways))
# Create the figure
fig, ax = plt.subplots(figsize=(7, 3))
# Plot NES: black bars
ax.barh(y_pos, downreg_pathways['nes'], color='blue', label='NES (downregulated)', align='center')
# Plot -log10(p-value): white bars with black edges
ax.barh(y_pos, downreg_pathways['log10_pval'], color='black', linewidth=1.5, label='-log10(p-value)', align='center')
# Y-axis: pathway names
ax.set_yticks(y_pos)
ax.set_yticklabels(downreg_pathways['Term'])
# Add vertical line at 0
ax.axvline(0, color='grey', linewidth=1)
# Axis label and title
ax.set_xlabel('NES ←        → -log10(p-value)')
ax.set_title('Downregulated Pathways: NES and Significance')
# Add legend (key)
ax.legend(loc='upper right')
plt.tight_layout()
plt.show()