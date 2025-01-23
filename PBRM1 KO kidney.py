import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import openpyxl

PBRM1_data = pd.read_excel('~/icr/all_pbrm1_results.xlsx')

PBRM1_data 

#selecting the CEN proteins

CEN_proteins = ['MIS18a','MIS18BP1','RSF1','KAT7','SMC2','SMC4','NCAPD3','NCAPG2','NCAPH2',
                'CENPC','CENPT','CENPS','CENPX','CENPO','CENPU','CENPL','CENPN','CENPH','CENPI','CENPK',
                'KNL1','ZWINT','NDC80','NUF2','SPC24','SPC25','DSN1','MIS12','NSL1','PMF1',
                'INCENP','BIRC5','CDCA8','AURKB',
                'CENPE','CENPF','CENPB','CENPV','CENPJ','CDCA2',
                'CBX1','SIRT6','UHRF2','SUV39H1','SETDB1','DAXX','SSRP1','SMARCAD1',
                'ATRX','BAZ1A','BAZ1B','CBX3','CBX5','CHRAC1','DNMT1','EZH2','HELLS','KDM4A','KMT5C','LRWD1','MACROH2A1','POLE3','SMARCA5']

#abundance data
#will only be using kidney cancer cell lines

#kidney_cls = ['769-P','786-0','A498','A704','ACHN','BB65-RCC','BFTC-909','CAKI-1','CAL-54',
              #'HA7-RCC','KMRC-1','KMRC-20','LB1047-RCC','LB2241-RCC','LB996-RCC',
              #'NCC010','NCC021','OS-RC-2','RCC10RGB','RCC-AB','RCC-ER','RCC-FG2','RCC-JF','RCC-JW','RCC-MF',
              #'RXF393','SN12C','SW156','TK10','U031','VMRC-RCW','VMRC-RCZ']

abundance = pd.read_csv('~/icr/simko/data/simko2_data/passport_prots.csv', index_col=0)
abundance.index = abundance.index.astype(str)
abundance_kid = abundance[kidney_cls]

#creating the diffs data frame
class_df_k = get_classes_by_mean_abundance(ko_proteins=['PBRM1'], abundance=abundance_kid, n=8)
control_diffs_k = get_control_differentials(abundance_kid, class_df_k, k=8)
diffs_k = get_ko_differentials(abundance=abundance_kid, class_df=class_df_k)
diffs_k
PBRM1_row_k = diffs_k[diffs_k["protein"] == "PBRM1"]
PBRM1_row_k
diff_stats_k = get_significance(diffs_k, control_diffs_k, n=10)
#-np.log10() on pvals for volcano plot
diff_stats_k['log10_pval']= -np.log10(diff_stats_k['adjusted_p'])
diff_stats_k['diff_copy'] = diff_stats_k['diff']
diff_stats_k = diff_stats_k.sort_values(by='log10_pval', ascending = False)
diff_stats_k[:10]

#divide up data frame based on what pathways the proteins are a part of
diff_stats_k['pathways'] = diff_stats_k['protein'].apply(lambda protein: [
    pathway
    for pathway, proteins in {
        "allograft_rejection": allograft_rejection,
        "apoptosis": apoptosis,
        "bile_acid_metabolism": bile_acid_metabolism,
        "cholesterol_homeostasis": cholesterol_homeostasis,
        "coagulation": coagulation,
        "emt": emt,
        "fatty_acid_metabolism": fatty_acid_metabolism,
        "heme_metabolism": heme_metabolism,
        "il2_stat5_signalling": il2_stat5_signalling,
        "interferon_gamma_response": interferon_gamma_response,
        "myogenesis": myogenesis,
        "tnfa_via_nfkb": tnfa_via_nfkb,
        "xenobiotic_metabolism": xenobiotic_metabolism,
        "dna_repair": dna_repair,
        "e2f_targets": e2f_targets,
        "g2m_checkpoint": g2m_checkpoint,
        "myc_targets": myc_targets
    }.items()
    if protein in proteins
])

#duplicating/tripling etc protein rows that appear in multiple pathways
diff_stats_k = diff_stats_k.explode("pathways").reset_index(drop=True)
#getting rid of nan rows
diff_stats_k = diff_stats_k.dropna(subset=["pathways"]).reset_index(drop=True)
diff_stats_k

diff_stats_sf_k = diff_stats_k.loc[diff_stats_k['p'] < 0.05]
diff_stats_sf_k

#boxplot
plt.figure(figsize=(15,10))
sns.boxplot(data=diff_stats_sf_k, x='pathways', y='diff_copy')
plt.title("Boxplot of Protein Abundance Changes by Pathway - signficant proteins", fontsize=14)
plt.xlabel("Pathway", fontsize=12)
plt.ylabel("Abudnance Change (LogFC)", fontsize=12)
plt.xticks(rotation=45, ha="right")
plt.axhline(0, color='black', linestyle='--', linewidth=2)
plt.tight_layout()
plt.show()












#diffs from the pbrm1 data set
PBRM1_data.rename(columns={PBRM1_data.columns[11]: "kidney_logfc"}, inplace=True)
PBRM1_data.rename(columns={PBRM1_data.columns[2]: "protein_name"}, inplace=True)
PBRM1_data

#scatter plot comparing the diffs
common_proteins = set(diff_stats_k['protein']).intersection(set(PBRM1_data['protein_name'])) #7106 proteins

# Extract the "diff" values for the common proteins in both datasets
simulated_diff = diff_stats_k.loc[diff_stats_k['protein'].isin(common_proteins), 
                                ['protein', 'diff']].set_index('protein')
pbrm1_diff = PBRM1_data.loc[PBRM1_data['protein_name'].isin(common_proteins),
                             ['protein_name', 'kidney_logfc']]
pbrm1_diff = pbrm1_diff.rename(columns={'protein_name': 'protein', 'kidney_logfc': 'diff'}).set_index('protein')

# Combine into a single DataFrame for comparison
comparison_df = simulated_diff.join(pbrm1_diff, lsuffix='_simulated', rsuffix='_pbrm1')

comparison_df

print(comparison_df.index.duplicated().sum())  # Number of duplicate index entries
print(comparison_df[comparison_df.index.duplicated()])  # Show rows with duplicate indexes

comparison_df = comparison_df[~comparison_df.index.duplicated(keep='first')]
comparison_df = comparison_df.loc[comparison_df.index.isin(CEN_proteins)]
comparison_df

plt.figure(figsize=(8, 6))
sns.scatterplot(
    x=comparison_df['diff_simulated'], 
    y=comparison_df['diff_pbrm1']
)

# Add plot labels
plt.title('Comparison of Diff Values Between Simulated and PBRM1 Data')
plt.xlabel('Diff (Simulated Data)')
plt.ylabel('Diff (PBRM1 Data)')
plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
plt.axvline(0, color='gray', linestyle='--', linewidth=0.8)
plt.grid(True)

# Show the plot
plt.show()