import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import openpyxl

#8457 proteins
abundancedata1 = pd.read_csv('~/icr/simko/data/simko2_data/passport_prots.csv', index_col=0)
abundancedata1.index = abundancedata1.index.astype(str)
abundancedata1


#doing simko on the abundance data
#note: the p values will slightly change as the control shuffles are always different 
#this initial analysis is to only get the significant proteins
class_df = get_classes_by_mean_abundance(ko_proteins=['PBRM1'], abundance=abundancedata1, n=30)
control_diffs = get_control_differentials(abundancedata1, class_df, k=30)
diffs = get_ko_differentials(abundance=abundancedata1, class_df=class_df)
diffs
PBRM1_row = diffs[diffs["protein"] == "PBRM1"]
PBRM1_row
diff_stats = get_significance(diffs, control_diffs, n=30)
diff_stats['log10_pval']= -np.log10(diff_stats['adjusted_p'])
diff_stats['diff_copy'] = diff_stats['diff']
diff_stats = diff_stats.sort_values(by='log10_pval', ascending = False)
diff_stats
#just checking what cell lines are in median and low class
class_df

#now we only want the significant proteins in the abundance data1 set
sf_proteins = diff_stats.loc[diff_stats['adjusted_p'] < 0.05]['protein']
sf_proteins

#now we want to reduce the abundance dataset to only have the significant proteins
#2783 significant proteins
abundancedata1_sf = abundancedata1.loc[sf_proteins]
abundancedata1_sf

### the cell lines in med and low are the same in the whole and reduced data set
#class_df_sf --> dont need to sort pathway assignments
class_df_sf = get_classes_by_mean_abundance(ko_proteins=['PBRM1'], abundance=abundancedata1_sf, n=30)
control_diffs_sf = get_control_differentials(abundancedata1_sf, class_df_sf, k=30)
control_diffs_sf
diffs_sf = get_ko_differentials(abundance=abundancedata1_sf, class_df=class_df_sf)
diffs_sf

#need to assign pathways to control_diffs_sf and diffs_sf
control_diffs_sf_p = get_pathway_assigned(control_diffs_sf)
control_diffs_sf_p
diffs_sf_p = get_pathway_assigned(diffs_sf)

#double grouping - finding the average pathway diff in eahc contorl_group
control_diff_grouped = control_diffs_sf_p.groupby(['ko', 'pathways'])['diff'].mean().reset_index()
control_diff_grouped

diff_stats_2 = get_significance(diffs_sf_p, control_diff_grouped, n=30)
diff_stats_2
#now we need to fidn the mean for each pathwat and put the mean for each in a new table
#maybe just the get stats function needs to do the groupby "pathway" instead

#need to thibk about how to separate into pathways and then do stats that way 