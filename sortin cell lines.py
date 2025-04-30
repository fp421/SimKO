import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import seaborn as sns
from scipy.stats import mannwhitneyu
import openpyxl
from simko_func import simko

abundancedata = pd.read_csv('~/icr/simko/data/simko2_data/passport_prots.csv', index_col=0)
abundancedata.index = abundancedata.index.astype(str)
abundancedata
#loading in cell line lists and tissue of origin
model_lists = pd.read_csv('~/icr/simko/data/simko2_data/model_list_20240110.csv', index_col=0)
model_lists['tissue'].value_counts()
model_lists
li_models = model_lists.loc[model_lists['tissue'] == 'Large Intestine', 'model_name']
li_models
abundancedata[abundancedata.columns.intersection(li_models)]


cls = pd.DataFrame({'model_name':abundancedata.columns})
cls = cls.merge(model_lists[['model_name', 'tissue']])
cls
cls['tissue'].value_counts().sum()



#need to make sure that the cell lines in model_list matched us with cell lines in abundance data
abundance_cell_lines = set(abundancedata.columns)
model_list_filt = model_lists[model_lists['model_name'].isin(abundance_cell_lines)]
unique_tissues = model_list_filt['tissue'].unique()
len(unique_tissues)
#column names to sue: model_name and tissues
# making a dictionary
grouped = model_list_filt.groupby('tissue')['model_name'].apply(list)
tissue_to_cell_lines = grouped.to_dict()
print(tissue_to_cell_lines)
#finding out how many cell lines for each tissue
for key, value in tissue_to_cell_lines.items():
    if isinstance(value, list):  # Check if the value is a list
        print(f"The length of the list for {key} is: {len(value)}")


specific_CL = tissue_to_cell_lines["Ovary"]
abundance_specific_CL = abundancedata[specific_CL]
abundance_specific_CL

#now we can do diff_stats on the specific cell lines
class_df_CL = simko.get_classes_by_mean_abundance(ko_proteins=['ARID1A'], abundance=abundance_specific_CL, n=10)
class_df_CL
control_diffs_CL = simko.get_control_differentials(abundance_specific_CL, class_df_CL, k=12)
diffs_CL = simko.get_ko_differentials(abundance=abundance_specific_CL, class_df=class_df_CL)
diffs_CL
PBRM1_row_CL = diffs_CL[diffs_CL["protein"] == "PBRM1"]
PBRM1_row_CL
diff_stats_CL = simko.get_significance(diffs_CL, control_diffs_CL, n=14)
#-np.log10() on pvals for volcano plot
diff_stats_CL['log10_pval']= -np.log10(diff_stats_CL['adjusted_p'])
diff_stats_CL['diff_copy'] = diff_stats_CL['diff']
diff_stats_CL = diff_stats_CL.sort_values(by='log10_pval', ascending = False)
diff_stats_CL[:15]

###if i want to use onyl significant values
#diff_stats_CL_sf = diff_stats_CL.loc[diff_stats_CL['adjusted_p'] < 0.05]
#diff_stats_CL_sf

#use function to assign pathways to proteins
diff_stats_cl_pathways = get_pathway_assigned(diff_stats_CL)

#boxplot
boxplots_for_pathways(diff_stats_cl_pathways, x='pathways', y='diff_copy',
                      title='Skin (n=55)',
                      figsize=(10,8),
                      ylim=(-8, 6),
                      xlabel=' ',
                      ylabel=' ',
                      save_path='skin_pathways')


#plotting distribution of the mdeiana dn low classes of proteins
class_df_CL
plt.figure(figsize=(10, 1))
sns.boxplot(class_df_CL, x='mean', hue='class', boxprops=dict(alpha=0.75), palette='hls')
#plt.title('Kidney', fontsize=18)
plt.xlim(-0.5, 4.5)
plt.ylabel(' ')
plt.xlabel('ARID1A Abundance Distribution')
plt.legend(title='Class (n=10)', loc='upper left')
#plt.savefig("skin_dist", dpi=300, bbox_inches='tight')
plt.show()

#statisitical test - is there a sig diference between median and low groups
class_df_CL
median_group = class_df_CL[class_df_CL['class'] == 'median']['mean']
low_group = class_df_CL[class_df_CL['class'] == 'low']['mean']

# Perform the Mann-Whitney U test
stat, p_value = mannwhitneyu(median_group, low_group, alternative='two-sided')

# Print the results
print("Mann-Whitney U Test Results:")
print(f"U-statistic: {stat}")
print(f"P-value: {p_value}")


#grid of bar plots - first the pathway changes
#first call path to the folder
CL_boxplot_folder = 'frans_fiigs/pathways_per_tissue'
#loading images from the folder
CL_boxplots = sorted([os.path.join(CL_boxplot_folder, img) for img in os.listdir(CL_boxplot_folder) if img.endswith(('.png'))])

#getting correct number of images - should be 13
num_images = len(CL_boxplots)
cols = 5
rows = (num_images + cols - 1) // cols

#creating the figure
fig, axes = plt.subplots(rows, cols, figsize=(30,20))
axes = axes.flatten()
#looping throught the images to plot them 
for i, img_path in enumerate(CL_boxplots):
    img = mpimg.imread(img_path)  # Read the image
    axes[i].imshow(img)
    axes[i].axis('off') 
#hiding any unsued axes    
for j in range(num_images, len(axes)):
    axes[j].axis('off')
fig.text(0.5, 0, 'Pathways', va='center', fontsize=40)  # x-axis label 
fig.text(-0.04, 0.5, 'Protein Abundance Changes (LogFC)', va='center', rotation='vertical', fontsize=40)  # y-axis label          
plt.tight_layout()
plt.show()



##grid of all the distribution between med and low for each cell line
CL_dist_folder = 'frans_fiigs/distributions'
#loading images from the folder
CL_dist_boxplots = sorted([os.path.join(CL_dist_folder, img) for img in os.listdir(CL_dist_folder) if img.endswith(('.png'))])

#getting correct number of images - should be 13
num_images1 = len(CL_dist_boxplots)
cols1 = 5
rows1 = (num_images1 + cols1 - 1) // cols1

#creating the figure
fig, axes = plt.subplots(rows1, cols1, figsize=(30,17))
axes = axes.flatten()
#looping throught the images to plot them 
for i, img_path in enumerate(CL_dist_boxplots):
    img = mpimg.imread(img_path)  # Read the image
    axes[i].imshow(img)
    axes[i].axis('off') 
#hiding any unsued axes    
for j in range(num_images1, len(axes)):
    axes[j].axis('off')

#adding x and y axis t the whole grid
fig.text(-0.04, 0.5, 'Mean Protein Abundances', va='center', rotation='vertical', fontsize=40)  # y-axis label   

plt.tight_layout()
plt.show()