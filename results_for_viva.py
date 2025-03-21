import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import seaborn as sns
from scipy.stats import mannwhitneyu
import openpyxl
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from simko_func import simko

#this code will perform simko on each tissue - and we will try to find if there are some proteins that always decrease with pbrm1

data = pd.read_csv('~/icr/simko/data/simko2_data/passport_prots.csv')
annotation = pd.read_csv('~/icr/simko/data/simko2_data/model_list_20240110.csv')

#creating a dictionary of tissue and the cell lines withn them
#need to make sure that the cell lines in model_list matched us with cell lines in abundance data
abundance_cell_lines = set(data.columns)
model_list_filt = annotation[annotation['model_name'].isin(abundance_cell_lines)]
unique_tissues = model_list_filt['tissue'].unique()
unique_tissues
#column names to sue: model_name and tissues
# making a dictionary
grouped = model_list_filt.groupby('tissue')['model_name'].apply(list)
tissue_to_cell_lines = grouped.to_dict()
print(tissue_to_cell_lines)
#finding out how many cell lines for each tissue
for key, value in tissue_to_cell_lines.items():
    if isinstance(value, list):  # Check if the value is a list
        print(f"The length of the list for {key} is: {len(value)}")

tissue_n_mapping = {
    'Adrenal Gland': 1,
    'Biliary Tract': 2,
    'Bladder': 4,
    'Bone': 10,
    'Breats': 12,
    'Central Nervous System': 14,
    'Cervix': 3,
    'Endometrium': 3,
    'Esophagus': 9,
    'Haematopoietic and Lymphoid': 30,
    'Head and Neck': 10,
    'Kidney': 8,
    'Large Intestine': 12,
    'Liver': 3,
    'Lung': 30,
    'Ovary': 10,
    'Pancreas': 8,
    'Peripheral Nervous System': 8,
    'Placenta': 1,
    'Prostate': 2,
    'Skin': 14,
    'Small Intestine': 1,
    'Soft Tissue': 4,
    'Stomach': 7,
}

results_dict = {}

#looping over each tissue to perform simko on their corresponding cell lines
# Loop over each tissue and its corresponding cell lines
for tissue, cell_lines in tissue_to_cell_lines.items():
    try:
        # Retrieve tissue-specific n, defaulting to 14 if not specified
        n_value = tissue_n_mapping.get(tissue, 14) 

        # Check that cell lines exist in your abundance data
        common_cell_lines = [cl for cl in cell_lines if cl in data.columns]
        if not common_cell_lines:
            print(f"No matching cell lines for tissue {tissue}")
            continue

        # Subset abundance data for the current tissue
        abundance_specific = data[cell_lines]
        
        # Run the analysis steps (replace with your function calls)
        class_df = simko.get_classes_by_mean_abundance(ko_proteins=['PBRM1'], abundance=abundance_specific, n=n_value)
        control_diffs = simko.get_control_differentials(abundance_specific, class_df, k=n_value)
        diffs = simko.get_ko_differentials(abundance=abundance_specific, class_df=class_df)
        diff_stats = simko.get_significance(diffs, control_diffs, n=n_value)
        
         # Check if diff_stats is a valid DataFrame
        if not isinstance(diff_stats, pd.DataFrame) or diff_stats.empty:
            print(f"No valid diff_stats for tissue {tissue}")
            continue

        # Calculate the additional columns
        diff_stats['diff_copy'] = diff_stats['diff']
        # Optionally sort (here sorting by significance)
        # Filter for significant proteins (e.g., adjusted_p < 0.05)
        diff_stats_sf = diff_stats.loc[diff_stats['adjusted_p'] < 0.05]

        # Filter for significant proteins (adjusted_p < 0.05)
        diff_stats_sf = diff_stats.loc[diff_stats['adjusted_p'] < 0.05]
        if diff_stats_sf.empty:
            print(f"No significant proteins for tissue {tissue}")
            continue        
        
        # Check if 'protein' column exists before setting the index
        if 'protein' not in diff_stats_sf.columns:
            print(f"No 'protein' column in diff_stats for tissue {tissue}")
            continue

        # Set the protein names as the index so that we can later merge by protein
        diff_stats_sf = diff_stats_sf.set_index('protein')
        
        # Save the diff_copy column for this tissue in the dictionary
        results_dict[tissue] = diff_stats_sf['diff_copy']
        
    except Exception as e:
        print(f"Error processing tissue {tissue}: {e}")

# Combine the results into a DataFrame where the index is protein names and each column is a tissue.
combined_df = pd.DataFrame(results_dict)

# Optionally, sort the DataFrame by protein names
combined_df = combined_df.sort_index()

# Display the combined DataFrame
print(combined_df)