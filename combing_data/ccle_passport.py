import pandas as pd
import numpy as np

passport = pd.read_csv('~/icr/simko/data/simko2_data/passport_prots.csv', index_col=0)
ccle = pd.read_csv('~/icr/simko/data/simko2_data/CCLE_all.csv', index_col=0)

#removing _ and tissue name from ccle data 
ccle.columns = [col.split('_')[0] for col in ccle.columns]
#removing - and capitalising CL names in passport data 
passport.columns = passport.columns.str.replace('-', '').str.upper()


#how many common cell lines are there
a = np.intersect1d(passport.columns, ccle.columns)


# Unique to passport
unique_to_passport = set(passport.columns) - set(a)
print(f"Cell lines unique to passport: {len(unique_to_passport)}")
# Unique to ccle
unique_to_ccle = set(ccle.columns) - set(a)
print(f"Cell lines unique to ccle: {len(unique_to_ccle)}")

#only looking at the shared cls for both data sets - how different are the protein values
#also only want to keep same proteins between data sets
common_proteins = passport.index.intersection(ccle.index) 
ccle_common = ccle.loc[common_proteins, a]
ccle_common = ccle_common.sort_index()
passport_common = passport.loc[common_proteins, a]

#checkikng for duplicates as lengths of dfs are different
print(ccle.index.duplicated().sum())
print(passport.index.duplicated().sum())

#removing proteins that have duplicate completely
duplicates = ccle_common.index[ccle_common.index.duplicated()].unique()
ccle_clean = ccle_common.drop(index=duplicates)
passport_clean = passport_common.drop(index=duplicates)

#also need to remove duplicated columns
dupe_cols = ccle_clean.columns[ccle_clean.columns.duplicated()].unique()
ccle_final = ccle_clean.drop(columns=dupe_cols)
passport_final = passport_clean.drop(columns=dupe_cols)

ccle_final
passport_final



# Calculate the absolute difference
difference = (passport_final - ccle_final).abs()

# Summary statistics per cell line
mean_diff_per_cellline = difference.mean()
print(mean_diff_per_cellline)