import pandas as pd
import numpy as np
from scipy.stats import t, shapiro, kstest
from sklearn.preprocessing import MinMaxScaler, StandardScaler
import matplotlib.pyplot as plt

def prepare_data():
    #loading data
    abundance = pd.read_csv('~/icr/simko/data/simko2_data/passport_prots.csv', index_col=0)
    abundance.index = abundance.index.astype(str)

    #removing cell lines with over 4000 nans
    nans_per_cl = abundance.isna().sum(axis=0)
    abundance_cl_filtered = abundance.loc[:, nans_per_cl<4000]

    #getting rid of protein with over 80% NaN (from the dataset filtered by CLs)
    prot_nan_count = abundance_cl_filtered.isna().sum(axis=1)
    prot_nan_percent = (prot_nan_count/abundance_cl_filtered.shape[1])*100

    abundance_filtered = abundance_cl_filtered[prot_nan_percent<80]

    #imputing witht the lower quartile average for each protein
    #set the protein names as the index - ignores it while we find the lower quartile

    def average_lower_quartile(x):
        sorted_abundances = x.dropna().sort_values()
        lower_qt_values = sorted_abundances.iloc[:int(len(sorted_abundances) * 0.25)]
        return lower_qt_values.mean()


    lower_qt_averages = abundance_filtered.apply(average_lower_quartile, axis=1)

    abundance_filtered_no_nan = abundance_filtered.apply(lambda x: x.fillna(lower_qt_averages[x.name]), axis=1)

    #transposing
    abundance_imputed = abundance_filtered_no_nan.T

    #scaling the imputed data
    scaler = StandardScaler()
    scaled_data = pd.DataFrame(scaler.fit_transform(abundance_imputed), index=abundance_imputed.index, columns=abundance_imputed.columns)

    return scaled_data

scaled_data = prepare_data()