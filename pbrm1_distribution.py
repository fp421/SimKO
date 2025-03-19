import pandas as pd
import numpy as np
from scipy.stats import t
import matplotlib.pyplot as plt
import seaborn as sns
import random
import math
import pandas as pd
import numpy as np
from scipy.stats import t
import matplotlib.pyplot as plt
import seaborn as sns
from simko_func import simko



abundance = pd.read_csv('data/simko2_data/passport_prots.csv', index_col=0)
abundance.index = abundance.index.astype(str)

#created data frame of top 'median' and bottom 'low' cell lines based on abundance of the "ko_protein"
class_df = simko.get_classes_by_mean_abundance(ko_proteins=['PBRM1'], abundance=abundance, n=30)
class_df

high_abundance = class_df[class_df['class'] == 'high']['mean']
low_abundance = class_df[class_df['class'] == 'low']['mean']

pbrm1_abundance = abundance.loc['PBRM1']

plt.figure(figsize=(10, 6))

# Plot the entire distribution of PBRM1 abundance using a KDE curve
sns.kdeplot(pbrm1_abundance, color='gray', label='Distribution of PBRM1 Abundance', shade=True, alpha=0.3)


# Adding a vertical line at 1.76 to indicate the threshold for the low group
plt.axvline(x=1.76, color='blue', linestyle='--', label="Low Group (n=30)")

# Adding vertical lines to indicate the range for the high group
plt.axvline(x=3.40, color='red', linestyle='--', label="High Group (n=30)")
plt.axvline(x=3.5, color='red', linestyle='--')

# Add labels and title
plt.title('PBRM1 Protein Abundance Distribution Across Cell Lines')
plt.xlabel('PBRM1 Abundance')
plt.ylabel('Frequency')

# Add a legend
plt.legend()

# Show the plot
plt.show()


#making a table for the high and low cell line