import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
from sklearn.impute import SimpleImputer
from sklearn.cluster import KMeans

abundance_data = pd.read_csv('~/icr/simko/data/simko2_data/passport_prots.csv', index_col=0)
abundance_data.index = abundance_data.index.astype(str)
abundance_data = abundance_data.T

imputer = SimpleImputer(strategy='mean')  # You can also try 'median'
abundance_imputed = imputer.fit_transform(abundance_data)
abundance_imputed = pd.DataFrame(abundance_imputed, columns=abundance_data.columns, index=abundance_data.index)
abundance_imputed

#trying clustering
kmeans = KMeans(n_clusters=17)
kmeans.fit(abundance_imputed)
abundance_imputed['cluster'] = kmeans.labels_

from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# Perform PCA for dimensionality reduction
pca = PCA(n_components=2)
pca_result = pca.fit_transform(abundance_imputed)

# Plot the clusters in 2D space
plt.figure()
plt.scatter(pca_result[:, 0], pca_result[:, 1], c=kmeans.labels_, cmap='viridis')
plt.title("PCA of Protein Abundance with K-Means Clusters")
plt.xlabel("Principal Component 1")
plt.ylabel("Principal Component 2")
plt.show()




protein_to_knockout = "PBRM1"
X = abundance_imputed
y = abundance_imputed.drop(columns=["PBRM1"])

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Train a multi-output regressor
from sklearn.multioutput import MultiOutputRegressor
from sklearn.ensemble import RandomForestRegressor

model = MultiOutputRegressor(RandomForestRegressor(random_state=42))
model.fit(X_train, y_train)

# Simulate a knockout by setting PBRM1 to 0
X_knockout = X_test.copy()
X_knockout["PBRM1"] = 0
y_simulated = model.predict(X_knockout)

