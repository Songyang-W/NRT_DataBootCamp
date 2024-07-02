# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_samples, silhouette_score

# Load data
jc105_tca = scipy.io.loadmat("jc105_tca.mat")
combined_average = scipy.io.loadmat("combined_average.mat")

# Assuming 'neurons' is the key for the relevant data in jc105_tca
neurons = jc105_tca['neurons']
combined_avg = combined_average['combined_avg']
combined_cellid = combined_average['combined_cellid']
new_idx = combined_average['new_idx']
condition = combined_average['condition']

# Preprocessing
neuron_norm = normalize(neurons)
pca = PCA()
neurons_pca = pca.fit_transform(neurons)
neurons_norm_pca = normalize(neurons_pca)
tsne = TSNE(n_components=2, n_iter=300)
neurons_tsne = tsne.fit_transform(neurons)

# Visualization Functions
def overlay_cluster(data, clusters):
    plt.scatter(data[:, 0], data[:, 1], c=clusters, cmap='viridis')
    plt.show()

def plot_cluster(data, condition, average_plot, exclude, flag, file_name):
    # Implement the function based on the MATLAB code functionality
    pass

# Visualize Data
cluster_single_group = np.ones(len(neurons_norm_pca))
overlay_cluster(neurons_norm_pca[:, :4], cluster_single_group)

plt.figure()
plt.scatter(neurons_tsne[:, 0], neurons_tsne[:, 1])
plt.show()

# Hierarchical Clustering
file_name = '/Users/songyangwang/Desktop/Chen Lab/Cluster_Result/hierarchical_normonly.pdf'

pair_dist_norm = pdist(neuron_norm)
cluster_tree_norm = linkage(neuron_norm, method="complete", metric="correlation")

plt.figure()
dendrogram(cluster_tree_norm, orientation='left')
plt.title('sample 25')
plt.show()

k = 8
clusters = fcluster(cluster_tree_norm, k, criterion='maxclust')
silhouette_avg = silhouette_score(neuron_norm, clusters)
print(f"Silhouette score for {k} clusters: {silhouette_avg}")

# Assuming the Avg_Cluster function returns these values
indextable_h, average_plot_h, avg_plot_cluster_h = Avg_Cluster(clusters, new_idx, combined_cellid, combined_avg, neurons)
exclude = []
plot_cluster(clusters, condition, average_plot_h, exclude, 1, file_name)
overlay_cluster(neurons_tsne, clusters)
overlay_cluster(neurons_norm_pca[:, :4], clusters)

# K-means Clustering
opts = {'max_iter': 500}
max_clusters = 20
sumd_array = []

for k in range(1, max_clusters + 1):
    kmeans = KMeans(n_clusters=k, max_iter=500)
    kmeans.fit(neuron_norm)
    sumd_array.append(kmeans.inertia_)

plt.figure()
plt.subplot(211)
plt.plot(range(1, max_clusters + 1), sumd_array, '-o')
plt.xlabel('Number of Clusters k')
plt.ylabel('Sum of Squared Distances')
plt.title('Elbow Method for Optimal k')
plt.grid(True)

silhouette_values = []

for k in range(2, max_clusters + 1):
    kmeans = KMeans(n_clusters=k, max_iter=500)
    cluster_labels = kmeans.fit_predict(neuron_norm)
    silhouette_avg = silhouette_score(neuron_norm, cluster_labels)
    silhouette_values.append(silhouette_avg)

plt.subplot(212)
plt.plot(range(2, max_clusters + 1), silhouette_values, '-o')
plt.xlabel('Number of Clusters k')
plt.ylabel('Average Silhouette Score')
plt.title('Silhouette Analysis for Optimal k')
plt.grid(True)
plt.show()

file_name_kmean = '/Users/songyangwang/Desktop/Chen Lab/Cluster_Result/kmean_normonly.pdf'

kmeans = KMeans(n_clusters=9, max_iter=500)
idx = kmeans.fit_predict(neuron_norm)
silhouette_avg = silhouette_score(neurons_pca, idx)
print(f"Silhouette score for 9 clusters: {silhouette_avg}")

# Assuming the Avg_Cluster function returns these values
indextable_k, average_plot_k, avg_plot_cluster_k = Avg_Cluster(idx, new_idx, combined_cellid, combined_avg, neurons)
exclude = [1, 3, 4, 6, 2, 8, 7]
plot_cluster(idx, condition, average_plot_k, exclude, 1, file_name_kmean)
overlay_cluster(neurons_tsne, idx)
overlay_cluster(neurons_norm_pca, idx)

# Gaussian Mixture Model (GMM)
rng = np.random.default_rng(3)
k = 10
gmm = GaussianMixture(n_components=k, max_iter=1000, random_state=3)
gmm.fit(neurons)
cluster_idx = gmm.predict(neurons)
aic = gmm.aic(neurons)
bic = gmm.bic(neurons)

max_clusters = 10
aic_scores = []
bic_scores = []

for k in range(1, max_clusters + 1):
    gmm = GaussianMixture(n_components=k, max_iter=1000, random_state=3, reg_covar=0.01)
    gmm.fit(neurons)
    aic_scores.append(gmm.aic(neurons))
    bic_scores.append(gmm.bic(neurons))

cluster_idx = gmm.predict(neurons)

# Assuming the Avg_Cluster function returns these values
indextable_k, average_plot_k, avg_plot_cluster_k = Avg_Cluster(idx, new_idx, combined_cellid, combined_avg, neurons)
exclude = []
plot_cluster(idx, condition, average_plot_k, exclude, 2)
plt.figure()
plt.plot(range(1, max_clusters + 1), aic_scores, 'b-', range(1, max_clusters + 1), bic_scores, 'r-')
plt.legend(['AIC', 'BIC'])
plt.title('AIC and BIC for Different Numbers of Clusters')
plt.show()