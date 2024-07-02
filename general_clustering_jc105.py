# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
# %% load packages
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
import pandas as pd
from sklearn.metrics import silhouette_samples, silhouette_score
import scipy.cluster.hierarchy as sch


# %% functions
def neuron_preprocess(neurons):
    '''
    Parameters
    ----------
    neurons : matrix
        any kinds of the matrix, but vector should be in row.

    Returns
    -------
    neuron_norm : matrix
        same size as neurons, should be normalized in row (row from 0 to 1).
    neurons_pca : matrix
        should be same size as neurons.
    neurons_tsne : matrix
        should be two vectors only.

    '''
    neuron_norm = normalize(neurons)
    pca = PCA()
    neurons_pca = pca.fit_transform(neurons)
    tsne = TSNE(n_components=2, n_iter=300)
    neurons_tsne = tsne.fit_transform(neurons)
    return neuron_norm, neurons_pca, neurons_tsne
    


def avg_cluster(clusters, new_idx, combined_cellid, combined_avg, neurons):
    """
    Inputs:
    clusters: array 
    new_idx: tca index in original average 2p plot
    combined_cellid: cell id for average 2p plot
    combined_avg: average 2p plot
    neurons: tca matrix

    Outputs:
    indextable: combination of the index
    average_plot: average plot for each cluster, row--cluster #, col--time, z--condition 
    avg_plot_cluster: concatenation of the average plot
    neurons_avg_bygroup: sort the tca vectors by cluster
    """
    # Reorder combined_cellid based on new_idx
    new_idx_flattened = new_idx.flatten()
    combined_cellid_reorder = [combined_cellid[i] for i in new_idx_flattened]

    # Create indextable DataFrame
    data = {
        'cluster': clusters,
        'cellid': np.array(combined_cellid_reorder),
        'avgplotindex': new_idx_flattened,
        'tcaindex': np.arange(1, len(new_idx) + 1)
    }
    
    indextable = pd.DataFrame(data, columns=['cluster', 'cellid', 'avgplotindex', 'tcaindex'])

    # Initialize avg_plot_cluster and tca_bygroup
    unique_clusters = np.unique(clusters)
    avg_plot_cluster = {cl: [] for cl in unique_clusters}
    data_bygroup = {cl: [] for cl in unique_clusters}

    # Populate avg_plot_cluster and tca_bygroup
    for neuron in range(len(new_idx)):
        cluster_number = indextable['cluster'][neuron]
        average_plot_index = indextable['avgplotindex'][neuron]
        avg_plot_cluster[cluster_number].append(combined_avg[average_plot_index, :, :])
        data_bygroup[cluster_number].append(neurons[neuron, :])

    # Compute average_plot and neurons_avg_bygroup
    average_plot = np.zeros((len(unique_clusters), combined_avg.shape[1], combined_avg.shape[2]))
    neurons_avg_bygroup = np.zeros((len(unique_clusters), neurons.shape[1]))

    for cl in unique_clusters:
        average_plot[cl - 1, :, :] = np.mean(avg_plot_cluster[cl], axis=0)
        neurons_avg_bygroup[cl - 1, :] = np.mean(neurons_avg_bygroup[cl-1], axis=0)

    return indextable, average_plot, avg_plot_cluster, neurons_avg_bygroup

# Example usage:
# clusters = ...  # array of clusters
# new_idx = ...  # tca index in original average 2p plot
# combined_cellid = ...  # cell id for average 2p plot
# combined_avg = ...  # average 2p plot
# neurons = ...  # tca matrix

# indextable, average_plot, avg_plot_cluster, neurons_avg_bygroup = avg_cluster(clusters, new_idx, combined_cellid, combined_avg, neurons)


def overlay_cluster(full_dataset, clusters):
    num_groups = len(np.unique(clusters))
    colors = plt.cm.get_cmap('tab20', num_groups)  # Using tab20 colormap for better visualization
    edge_colors = plt.cm.get_cmap('tab10', 7)  # Using tab10 colormap for edge colors
    legend_entries = [f'Cluster {i}' for i in range(1, num_groups + 1)]

    fig = plt.figure()
    plt.title('Clustering Result Overlay On TCA')
    plt.xlabel('1st PCA Component')
    plt.ylabel('2nd PCA Component')
    plt.grid(True)
    plt.gca().tick_params(labelsize=14)
    
    if full_dataset.shape[1] > 2:
        ax = fig.add_subplot(111, projection='3d')
        ax.set_zlabel('3rd PCA Component')
        for i in range(1, num_groups + 1):
            idx = clusters == i
            ax.scatter(full_dataset[idx, 0], full_dataset[idx, 1], full_dataset[idx, 2],
                       color=colors(i - 1), edgecolor=edge_colors(i % 7), linewidth=1.5, label=f'Cluster {i}')
    else:
        for i in range(1, num_groups + 1):
            idx = clusters == i
            plt.scatter(full_dataset[idx, 0], full_dataset[idx, 1],
                        color=colors(i - 1), edgecolor=edge_colors(i % 7), linewidth=1.5, label=f'Cluster {i}')
    
    plt.legend(legend_entries, loc='best')
    plt.show()

# Example usage:
# full_dataset = np.random.rand(100, 3)  # Replace with your actual data
# clusters = np.random.randint(1, 6, size=100)  # Replace with your actual cluster assignments
# overlay_cluster(full_dataset, clusters)

def combined_cell_convert(combined_cellid):
    combined_cellid_flat = combined_cellid.flatten()
    combined_cellid_flat = [str(cell[0]) for cell in combined_cellid_flat]
    return combined_cellid_flat

def find_indices(combined_cellid, selected_signals):
    '''
    
    Parameters
    ----------
    combined_cellid : numpy.ndarray
        the combined cell ID.
    selected_signals : list
        the selected cell ID.

    Returns
    -------
    ind_array : list
        index of selected_signal ID in combined_cellid 
        (should be same len as selected_signal)
    indices : dict

    '''
    
    
    ind_array = []
    indices = {}
    for signal in selected_signals:
        try:
            index = combined_cellid.index(signal)
            indices[signal] = index
            ind_array.append(index)
        except ValueError:
            indices[signal] = None  # Signal not found in combined_cellid
            
    return ind_array, indices

def cluster_single_group(array):
    return np.ones(len(array))
    


# %%load data
data_path = '/Users/songyangwang/Desktop/Chen Lab/'
jc105_tca = scipy.io.loadmat(data_path+"jc105_tca.mat")
combined_average = scipy.io.loadmat(data_path+"combined_average.mat")

# Assuming 'neurons' is the key for the relevant data in jc105_tca
neurons = jc105_tca['neurons']
combined_avg = combined_average['combined_avg']
combined_cellid = combined_cell_convert(combined_average['combined_cellid'])
new_idx = jc105_tca['new_idx']-1 #python index vs matlab index
condition = combined_average['condition']
# %% Preprocess
[neuron_norm,neurons_pca,neurons_tsne] = neuron_preprocess(neurons)
#%%
# Visualize Data
neurons_norm_pca = normalize(neurons_pca)
overlay_cluster(neurons_norm_pca[:, :4], cluster_single_group(neurons_norm_pca))


# %% selected neurons
with open('selected_signal.txt', 'r') as file:
    selected_signals = file.read().splitlines()
    
[ind_array,indices] = find_indices(combined_cellid, selected_signals)
overlay_cluster(neurons_norm_pca[ind_array, :4], cluster_single_group(neurons_norm_pca[ind_array]))

# %% TODO
# Hierarchical Clustering
file_name = '/Users/songyangwang/Desktop/Chen Lab/Cluster_Result/hierarchical_normonly.pdf'

pair_dist_norm = pdist(neuron_norm)
cluster_tree_norm = linkage(neuron_norm, method="complete", metric="correlation")

plt.figure()
dendrogram(cluster_tree_norm, orientation='left')
plt.title('sample 25')
plt.show()

k = 8
clusters = sch.fcluster(cluster_tree_norm, k, criterion='maxclust')
silhouette_avg = silhouette_score(neuron_norm, clusters)
print(f"Silhouette score for {k} clusters: {silhouette_avg}")

indextable_h, average_plot_h, avg_plot_cluster_h,_ = avg_cluster(clusters, new_idx, combined_cellid, combined_avg, neurons)
exclude = []
plot_cluster(clusters, condition, average_plot_h, exclude, 1, file_name)
overlay_cluster(neurons_tsne, clusters)
overlay_cluster(neurons_norm_pca[ind_array, :4], clusters[ind_array])

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