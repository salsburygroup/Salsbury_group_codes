import argparse
from AH_clustering import clustering
import numpy as np
import mdtraj as md
import sklearn.metrics
import sklearn.cluster
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('bmh')
import copy
from collections import Counter

parser = argparse.ArgumentParser(description='Compare Minkowski weights', add_help=False)

inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('trajectory', type=str, help='Path to the trajectory file.')
inputs.add_argument('topology', type=str, help='Path to the topology file.')
inputs.add_argument('prefix', type=str, help='Prefix for output file names.')

UserInput = parser.parse_args()

trajectory = UserInput.trajectory
topology = UserInput.topology
prefix = UserInput.prefix

t = md.load(trajectory, top=topology)
sel = t.topology.select('not element H')
t = t.atom_slice(sel)

temp = t.xyz
frames = t.xyz.shape[0]
atoms = t.xyz.shape[1]
original_data = temp.reshape((frames, atoms * 3))
original_data = original_data.astype('float64')
temp = []

np.seterr(all='raise')
cl = clustering.Clustering()
if frames > 10000:
    sample_size = 10000
else:
    sample_size = None

original_data = cl.my_math.standardize(original_data)

optimal_p = 2

[labels, centroids, weights, ite, dist_tmp] = cl.imwk_means(original_data, p=optimal_p)
maxk = max(labels) + 1
silhouette_averages = np.zeros(maxk - 1)
print("Rescale iMWK trials")
k_to_try = np.arange(2, maxk + 1)
for k in k_to_try:
    print('Testing k=' + str(k) + ' of ' + str(maxk))
    cl = []
    labels = []
    weights = []
    centroids = []
    cl = clustering.Clustering()
    data = copy.copy(original_data)
    [labels, centroids, weights, ite, dist_tmp] = cl.imwk_means(data, p=optimal_p, k=k)

    for k1 in np.arange(0, max(labels) + 1):
        data[labels == k1] = np.multiply(data[labels == k1], np.tile(weights[k1], (np.sum(labels == k1), 1)))
        centroids[k1] = np.multiply(centroids[k1], weights[k1])

    kmeans_clusterer = sklearn.cluster.KMeans(n_clusters=k, n_init=100)
    kmeans_clusters = kmeans_clusterer.fit(data)
    labels = kmeans_clusters.labels_
    centroids = kmeans_clusters.cluster_centers_
    silhouette_averages[k - 2] = sklearn.metrics.silhouette_score(data, labels, sample_size=sample_size)

optimal_k = k_to_try[np.argmax(silhouette_averages)]

cl = []
labels = []
weights = []
centroids = []
cl = clustering.Clustering()
data = copy.copy(original_data)
[labels, centroids, weights, ite, dist_tmp] = cl.imwk_means(data, p=optimal_p, k=optimal_k)

for k1 in np.arange(0, max(labels) + 1):
    data[labels == k1] = np.multiply(data[labels == k1], np.tile(weights[k1], (np.sum(labels == k1), 1)))
    centroids[k1] = np.multiply(centroids[k1], weights[k1])

kmeans_clusterer = sklearn.cluster.KMeans(n_clusters=optimal_k, n_init=100)
kmeans_clusters = kmeans_clusterer.fit(data)
labels = kmeans_clusters.labels_
centroids = kmeans_clusters.cluster_centers_
silhouette_score = sklearn.metrics.silhouette_score(data, labels, sample_size=sample_size)

with open(prefix + '_AH_frame_clusters.txt', 'w') as f:
    for i, label in enumerate(labels):
        f.write(f'Frame {i}: Cluster {label}\n')

with open(prefix + '_AH_silhouette_score.txt', 'w') as f:
    f.write("Silhouette score: {0}\nOptimal p: {1}\n".format(silhouette_score, optimal_p))

occupancies = Counter(labels)
sorted_occupancies = sorted(occupancies.items(), key=lambda x: x[1], reverse=True)

with open(prefix + '_AH_cluster_occupancies.txt', 'w') as f:
    for label, count in sorted_occupancies:
        f.write(f'Cluster {label}: Occupancy {count}\n')

clusters, counts = zip(*sorted_occupancies)
plt.bar(clusters, counts)
plt.xlabel('Cluster')
plt.ylabel('Occupancy')
plt.title('Cluster Occupancies')
plt.savefig(prefix + '_AH_occupancies.png')
plt.clf()

cumulative_counts = np.cumsum(counts)
plt.plot(clusters, cumulative_counts)
plt.xlabel('Cluster')
plt.ylabel('Cumulative Occupancy')
plt.title('Cumulative Cluster Occupancies')
plt.savefig(prefix + '_AH_cumulative_occupancies.png')
plt.clf()

plt.figure()
plt.scatter(np.arange(frames), labels, marker='+')
plt.xlabel('Frame')
plt.ylabel('Cluster')
plt.title('iMWK-means with Explicit Rescaling and Kmeans')
plt.savefig(prefix + '_AH_cluster_vs_time.png')
plt.clf()

