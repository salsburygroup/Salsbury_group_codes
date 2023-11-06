import hdbscan
from . import Unsupervised, Scorer, Optimizer
from sklearn.cluster import MiniBatchKMeans
from sklearn import mixture
from sklearn import cluster
import numpy as np
import copy

class Clusterer:
    """
    Base class for clustering algorithms.
    """
    def __init__(self, trajectory_2d):
        """
        Initialize the Clusterer.

        Parameters:
        - trajectory_2d : array-like, shape [n_samples, n_features]
        """
        self.trajectory_2d = trajectory_2d
        self.labels = []

    def fit(self):
        """
        Fit the model to the data. Must be implemented by subclasses.
        """
        raise NotImplementedError("Subclasses must implement fit method")

class HDBSCAN(Clusterer):
    def __init__(self, trajectory_2d, minimum_membership=5, min_samples=5, core_distance = None):
        self.core_distance = core_distance
        self.minimum_membership = minimum_membership
        self.min_samples = min_samples
        super().__init__(trajectory_2d)  # inherits constructor

    def fit(self):
        framework = hdbscan.HDBSCAN(min_cluster_size=self.minimum_membership, min_samples=self.min_samples, core_distance=self.core_distance)
        self.labels = framework.fit_predict(self.trajectory_2d)
        return self.labels

class KMeans(Clusterer):
    def __init__(self, trajectory_2d):
        self.centers = []
        super().__init__(trajectory_2d)

    def fit(self):
        cl = Unsupervised.clustering.Clustering()
        labels, _, _, _, _, = cl.ik_means(data=self.trajectory_2d)
        max_clusters = max(labels) + 1
        k_to_try = range(2, max_clusters + 2)  # Have to go +2 for slope
        scores = np.zeros(max_clusters)
        for k in k_to_try:
            clusterer = MiniBatchKMeans(
                n_clusters=k, n_init=100)
            labels = clusterer.fit_predict(self.trajectory_2d)
            scores[k - 2] = Scorer.Silhouette(data=self.trajectory_2d, labels=labels).evaluate()
        optimizer = Optimizer.Slope(scores, k_to_try)
        num_clusters = optimizer.minimize()
        clusterer = MiniBatchKMeans(n_clusters=num_clusters)
        self.labels = clusterer.fit_predict(self.trajectory_2d)
        self.centers = clusterer.cluster_centers_
        return self.labels, self.centers


class GMM(Clusterer):
    def __init__(self, trajectory_2d, max_clusters=100):
        self.max_clusters = max_clusters
        self.centers = []
        super().__init__(trajectory_2d)

    def fit(self):
        aic = np.zeros(self.max_clusters - 1)
        k_to_try = range(2, self.max_clusters + 1)
        for k in k_to_try:
            clusterer = mixture.GMM(
                n_components=k,
                covariance_type='tied',
                n_init=5
            )
            cluster = clusterer.fit(self.trajectory_2d)
            clusterer.predict(self.trajectory_2d)
            aic[k - 2] = cluster.aic(self.trajectory_2d)
        optimizer = Optimizer.Optimizer(aic, k_to_try)
        num_clusters = optimizer.minimize()
        clusterer = mixture.GMM(
            n_components=num_clusters,
            covariance_type='tied'
        )
        self.labels = clusterer.fit_predict(self.trajectory_2d)
        self.centers = clusterer.means_
        return self.labels, self.centers


class IMWKRescaled(Clusterer):
    def __init__(self, trajectory_2d, minkowski_weight=2):
        self.minkowski_weight = minkowski_weight
        self.centers = []
        self.data = []
        self.weights = []
        self.optimal_k = []
        super().__init__(trajectory_2d)

    def fit(self):
        cl = Unsupervised.clustering.Clustering()
        # Get maximum number of clusters
        [labels, _, _, _, _] = cl.ik_means(self.trajectory_2d, p=self.minkowski_weight)
        max_clusters = max(labels) + 1
        silhouette_averages = np.zeros(max_clusters - 1)
        # Tru all less than or equal
        k_to_try = np.arange(2, max_clusters + 1)
        for k in k_to_try:
            cl = Unsupervised.clustering.Clustering()
            data = copy.copy(self.trajectory_2d)
            [labels, centroids, weights, _, _] = cl.imwk_means(data, p=self.minkowski_weight, k=k)
        # Rescale the data
            populations = np.bincount(labels)
            for k1 in np.arange(0, max(labels)+1):
                data[labels == k1] = np.multiply(
                    self.trajectory_2d[labels == k1], np.tile(weights[k1], (populations[k1], 1))
                )
                centroids[k1] = np.multiply(centroids[k1], weights[k1])
        # Apply Euclidean KMeans
            kmeans_clusterer = MiniBatchKMeans(n_clusters=k, n_init=5)
            kmeans_clusters = kmeans_clusterer.fit(data)
            labels = kmeans_clusters.labels_
            silhouette_averages[k - 2] = Scorer.Silhouette(labels=labels, data=data).evaluate()

        self.optimal_k = Optimizer.Optimizer(scores=silhouette_averages, parameter_list=k_to_try).maximize()
        # Do optimal clustering
        self.data = copy.copy(self.trajectory_2d)
        [self.labels, self.centers, self.weights, _, _] = cl.imwk_means(
            self.data, p=self.minkowski_weight, k=self.optimal_k)
        # Rescale the data
        for k1 in np.arange(0, max(self.labels)+1):
            populations = np.bincount(self.labels)
            self.data[self.labels == k1] = np.multiply(self.data[self.labels == k1],
                                                np.tile(self.weights[k1], (populations[k1], 1)))
            self.centers[k1] = np.multiply(self.centers[k1], self.weights[k1])
        return self.labels, self.centers, self.data, self.weights, self.optimal_k


class AmorimHennig(IMWKRescaled):
    def fit(self):
        # Apply Euclidean KMeans
        _, _, self.data, _, self.optimal_k = IMWKRescaled(
            trajectory_2d=self.trajectory_2d, minkowski_weight=self.minkowski_weight).fit()
        kmeans_clusterer = MiniBatchKMeans(n_clusters=self.optimal_k, n_init=5)
        kmeans_clusters = kmeans_clusterer.fit(self.data)
        self.labels = kmeans_clusters.labels_
        self.centers = kmeans_clusterer.cluster_centers_
        return self.labels, self.centers


class VBGMM(Clusterer):
    def __init__(self, trajectory_2d, max_clusters=100):
        self.max_clusters = max_clusters
        self.centers = []
        super().__init__(trajectory_2d)

    def fit(self):
        aic = np.zeros(self.max_clusters - 1)
        k_to_try = range(2, self.max_clusters + 1)
        for k in k_to_try:
            clusterer = mixture.VBGMM(
                n_components=k,
                covariance_type='tied',
                )
            cluster = clusterer.fit(self.trajectory_2d)
            _ = clusterer.predict(self.trajectory_2d)
            aic[k - 2] = cluster.aic(self.trajectory_2d)
        optimizer = Optimizer.Optimizer(aic, k_to_try)
        num_clusters = optimizer.minimize()
        clusterer = mixture.VBGMM(
            n_components=num_clusters,
            covariance_type='tied'
        )
        self.labels = clusterer.fit_predict(self.trajectory_2d)
        self.centers = clusterer.means_
        return self.labels, self.centers


class MeanShift(Clusterer):
    def __init__(self, trajectory_2d, bandwidth=None, cluster_all=True):
        self.bandwidth = bandwidth
        self.cluster_all = cluster_all
        self.centers = []
        super().__init__(trajectory_2d)

    def fit(self):
        clusterer = cluster.MeanShift(bandwidth=self.bandwidth, n_jobs=-1, cluster_all=self.cluster_all)
        self.labels = clusterer.fit_predict(self.trajectory_2d)
        self.centers = clusterer.cluster_centers_
        return self.labels, self.centers


class AffinityPropagation(Clusterer):
    def __init__(self, trajectory_2d, damping=0.5, preference=None):
        self.damping = damping
        self.preference = preference
        self.centers = []
        super().__init__(trajectory_2d)

    def fit(self):
        clusterer = cluster.AffinityPropagation(damping=self.damping, preference=self.preference)
        self.labels = clusterer.fit_predict(self.trajectory_2d)
        self.centers = clusterer.cluster_centers_
        return self.labels, self.centers


class QualityThreshold:
    def __init__(self, distances, cutoff):
        assert isinstance(distances, np.ndarray)
        self.distances = distances
        self.cutoff = cutoff
        self.labels = []
        self.centers = []

    def fit(self):
        cutoff_mask = self.distances < self.cutoff
        cluster = 0
        self.labels = np.empty(self.distances.shape[0])
        self.labels.fill(np.NAN)

        while cutoff_mask.any():
            membership = cutoff_mask.sum(axis=1)
            center = np.argmax(membership)
            members = np.where(cutoff_mask[center, :] == True)
            if max(membership) == 1:
                self.labels[np.where(np.isnan(self.labels))] = -1
                break
            self.labels[members] = cluster
            self.centers.append(center)
            cutoff_mask[members, :] = False
            cutoff_mask[:, members] = False
            cluster = cluster + 1

        return self.labels, self.centers
