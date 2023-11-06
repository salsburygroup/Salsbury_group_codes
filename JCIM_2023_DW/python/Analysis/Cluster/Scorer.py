import numpy
import sklearn.metrics
import warnings


class Scorer:
    def __init__(self, labels, data):
        self.labels = labels
        self.data = data
        self.score = []

    def evaluate(self):
        raise NotImplementedError


class CalinskiHarabasz(Scorer):
    # based on a partial implementation in scikit-learn pull requests
    def __init__(self, labels, data, centers):
        self.centers = centers
        super().__init__(labels, data)

    def evaluate(self):
        mean = numpy.mean(self.data, axis=0)
        between_cluster_variance = numpy.sum([(center - mean)**2 for center in self.centers])
        within_cluster_variance = numpy.sum([(x - self.centers[self.labels[i]])**2 for i, x in enumerate(self.data)])
        num_clusters = len(self.centers)
        samples = len(self.data)
        self.score = (between_cluster_variance / (num_clusters-1))/(within_cluster_variance / (samples - num_clusters))
        return self.score


class DaviesBouldin(Scorer):
    def __init__(self, labels, data, centers):
        self.centers = centers
        super().__init__(labels, data)

    def evaluate(self):
        # based on https://github.com/jqmviegas/jqm_cvi/blob/master/base.py
        db = 0
        num_clusters = len(self.centers)
        within_cluster_scatter = numpy.empty(num_clusters, dtype=float)
        clusters = numpy.unique(self.labels)
        for i in range(num_clusters):
            cluster = clusters[i]
            points_in_cluster = self.data[numpy.where(self.labels == cluster)]
            norm = 1 / points_in_cluster.shape[0]
            within_cluster_scatter[i] = norm * numpy.sum(
                [numpy.linalg.norm(points_in_cluster[j] - self.centers[i])
                 for j in range(points_in_cluster.shape[0])]
            )
        cluster_separation = numpy.empty([num_clusters, num_clusters], dtype=float)
        for i in range(num_clusters):
            for j in range(num_clusters):
                cluster_separation[i, j] = numpy.linalg.norm(self.centers[i] - self.centers[j])
        for i in range(num_clusters):
            d = numpy.empty(num_clusters-1, dtype=float)
            for j in range(0, i):
                d[j] = (within_cluster_scatter[i] + within_cluster_scatter[j]) / cluster_separation[i, j]
            for j in range(i+1, num_clusters):
                d[j-1] = (within_cluster_scatter[i] + within_cluster_scatter[j]) / cluster_separation[i, j]
            db += numpy.max(d)
        self.score = db/num_clusters
        return self.score


class Silhouette(Scorer):
    def evaluate(self):
        if self.data.shape[0] > 10000:
            sample_size = 10000
            warnings.warn("More than 10,000 samples. "
                          "Silhouette will score on random sampling of 10,000."
                          "To override, score directly with scikit-learn")
        else:
            sample_size = None
        self.score = sklearn.metrics.silhouette_score(self.data, self.labels, sample_size=sample_size)
        return self.score
