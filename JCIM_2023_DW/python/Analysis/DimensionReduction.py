import numpy
import sklearn.decomposition


class Reducer:
    def __init__(self, coordinates):
        self.coordinates = coordinates

    def reduce(self):
        raise NotImplementedError


class Projector:
    def __init__(self, coordinates, basis_set):
        self.coordinates = coordinates
        self.basis_set = basis_set
        self.projection = []

    def project(self):
        self.projection = numpy.zeros_like(self.coordinates)
        for point_idx in range(self.coordinates.shape[1]):
            for vec_idx in range(self.basis_set.shape[0]):
                self.projection[vec_idx, point_idx] = \
                    numpy.dot(self.basis_set[:, vec_idx], self.coordinates[:, point_idx])
        return self.projection


class PCA(Reducer):
    def __init__(self, coordinates):
        self.projection = []
        super().__init__(coordinates)

    def reduce(self):
        model = sklearn.decomposition.PCA().fit(self.coordinates)
        self.projection = sklearn.decomposition.PCA().fit_transform(self.coordinates)
        return self.projection, model.components_, model.explained_variance_ratio_


class TICA(Reducer):
    def __init__(self, coordinates):
        self.projection = []
        self.model = []
        super().__init__(coordinates)

    def reduce(self):
        raise NotImplementedError
