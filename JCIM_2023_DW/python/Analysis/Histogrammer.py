import numpy


class TwoDimensions:
    def __init__(self, vector1, vector2):
        self.vector1 = vector1
        self.vector2 = vector2

    def bin(self):
        raise NotImplementedError


class FreedmanDiaconis(TwoDimensions):
    def bin(self):
        temp = numpy.column_stack((self.vector1, self.vector2))
        q75 = numpy.percentile(temp, 75)
        q25 = numpy.percentile(temp, 25)
        irq = q75 - q25
        bin_width = 2 * irq / (len(self.vector1) ** (1. / 3))
        bins = numpy.ceil((numpy.max(temp[:, 0]) - numpy.min(temp[:, 0])) / bin_width)
        histogram, x_edges, y_edges = numpy.histogram2d(temp[:, 0], temp[:, 1], bins=bins)
        return histogram, x_edges, y_edges


class Rice(TwoDimensions):
    def bin(self):
        bins = numpy.floor(2 * len(self.vector1) **(1./3))
        histogram, x_edges, y_edges = numpy.histogram2d(self.vector1, self.vector2, bins=bins)
        return histogram, x_edges, y_edges
