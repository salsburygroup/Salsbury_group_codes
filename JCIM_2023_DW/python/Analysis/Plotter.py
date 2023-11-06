import matplotlib
matplotlib.use('Agg') # For use on DEAC cluster
import matplotlib.pyplot
import matplotlib.colors
import seaborn
import numpy


class Plotter:
    def __init__(self, y, out_name, x_label=' ', y_label=' ', title=' '):
        self.y = y
        self.out_name = out_name
        self.x_label = x_label
        self.y_label = y_label
        self.title = title

    def plot(self):
        raise NotImplementedError


class Y(Plotter):
    def plot(self):
        matplotlib.pyplot.figure()
        matplotlib.pyplot.plot(self.y)
        matplotlib.pyplot.xlabel(self.x_label)
        matplotlib.pyplot.ylabel(self.y_label)
        matplotlib.pyplot.title(self.title)
        matplotlib.pyplot.savefig(self.out_name)
        matplotlib.pyplot.close()


class XY(Plotter):
    def __init__(self, x, y, out_name, x_label=' ', y_label=' ', title=' '):
        self.x = x
        super().__init__(y, out_name, x_label, y_label, title)

    def plot(self):
        matplotlib.pyplot.figure()
        matplotlib.pyplot.plot(self.x, self.y)
        matplotlib.pyplot.xlabel(self.x_label)
        matplotlib.pyplot.ylabel(self.y_label)
        matplotlib.pyplot.title(self.title)
        matplotlib.pyplot.savefig(self.out_name)
        matplotlib.pyplot.close()


class SimplePColor(Plotter):
    def plot(self):
        matplotlib.pyplot.pcolor(self.y, cmap='jet')
        matplotlib.pyplot.colorbar()
        matplotlib.pyplot.ylim(0, len(self.y))
        matplotlib.pyplot.xlim(0, len(self.y))
        matplotlib.pyplot.xlabel(self.x_label)
        matplotlib.pyplot.ylabel(self.y_label)
        matplotlib.pyplot.title(self.title)
        matplotlib.pyplot.savefig(self.out_name)
        matplotlib.pyplot.close()

class UnityPColor(Plotter):
    def plot(self):
        matplotlib.pyplot.pcolor(self.y, cmap='jet', vmin=-1, vmax=1)
        matplotlib.pyplot.colorbar()
        matplotlib.pyplot.ylim(0, len(self.y))
        matplotlib.pyplot.xlim(0, len(self.y))
        matplotlib.pyplot.xlabel(self.x_label)
        matplotlib.pyplot.ylabel(self.y_label)
        matplotlib.pyplot.title(self.title)
        matplotlib.pyplot.savefig(self.out_name + '.png', pad_inches=0.03, bbox_inches='tight', dpi=400)
        matplotlib.pyplot.savefig(self.out_name + '.tiff', pad_inches=0.03, bbox_inches='tight', dpi=600)
        matplotlib.pyplot.savefig(self.out_name + '.pdf', bbox_inches='tight')
        matplotlib.pyplot.close()

class MeshPColor(Plotter):
    def __init__(self, y, x_edges, y_edges, out_name, x_label=' ', y_label=' ', title=' '):
        self.x_edges = x_edges
        self.y_edges = y_edges
        super().__init__(y, out_name, x_label, y_label, title)

    def plot(self):
        masked_y = numpy.ma.masked_where(numpy.isnan(self.y), self.y)
        matplotlib.pyplot.figure()
        matplotlib.pyplot.pcolormesh(self.x_edges, self.y_edges, masked_y.T)
        matplotlib.pyplot.colorbar()
        matplotlib.pyplot.xlabel(self.x_label)
        matplotlib.pyplot.ylabel(self.y_label)
        matplotlib.pyplot.title(self.title)
        matplotlib.pyplot.savefig(self.out_name)
        matplotlib.pyplot.close()


class MeshContour(Plotter):
    def __init__(self, y, x_edges, y_edges, out_name, x_label=' ', y_label=' ', title=' '):
        self.x_edges = x_edges
        self.y_edges = y_edges
        super().__init__(y, out_name, x_label, y_label, title)

    def plot(self):
        masked_y = numpy.ma.masked_where(numpy.isnan(self.y), self.y)
        extent = (self.y_edges[0], self.y_edges[-1], self.x_edges[0], self.x_edges[-1])
        matplotlib.pyplot.figure()
        matplotlib.pyplot.contourf(masked_y, extent=extent, cmap=matplotlib.pyplot.cm.get_cmap("nipy_spectral"))
        matplotlib.pyplot.colorbar().ax.set_ylabel('Free energy (kT)')
        matplotlib.pyplot.xlabel(self.x_label)
        matplotlib.pyplot.ylabel(self.y_label)
        matplotlib.pyplot.title(self.title)
        matplotlib.pyplot.savefig(self.out_name)
        matplotlib.pyplot.close()
