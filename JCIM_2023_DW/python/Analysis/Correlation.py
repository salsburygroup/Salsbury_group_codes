import numpy
import mdtraj
import hdbscan
import os
import pandas
import signal
import subprocess


class Correlation:
    def __init__(self, trajectory):
        assert isinstance(trajectory, mdtraj.Trajectory)
        self.trajectory = trajectory
        self.correlation_matrix = []

    def calculate(self):
        raise NotImplementedError


class Pearson(Correlation):
    def calculate(self):
        average = numpy.average(self.trajectory.xyz, axis=0)
        fluctuations = self.trajectory.xyz - average[numpy.newaxis, :]
        del average
        dots = numpy.zeros((self.trajectory.n_atoms, self.trajectory.n_atoms))
        for i in range(self.trajectory.n_frames):
            dot = numpy.dot(fluctuations[i, :, :], numpy.transpose(fluctuations[i, :, :]))
            dots = dots + dot
        del fluctuations
        dots = numpy.divide(dots, self.trajectory.n_frames)
        diagonal = numpy.diag(dots)
        normalization_matrix = numpy.outer(diagonal, diagonal)
        normalization_matrix = numpy.sqrt(normalization_matrix)
        self.correlation_matrix = numpy.divide(dots, normalization_matrix)
        return self.correlation_matrix


class TimeLagged(Correlation):
    def __init__(self, trajectory, covariance_tau):
        assert isinstance(covariance_tau, int)
        self.covariance_tau=covariance_tau
        self.normalization_matrix = []
        super().__init__(trajectory)
    def calculate(self):
        average = numpy.average(self.trajectory.xyz, axis=0)
        fluctuations = self.trajectory.xyz - average[numpy.newaxis, :]
        del average
        dots = numpy.zeros((self.trajectory.n_atoms, self.trajectory.n_atoms))
        for i in range(self.trajectory.n_frames - self.covariance_tau):
            dot = numpy.dot(fluctuations[i, :, :], numpy.transpose(fluctuations[i + self.covariance_tau, :, :]))
            dots = dots + dot
        del fluctuations
        dots = numpy.divide(dots, self.trajectory.n_frames)
        diagonal = numpy.diag(dots)
        self.normalization_matrix = numpy.outer(diagonal, diagonal)
        self.normalization_matrix = numpy.sqrt(numpy.absolute(self.normalization_matrix))
        self.correlation_matrix = numpy.divide(dots, self.normalization_matrix)
        return self.correlation_matrix


class MutualInformation(Correlation):
    def calculate(self):
        raise NotImplementedError


class Propagator(Correlation):
    def __init__(self, trajectory, tau):
        self.tau = tau
        super().__init__(trajectory)

    def calculate(self):
        delta_sum = numpy.zeros([self.trajectory.topology.n_atoms, 3], dtype=float)
        dot_sum = numpy.zeros([self.trajectory.topology.n_atoms, self.trajectory.topology.n_atoms], dtype=float)
        for frame in numpy.arange(self.trajectory.n_frames - self.tau):
            delta_temp = self.trajectory.xyz[frame] - self.trajectory.xyz[frame + self.tau]
            delta_sum = delta_sum + delta_temp
            dot_sum = dot_sum + numpy.inner(delta_temp, delta_temp) # same as dot o f v and v' where v is <n_atoms, 3>
        average_delta = delta_sum / (self.trajectory.n_frames - self.tau)
        average_dot = dot_sum / (self.trajectory.n_frames - self.tau)
        dot_average_delta = numpy.inner(average_delta, average_delta)
        #_, normalization_matrix = Pearson(self.trajectory).calculate()
        diagonal = numpy.diag(average_dot)
        normalization_matrix = numpy.outer(diagonal, diagonal)
        normalization_matrix = numpy.sqrt(normalization_matrix)
        normalized_average_dot = numpy.divide(average_dot, normalization_matrix)
        return normalized_average_dot, average_delta, dot_average_delta,


class Clustering:
    @staticmethod
    def cluster(correlation_matrix, input_type='correlation', minimum_membership=None):
        number_residues = len(correlation_matrix)
        three_percent = int(numpy.ceil(number_residues * 0.03))
        if minimum_membership:
            min_cluster_size = minimum_membership
        elif three_percent >= 2:
            min_cluster_size = three_percent
        else:
            min_cluster_size = 2

        if input_type == 'similarity':
            clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, metric='precomputed')
            distance = 1 - numpy.abs(correlation_matrix)
            labels = clusterer.fit_predict(distance)
        elif input_type == 'dots':
            clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, metric='precomputed')
            dots = correlation_matrix.dot(correlation_matrix.T)
            distance = 1 - dots / numpy.max(numpy.abs(dots))
            labels = clusterer.fit_predict(distance)
        else:
            clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size)
            labels = clusterer.fit_predict(correlation_matrix)

        return labels

    @staticmethod
    def visualize(labels, pdb_file, out_name):
        num_residues = len(labels)
        df = pandas.DataFrame(columns=['residue', 'cluster'])
        df['residue'] = numpy.arange(num_residues)
        df['cluster'] = labels
        clusters = numpy.unique(df.loc[df['cluster'] != -1].cluster.values)
        noise_string = ' '.join(
            ['%d' % num for num in df.loc[df['cluster'] == -1].residue.values]
        )

        with open(out_name + '_image_script.vmd', 'w+') as vmd_file:
            vmd_file.write(
                'mol new ' + pdb_file + '\n' + 'mol delrep 0 top\n')
            if len(df.loc[df['cluster'] == -1]) > 0:
                vmd_file.write(
                    'mol representation NewCartoon\n'
                    + 'mol selection {residue ' + noise_string + '}\n' + 'mol material Ghost\n' + 'mol addrep top\n')
            vmd_file.write(
                'mol representation NewCartoon\n' + 'mol material AOChalky\n' + 'display ambientocclusion on\n')
            for cluster in clusters:
                cluster_string = ' '.join(
                    ['%d' % num for num in df.loc[df['cluster'] == cluster].residue.values]
                )
                vmd_file.write(
                    'mol color ColorID ' + str(cluster) + '\n' + 'mol selection {residue ' + cluster_string + '}\n'
                    + 'mol addrep top\n'
                )
            vmd_file.write(
                'display resize 1920 1080\n' + 'display resetview\n' + 'render TachyonInternal ' + out_name + '\nexit')

        # Posix-compliant way to exit shell
        signal.signal(signal.SIGTTOU, signal.SIG_IGN)
        # Now, let's make some pretty pictures
        vmd_render_cmd = (
            'vmd '
            + ' -dispdev text -e '
            + out_name + '_image_script.vmd'
        )
        subprocess.call([os.getenv('SHELL'), '-i', '-c', vmd_render_cmd])
        os.tcsetpgrp(0, os.getpgrp())
