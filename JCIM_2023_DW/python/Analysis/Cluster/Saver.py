import numpy
import pandas
import mdtraj
import os
from .. import AtomSelection, Featurizer
import scipy.spatial.distance
import subprocess
import PIL.Image
import glob
import copy


class Saver:
    def __init__(self, out_name):
        self.out_name = out_name

    def save(self):
        raise NotImplementedError


class TimeSeries(Saver):
    def __init__(self, out_name, labels):
        self.labels = labels
        super().__init__(out_name)

    def save(self):
        numpy.savetxt(self.out_name, self.labels)


class Score(Saver):
    def __init__(self, out_name, score):
        self.score = score
        super().__init__(out_name)

    def save(self):
        with open(self.out_name, 'w') as f:
            f.write(str(self.score))


class Scores(Saver):
    def __init__(self, out_name, scores):
        self.scores = scores
        super().__init__(out_name)

    def save(self):
        numpy.savetxt(self.out_name, self.scores)


class ClusterFrames(Saver):
    def __init__(self, out_name, labels):
        self.labels = labels
        super().__init__(out_name)
        num_frames = len(self.labels)
        self.clusters = int(max(self.labels)) + 1
        self.labeled_traj = pandas.DataFrame(columns=['frame', 'cluster'])
        self.labeled_traj['frame'] = numpy.arange(num_frames)
        self.labeled_traj['cluster'] = self.labels

    def save(self):
        with open(self.out_name, 'w') as f:
            for i in range(0, self.clusters):
                cluster_string = ' '.join(
                    ['%d' % num for num in self.labeled_traj.loc[self.labeled_traj['cluster'] == i].frame.values]
                )
                f.write(cluster_string + '\n')


class PDB(ClusterFrames):
    #def __init__(self, out_name, labels, trajectory, atom_selection='all'):
    def __init__(self, out_name, labels, trajectory):
        assert isinstance(trajectory, mdtraj.Trajectory)
        self.trajectory = trajectory
        #self.atom_selection = atom_selection
        super().__init__(out_name, labels)

    def save(self):
        #trajectory_temp_slice = AtomSelection.Slice(trajectory=copy.copy(self.trajectory), atom_selection=self.atom_selection).select()
        #trajectory_2d = Featurizer.XYZ(trajectory=trajectory_temp_slice).extract()
        trajectory_2d = Featurizer.XYZ(trajectory=self.trajectory).extract()
        for i in range(0, int(max(self.labels)) + 1):
            directory = self.out_name + '/cluster' + str(i)
            os.makedirs(directory)
            cluster_string = self.labeled_traj.loc[self.labeled_traj['cluster'] == i].frame.values
            cluster_coords = trajectory_2d[cluster_string]
            mean = cluster_coords.mean(axis=0)
            distance = [scipy.spatial.distance.euclidean(row, mean) for row in cluster_coords]
            rep = cluster_string[numpy.argmin(distance)]
            self.trajectory[rep].save_pdb(directory + '/rep.pdb')
            #self.trajectory[cluster_string].save_pdb(directory + '/all.pdb')
            self.trajectory[cluster_string].save_dcd(directory + '/all.dcd')


class Shadows(Saver):
    def __init__(self, out_name, middle, shadow, rep='NewCartoon'):
        assert isinstance(middle, str)
        assert isinstance(shadow, str)
        assert isinstance(rep, str)
        self.middle = middle
        self.shadow = shadow
        self.rep = rep
        super().__init__(out_name)

    def save(self):
        directory = os.path.dirname(__file__)
        print(directory)
        shadow_helper = os.path.join(directory, 'resources', 'generate_shadow.vmd')
        middle_helper = os.path.join(directory, 'resources', 'generate_middle.vmd')

        vmd_render_shadow_cmd = ('vmd ' + self.middle + ' ' + self.shadow + ' -dispdev text -e ' + shadow_helper + ' -args ' + ' -rep ' + self.rep)
        #vmd_render_shadow_cmd = ('vmd ' + self.shadow + ' -dispdev text -e ' + shadow_helper + ' -args ' + ' -rep ' + self.rep)
        vmd_render_middle_command = ('vmd ' + self.middle + ' -dispdev text -e ' + middle_helper + ' -args ' + ' -rep ' + self.rep + ' -outfile ' + self.out_name + '/middle.tga')
        print(self.out_name)
        print(vmd_render_middle_command)
        print(os.getenv('SHELL'), '-i', '-c', 'cd ' + self.out_name + ' && ls && ' + vmd_render_shadow_cmd + '; ' + vmd_render_middle_command + '; exit')

        subprocess.call(
            [os.getenv('SHELL'), '-i', '-c', 'cd ' + self.out_name + ' && ls && ' + vmd_render_middle_command + '; ' + vmd_render_shadow_cmd + '; exit'],
            cwd=self.out_name
        )
        shadow_pattern = self.out_name + "/shadow.*.tga"
        if os.path.isfile(self.out_name + '/shadow.png'):
            os.remove(self.out_name + '/shadow.png')

        for file in glob.glob(shadow_pattern):
            shadow_img = PIL.Image.open(file)
            shadow_img = shadow_img.convert("RGBA")
            shadow_data = shadow_img.getdata()

            new_shadow_data = []
            for item in shadow_data:
                if item[0] == 255 and item[1] == 255 and item[2] == 255:
                    new_shadow_data.append((255, 255, 255, 0))
                else:
                    new_shadow_data.append((item[0], item[1], item[2], 26))

            shadow_img.putdata(new_shadow_data)

            if os.path.isfile(self.out_name + '/shadow.png'):
                shadow_old = PIL.Image.open(self.out_name + '/shadow.png')
                shadow_new = PIL.Image.alpha_composite(shadow_img, shadow_old)
                shadow_new.save(self.out_name + '/shadow.png')
            else:
                shadow_img.save(self.out_name + '/shadow.png', "PNG")

        for file in glob.glob(shadow_pattern):
            os.remove(file)

        # Let's get rid of the white pixels and convert the TGAs to PNGs
        middle_image = PIL.Image.open(self.out_name + '/middle.tga')
        middle_image = middle_image.convert("RGBA")
        middle_data = middle_image.getdata()

        new_middle_data = []
        for item in middle_data:
            if item[0] == 89 and item[1] == 89 and item[2] == 89:
                new_middle_data.append((89, 89, 89, 0))
            else:
                new_middle_data.append(item)

        middle_image.putdata(new_middle_data)
        middle_image.save(self.out_name + '/middle.png', "PNG")

        # Now, let's layer them together
        layered_img = PIL.Image.alpha_composite(shadow_new, middle_image)
        layered_img.save(self.out_name + '/layered.png', "PNG")

        blended_img = PIL.Image.blend(middle_image, shadow_new, 0.5)
        blended_img.save(self.out_name + '/blended.png', "PNG")
