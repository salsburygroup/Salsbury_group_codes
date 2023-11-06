import numpy
import mdtraj
from .AtomSelection import Slice


class Distance:
    def __init__(self, trajectory, atom_selection):
        assert isinstance(trajectory, mdtraj.Trajectory)
        self.trajectory = trajectory
        self.sel = atom_selection

    def calculate(self):
        raise NotImplementedError

class fluctuation(Distance):
    def calculate(self):
        sub_trajectory = Slice(trajectory=self.trajectory, atom_selection=self.sel).select()
        assert isinstance(sub_trajectory, mdtraj.Trajectory)
        reference_positions = sub_trajectory.xyz.mean(0)
        difference = sub_trajectory.xyz - reference_positions
        #sum_squares = numpy.sum(numpy.sum(numpy.square(difference), axis=2), axis=0)
        #rmsf = (sum_squares/sub_trajectory.n_frames)**0.5
        return difference

class RMSF(Distance):
    def calculate(self):
        sub_trajectory = Slice(trajectory=self.trajectory, atom_selection=self.sel).select()
        assert isinstance(sub_trajectory, mdtraj.Trajectory)
        reference_positions = sub_trajectory.xyz.mean(0)
        difference = sub_trajectory.xyz - reference_positions
        sum_squares = numpy.sum(numpy.sum(numpy.square(difference), axis=2), axis=0)
        # Code at https://github.com/schilli/Tools/blob/master/rmsf.py uses 3*n_frames. Seems wrong...
        rmsf = (sum_squares/sub_trajectory.n_frames)**0.5
        return rmsf

class RMSF_separate(Distance):
    def calculate(self):
        sub_trajectory = Slice(trajectory=self.trajectory, atom_selection=self.sel).select()
        assert isinstance(sub_trajectory, mdtraj.Trajectory)
        reference_positions = sub_trajectory.xyz.mean(0)
        difference = sub_trajectory.xyz - reference_positions
        difference1=[]
        difference2=[]
        difference3=[]
        difference4=[]
        difference5=[]
        each_length=int(len(difference)/5)
        for i in range(each_length):
            difference1.append(difference[i])
            difference2.append(difference[i+each_length])
            difference3.append(difference[i+2*each_length])
            difference4.append(difference[i+3*each_length])
            difference5.append(difference[i+4*each_length])
#        print('difference')
#        print(difference[0])
#        print('difference1+difference2+difference3+difference4+difference5')
        #print(difference1+difference2+difference3+difference4+difference5)
#        print(difference1[0])
#        print('len:')
#        print(len(difference))
#        print(len(difference1+difference2+difference3+difference4+difference5))

        sum_squares = numpy.sum(numpy.sum(numpy.square(difference), axis=2), axis=0)
        sum_squares1 = numpy.sum(numpy.sum(numpy.square(difference1), axis=2), axis=0)
        sum_squares2 = numpy.sum(numpy.sum(numpy.square(difference2), axis=2), axis=0)
        sum_squares3 = numpy.sum(numpy.sum(numpy.square(difference3), axis=2), axis=0)
        sum_squares4 = numpy.sum(numpy.sum(numpy.square(difference4), axis=2), axis=0)
        sum_squares5 = numpy.sum(numpy.sum(numpy.square(difference5), axis=2), axis=0)
#        print('sum_squares:')
#        print(sum_squares)
#        print(sum_squares1)
#        print(sum_squares2)
#        print(sum_squares3)
#        print(sum_squares4)
#        print(sum_squares5)
        # Code at https://github.com/schilli/Tools/blob/master/rmsf.py uses 3*n_frames. Seems wrong...
        rmsf = (sum_squares/sub_trajectory.n_frames)**0.5
        rmsf1 = (5*sum_squares1/sub_trajectory.n_frames)**0.5
        rmsf2 = (5*sum_squares2/sub_trajectory.n_frames)**0.5
        rmsf3 = (5*sum_squares3/sub_trajectory.n_frames)**0.5
        rmsf4 = (5*sum_squares4/sub_trajectory.n_frames)**0.5
        rmsf5 = (5*sum_squares5/sub_trajectory.n_frames)**0.5
#        print('rmsf:')
#        print(rmsf)
#        print('average_rmsf:')
#        print(((rmsf1**2+rmsf2**2+rmsf3**2+rmsf4**2+rmsf5**2)/5)**0.5)
#        print('num_of_frames:')
#        print(sub_trajectory.n_frames)
        return [rmsf1, rmsf2, rmsf3, rmsf4, rmsf5]

class RMSD(Distance):
    #def __init__(self, trajectory, atom_selection, reference_frame=0):
    def __init__(self, trajectory, atom_selection, reference_frame):
        self.reference_frame = reference_frame
        super().__init__(trajectory,atom_selection)
    def calculate(self):
        sub_trajectory = Slice(trajectory=self.trajectory, atom_selection=self.sel).select()
        assert isinstance(sub_trajectory, mdtraj.Trajectory)
        rmsd = mdtraj.rmsd(self.trajectory, self.trajectory, frame=self.reference_frame)
        return rmsd
