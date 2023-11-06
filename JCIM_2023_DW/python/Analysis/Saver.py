import numpy



class Saver:
    def __init__(self, out_name):
        self.out_name = out_name

    def save(self):
        raise NotImplementedError

class Array(Saver):
    def __init__(self, array, out_name):
        self.array = array
        super().__init__(out_name)
    def save(self):
        numpy.savetxt(self.out_name, self.array)

