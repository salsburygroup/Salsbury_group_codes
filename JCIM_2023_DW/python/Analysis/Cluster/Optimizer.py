import numpy


class Optimizer:
    def __init__(self, scores, parameter_list):
        self.scores = scores
        self.parameter_list = parameter_list

    def minimize(self):
        optimal_parameter = self.parameter_list[numpy.argmin(self.scores)]
        return optimal_parameter

    def maximize(self):
        optimal_parameter = self.parameter_list[numpy.argmax(self.scores)]
        return optimal_parameter


class Slope(Optimizer):
    def __init__(self, scores, parameter_list, exponent=1):
        self.exponent = exponent
        super().__init__(scores, parameter_list)

    def minimize(self):
        slope_indices = numpy.zeros(max(self.parameter_list) + 1)
        for s in range(self.scores.shape[0] - 1):
            slope_index = -(self.scores[s + 1] - self.scores[s]) \
                          * self.scores[s] ** self.exponent
            slope_indices[s] = slope_index
        optimal_parameter = self.parameter_list[numpy.argmin(slope_indices)]
        return optimal_parameter

    def maximize(self):
        raise NotImplementedError
