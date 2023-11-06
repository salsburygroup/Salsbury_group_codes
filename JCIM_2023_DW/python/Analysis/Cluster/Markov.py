import pyemma


class Model:
    def __init__(self, labels):
        self.labels = labels
        self.model = None

    def estimate(self):
        self.model = pyemma.msm.estimate_markov_model(dtrajs=self.labels, lag=1)
        return self.model
