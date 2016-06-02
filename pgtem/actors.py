__author__ = "Paul Schultz"

import numpy as np


class Creator(object):

    def __init__(self, type, **kwargs):
        assert type in ["sink", "source"]
        self.type = type

    def __str__(self):
        return "Creator class."

    @classmethod
    def source(cls, **kwargs):
        return cls(type="source", **kwargs)

    @classmethod
    def sink(cls, **kwargs):
        return cls(type="sink", **kwargs)

    def new_node(self):
        return np.random.random(2)

class Coordinator(object):

    def __init__(self, **kwargs):
        pass

    def __str__(self):
        return "Coordinator class."

    def proposition(self, new_loc):
        assert np.shape(new_loc) == (2,)
        self.new_loc = np.aallowance = True
        possible_extensions = [
            [(0, 1), (1, 2)],
            [(0, 1), (1, 3)]
        ]
        return possible_extensions

class Provider(object):

    def __init__(self, type="HV", **kwargs):
        assert type in ["LV", "MV", "HV"]
        self.type = type

    def __str__(self):
        return "Provider class."

    @classmethod
    def LV(cls, **kwargs):
        return cls(type="LV", **kwargs)

    @classmethod
    def MV(cls, **kwargs):
        return cls(type="MV", **kwargs)

    @classmethod
    def HV(cls, **kwargs):
        return cls(type="HV", **kwargs)

    def get_costs(self, grid_extensions):
        return [1e6 for edges in grid_extensions]



class Regulator(object):

    def __init__(self, type="source", **kwargs):
        pass

    def __str__(self):
        return "Regulator class."

    def evaluate(self, grid_extensions):
        return [True for edges in grid_extensions]

if __name__ == '__main__':
    raise RuntimeWarning("Nothing will happen...")