__author__ = "Paul Schultz"
__date__ = "Mar 30, 2016"
__version__ = "v3.1"

# This file is based on the network creation algorithm published in:
#
# A Random Growth Model for Power Grids and Other Spatially Embedded Infrastructure Networks
# Paul Schultz, Jobst Heitzig, and Juergen Kurths
# Eur. Phys. J. Special Topics on "Resilient power grids and extreme events" (2014)
# DOI: 10.1140/epjst/e2014-02279-6
#
# The single-node basin stability predictor has appeared here:
#
# Detours around basin stability in power networks
# Paul Schultz, Jobst Heitzig, and Juergen Kurths
# New J. Phys. 16, 125001 (2014).
# DOI: 10.1088/1367-2630/16/12/125001


"""
layer parameters:
    - N^{i,L} final number of nodes in connected component L belonging to layer i
    - N_0^{i,L} number of seed nodes in connected component L belonging to layer i
    - p^i, q^i, s^i and r^i model parameters adapted to each layer i
    - \sigma^i fraction of connector nodes to other layers (node type allocation)
    - \epsilon^i link length scale, sublayer i grows in boxes around connectors of width \epsilon^i
"""


import numpy as np
from scipy.sparse import dok_matrix
from igraph import Graph, plot, palettes, rescale
import os

from rpgm.rpgm_core import RPG

class NeoNet(RPG):
    def __init__(self):
        super(NeoNet, self).__init__()

    def place_toplayer(self):
        # TODO generate transmission grid
        raise NotImplementedError()

    def place_sublayer(self):
        # TODO: use this to place a network instance for each root node of a distribution grid
        raise NotImplementedError()

    ###############################################################################
    # ##                       PRIVATE FUNCTIONS                               ## #
    ###############################################################################

    def _get_coords(self):
        if self.sampling == "uniform":
            return self._uniformunitsquare()
        else:
            print "ERROR: Not implemented yet."
            exit(1)

    @staticmethod
    def _uniformsquare(epsilon=1., point=None):
        """
        return point drawn uniformly at random
        from unit square of width epsilon
        centered around a given point
        -> 2D coordinates
        """
        if point is None:
            return np.random.uniform(size=2) * epsilon
        else:
            return np.random.uniform(size=2) * epsilon - epsilon * np.ones(2) / 2. + np.array(point)



#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


def main():
    g = NeoNet()

    g.set_params(n0=10, n=20, p=0.1, q=0.2)
    g.initialise()
    g.grow()

    g.plot_net()

    print g


if __name__ == "__main__":
    main()






