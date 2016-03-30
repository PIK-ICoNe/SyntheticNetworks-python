__author__ = "Paul Schultz"
__date__ = "Mar 30, 2016"
__version__ = "v3.0"

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


import numpy as np
from scipy.sparse import dok_matrix
from igraph import Graph, plot, palettes, rescale
import os

from rpgm_algo import RpgAlgorithm

class NeoNet(RpgAlgorithm):
    def __init__(self):
        super(NeoNet, self).__init__()

    ###############################################################################
    # ##                       PRIVATE FUNCTIONS                               ## #
    ###############################################################################

    def _get_coords(self):
        if self.sampling == "uniform":
            return self._uniformunitsquare()
        else:
            print "ERROR: Not implemented yet."
            exit(1)

    def _uniformsquare(self, epsilon, point=None):
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

    m = np.array([2, 4, 8, 16, 32, 64])
    epsilon = np.array([1., 1. / 3., 1. / 9., 1. / 27., 1. / 81., 1. / 243.])
    points = []
    for layer in xrange(len(m)):
        points.append(np.empty([m[layer], 2]))

    for j in xrange(m[0]):
        points[0][j, :] = g._uniformsquare(epsilon[0])

    for layer in xrange(1, len(m)):
        for node in xrange(m[layer]):
            points[layer][node, :] = g._uniformsquare(epsilon[layer], points[layer - 1][node % m[layer - 1], :])

    import matplotlib.pyplot as plt
    fig = plt.figure()
    currentAxis = fig.add_subplot(111,aspect='equal')
    size = np.array([60, 50, 40, 30, 20, 10])
    col = np.array(["b", "r", "k", "g", "y", "gray"])
    for layer in xrange(len(m)):
        plt.scatter(points[layer][:, 0], points[layer][:, 1], s=size[layer], c=col[layer])
    for layer in xrange(len(m) - 1):
        for j in xrange(m[layer]):
            currentAxis.add_patch(plt.Rectangle(points[layer][j, :] - epsilon[layer + 1] * np.ones(2) / 2.,
                                                epsilon[layer + 1], epsilon[layer + 1], alpha=1, facecolor='none'))

    plt.show()



if __name__ == "__main__":
    main()






