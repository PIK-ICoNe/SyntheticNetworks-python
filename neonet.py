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

from rpgm_core import RPG




#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


def main():

    g = RPG()
    assert(isinstance(g, RPG))

    #g.debug = True
    g.set_params(n=100, n0=1, r=1./3.)
    g.initialise()
    g.grow()

    print g

    g.stats


if __name__ == "__main__":
    main()






