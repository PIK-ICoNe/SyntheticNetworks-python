#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = "Paul Schultz"
__date__ = "Jun 21, 2016"
__version__ = "v0.1"
__copyright__ = "Copyright (C) 2016 Paul Schultz, GNU GPL 3"

"""
Module contains class PowerGrid as well as the classes Creator, Coordinator, Provider and Regulator.
"""

#TODO: description

# Import NumPy for the array object and fast numerics.
import numpy as np
# Import Pandas for data structures.
import pandas as pd
# Import SciPy sparse DOK matrices to handle large adjacency matrices.
# The DOK format is optimal when new entries are sequentially added.
from scipy.sparse import dok_matrix


class Creator(object):
    """
    Creator class

    The class Creator contains methods for a "Creator"-agent
    (an instance) to suppose new sites for power generation and loads.

    This class can either be instantiated as a sink (load) or source (generator).

    **Examples:**

    >>> print Creator.sink()
    Creator class.

    """

    def __init__(self, type, **kwargs):
        """
        Initialise an instance of Creator.

        Parameters
        ----------
        type: string
            flag to determine how Creator is instantiated

        """
        assert type in ["sink", "source"]
        self.type = type

    def __str__(self):
        """
        Output of passing a Creator instance to the print command.

        Returns
        -------
        text: string

        **Examples:**

        >>> print Creator.sink()
        Creator class.

        """
        text = "Creator class."
        return text

    @classmethod
    def source(cls, **kwargs):
        return cls(type="source", **kwargs)

    @classmethod
    def sink(cls, **kwargs):
        return cls(type="sink", **kwargs)

    def new_node(self, BusID, time):
        """
        New node at a random location. Node type is determined by class instance.

        Parameters
        ----------
        BusID: int
            a unique node identifier
        time: int/float
            a timestamp attribute

        Returns
        -------
        subgraph: DataFrame object
            Table containing the following node attributes: BusID, lat, lon, time, type

        **Examples:**
        >>> c = Creator.source()
        >>> print c.new_node(BusID=0, time=0) # doctest:+ELLIPSIS
           BusID       lat       lon  time    type
        ...

        """
        
        loc = np.random.random(2)
        return pd.DataFrame(data={"BusID": int(BusID), "lat": loc[0], "lon": loc[1], "time": time, "type":self.type},
        index=range(1))

class Coordinator(object):
    """
        Coordinator class

        The class Coordinator contains methods for a "Coordinator"-agent
        (an instance) to do two things:
        a) propose possible grid extensions to connect new sites.
        b) negotiate between Grid Providers and the Regulator instance to find cheapest allowed extension.

        **Examples:**

        >>> print Coordinator()
        Coordinator class.

        """

    def __init__(self, **kwargs):
        """
        Initialise an instance of Coordinator.

        """
        pass

    def __str__(self):
        """
        Output of passing a Coordinator instance to the print command.

        Returns
        -------
        text: string

        **Examples:**

        >>> print Coordinator()
        Coordinator class.

        """
        text = "Coordinator class."
        return text

    def proposition(self, node, grid):
        assert isinstance(grid, PowerGrid)
        assert isinstance(node, pd.DataFrame)

        assert grid.number_of_nodes > 0

        new_loc = np.squeeze([node.lat.values, node.lon.values])
        closest_node = grid.get_closest_node(new_loc)
        possible_extensions = [[(np.squeeze(node.BusID.values) - 1, closest_node),],]
        return possible_extensions

    def negotiation(self, possible_extensions, providers, regulator):

        assert [isinstance(provider, Provider) for provider in providers]
        assert isinstance(regulator, Regulator)

        # check for empty proposition
        assert possible_extensions

        # regulator filters possible extensions
        allowed = regulator.evaluate(possible_extensions)

        if not True in allowed:
            # stop if nothing is allowed
            return False, None
        else:
            # get costs for allowed extensions
            for i in xrange(len(allowed)):
                if not allowed[i]:
                    possible_extensions.pop(i)

            possible_extensions = np.array(possible_extensions)

            costs = [provider.get_costs(possible_extensions) for provider in providers]

            min_costs = np.min(costs, axis=0)

            # rank extensions by costs
            return True, {i: v for i,v in enumerate(possible_extensions[np.argsort(min_costs)])}



class Provider(object):
    """
    Provider class

    The class Provider contains methods for a "Provider"-agent
    (an instance) to put price tags on the grid extensions proposed
    by the Coordinator instance.

    This class can be instantiated as an LV, MV or HV (incl. UHV) grid provider.

    The role of this agent might be extended in the future.

    **Examples:**

    >>> print Provider.HV()
    Provider class.

    """

    def __init__(self, type="HV", **kwargs):
        """
        Initialise an instance of Provider.

        Parameters
        ----------
        type: string
            flag to determine how Provider is instantiated

        """
        assert type in ["LV", "MV", "HV"]
        self.type = type

    def __str__(self):
        """
        Output of passing a Provider instance to the print command.

        Returns
        -------
        text: string

        **Examples:**

        >>> print Provider.HV()
        Provider class.

        """

        text = "Provider class."
        return text

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
    """
    Regulator class

    The class Regulator contains methods for a "Regulator"-agent
    (an instance) to evaluate wether grid extensions proposed by
    a Coordinator agent comply with grid codes (== are allowed).

    **Examples:**

    >>> print Regulator()
    Regulator class.

    """

    def __init__(self, type="source", **kwargs):
        """
        Initialise an instance of Regulator.

        """
        pass

    def __str__(self):
        """
        Output of passing a Regulator instance to the print command.

        Returns
        -------
        text: string

        **Examples:**

        >>> print Regulator()
        Regulator class.

        """

        text = "Regulator class."
        return text

    def evaluate(self, grid_extensions):
        return [True for edges in grid_extensions]



class PowerGrid(object):
    """
    PowerGrid class

    The class PowerGrid contains methods for a "PowerGrid"-agent
    (an instance) to suppose new sites for power generation and loads.

    **Examples:**

    >>> print PowerGrid(n_final=10)
    An instance of the PowerGrid class.
    Final number of nodes: 10

    """

    def __init__(self, n_final):
        """
        Initialise an instance of PowerGrid.

        Parameters
        ----------
        n_final: int
            determines final number of nodes

        """

        self.n_final = n_final

        # empty adjacency
        self.adjacency = dok_matrix(np.zeros([self.n_final, self.n_final]))

        # set internal counters for nodes and edges
        self.number_of_nodes = 0
        self.number_of_edges = 0

        # data containers for nodes and edges
        self.nodes = pd.DataFrame(columns=["BusID", "type", "lon", "lat", "time"])
        self.edges = pd.DataFrame(columns=["BranchID", "source", "target", "type", "time"])

        # internal settings (do not change)
        self.distance_measure = "euclidean"

    def __str__(self):
        """
        Output of passing a PowerGrid instance to the print command.

        Returns
        -------
        text: string

        **Examples:**

        >>> print PowerGrid(n_final=10)
        An instance of the PowerGrid class.
        Final number of nodes: 10

        """

        text = "An instance of the PowerGrid class.\n"
        text += "Final number of nodes: " + str(self.n_final)
        return text

    def add_node(self, new_node):
        # entry in node data table
        self.nodes = pd.merge(self.nodes, new_node, how="outer")
        # update counter
        self.number_of_nodes += 1
        #TODO: extend adjacency if number_of_nodes > len(adjacency)


    def add_edge(self, source, target, type, time):
        new_idx = self.number_of_edges + 1
        # entry in edge data table
        length = self._get_distance(int(source), int(target))
        if source < target:
            x, y = source, target
        else:
            x, y = target, source
        new_edge = pd.DataFrame(data={"BranchID": new_idx, "source": int(x), "target": int(y),
                                      "length":length, "type": type, "time": time}, index=range(1))
        self.edges = pd.merge(self.edges, new_edge, how="outer")
        # add edge in adjacency matrix
        self.adjacency[int(source), int(target)] = \
            self.adjacency[int(target), int(source)] = length
        # update counter
        self.number_of_edges += 1
        return new_idx

    def import_from_rpg(self, g=None, **kwargs):
        #TODO: recursive import ?!
        from rpgm.rpgm_algo import RpgAlgorithm

        if g is None:
            g = RpgAlgorithm()
        else:
            assert isinstance(g, RpgAlgorithm)

        g.set_params(**kwargs)
        g.initialise()
        g.grow()

        # update node data table
        new_nodes = pd.DataFrame(data={"BusID": 1 + np.arange(g.added_nodes),
                                       "lat": g.lat,
                                       "lon": g.lon,
                                       "time": np.repeat(-1, g.added_nodes),
                                       "type": np.repeat("rpg_node", g.added_nodes)},
                                 index=range(g.added_nodes))
        self.nodes = pd.merge(self.nodes, new_nodes, how="outer")
        # update node counter
        self.number_of_nodes += g.added_nodes

        for edge in g.edge:
            if edge[0] < edge[1]:
                x, y = edge[0], edge[1]
            else:
                x, y = edge[1], edge[0]
            length = g.distance[edge]
            # update edge data table
            new_edge = pd.DataFrame(data={"BranchID": 1 + self.number_of_edges,
                                          "source": int(x),
                                          "target": int(y),
                                          "type": "rpg_edge",
                                          "length": length,
                                          "time": -1},
                                    index=range(1))
            self.edges = pd.merge(self.edges, new_edge, how="outer")
            # add edge to adjacency matrix
            self.adjacency[edge[0], edge[1]] = self.adjacency[edge[1], edge[0]] = length
            # update edge counter
            self.number_of_edges += 1

    def convert_to_igraph(self):
        import igraph as ig
        pass

    def update(self, new_node=None, edgelist=None):
        # add new node and all edges in edgelist

        if not new_node is None:
            self.add_node(new_node)

        if not edgelist is None:
            for edge in edgelist:
                self.add_edge(edge[0], edge[1], "line", np.squeeze(new_node.time))

    def save(self, path, raw=False):
        pd.to_pickle(self.nodes, path + ".nodes")
        pd.to_pickle(self.edges, path + ".edges")
        np.save(path + ".adj", self.adjacency)


    def get_closest_node(self, source):
        # TODO: rework this
        # vertices need to be properly ordered for this to work, i.e. nodes in the connected component
        # should be labeled from 0 to connected-1
        min = np.inf
        target = -1
        for node in xrange(self.number_of_nodes):
            y = np.array([self.nodes.lat.values[node], self.nodes.lon.values[node]])
            d = self._euclidean(source, y)
            if d < min:
                min = d
                target = node
        return target

    def _get_distance(self, u, v):
        x = np.array([self.nodes.lat.values[int(u)], self.nodes.lon.values[int(u)]])
        y = np.array([self.nodes.lat.values[int(v)], self.nodes.lon.values[int(v)]])
        if self.distance_measure == "euclidean":
            return self._euclidean(x, y)
        else:
            print "ERROR: Not implemented yet."
            exit(1)

    def _euclidean(self, x, y):
        """
        return euclidean distance between x and y
        """
        return np.sqrt(sum((x - y) ** 2))



if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags="ELLIPSIS")