__author__ = "Paul Schultz"

import numpy as np
import pandas as pd
from scipy.sparse import dok_matrix


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

    def new_node(self, BusID, time):
        # draw random location
        loc = np.random.random(2)
        return pd.DataFrame(data={"BusID": int(BusID), "lat": loc[0], "lon": loc[1], "time": time, "type":self.type},
        index=range(1))

class Coordinator(object):

    def __init__(self, **kwargs):
        pass

    def __str__(self):
        return "Coordinator class."

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



class PowerGrid(object):
    """
    derive from igraph?
    """

    def __init__(self, n_final):
        self.n_final = n_final
        self.adjacency = dok_matrix(np.zeros([self.n_final, self.n_final]))
        self.number_of_nodes = 0
        self.number_of_edges = 0

        # data containers for nodes and edges
        self.nodes = pd.DataFrame(columns=["BusID", "type", "lon", "lat", "time"])
        self.edges = pd.DataFrame(columns=["BranchID", "source", "target", "type", "time"])

        # internal settings
        self.distance_measure = "euclidean"

    def __str__(self):
        pass

    def add_node(self, new_node):
        # entry in node data table
        self.nodes = pd.merge(self.nodes, new_node, how="outer")
        # update counter
        self.number_of_nodes += 1
        #TODO: extend adjacency if number_of_nodes > len(adjacency)


    def add_edge(self, source, target, type, time):
        new_idx = self.number_of_edges + 1
        # entry in edge data table
        new_edge = pd.DataFrame(data={"BranchID": new_idx, "source": int(source), "target": int(target), "type": type,
                                      "time": time}, index=range(1))
        self.edges = pd.merge(self.edges, new_edge, how="outer")
        # add edge in adjacency matrix
        self.adjacency[int(source), int(target)] = \
            self.adjacency[int(target), int(source)] = \
            self._get_distance(int(source), int(target))
        # update counter
        self.number_of_edges += 1
        return new_idx

    def convert_to_igraph(self):
        import igraph as ig
        pass

    def update(self, new_node, edgelist):
        # add new node and all edges in edgelist

        self.add_node(new_node)

        if len(edgelist) == 0:
            pass
        else:
            for edge in edgelist:
                self.add_edge(edge[0], edge[1], "line", np.squeeze(new_node.time))

    def save(self, path, raw=False):
        import igraph as ig

        edgelist = zip(np.array(self.edges.source.values, dtype=np.int),
                       np.array(self.edges.target.values, dtype=np.int))

        g = ig.Graph()
        g.__init__(n=self.number_of_nodes, edges=edgelist, directed=False)

        if not raw:
            g.vs["BusID"] = self.nodes.BusID.values
            g.vs["time"] = self.nodes.time.values
            g.vs["type"] = self.nodes.type.values
            g.vs["lon"] = self.nodes.lon.values
            g.vs["lat"] = self.nodes.lat.values
            g.es["time"] = self.edges.time.values
            g.es["type"] = self.edges.type.values
            g.es["BranchID"] = self.edges.BranchID.values

        g.write_picklez(path)

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
    raise RuntimeWarning("Nothing will happen...")