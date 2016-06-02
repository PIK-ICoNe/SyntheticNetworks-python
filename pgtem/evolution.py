import pandas as pd
from actors import *

class PowerGrid(object):
    """
    derive from igraph?
    """

    def __init__(self):
        self.nodes = pd.DataFrame(columns=["BusID", "type", "loc", "time"])
        self.edges = pd.DataFrame(columns=["BranchID", "source", "target", "type", "time"])
        self.number_of_nodes = 0
        self.number_of_edges = 0

    def __str__(self):
        pass

    def add_node(self, type, loc, time):
        new_idx = self.number_of_nodes + 1
        self.nodes.loc[new_idx] = [int(new_idx), type, loc, time]
        self.number_of_nodes += 1
        return new_idx

    def add_edge(self, source, target, type, time):
        new_idx = self.number_of_edges + 1
        self.edges.loc[new_idx] = [int(new_idx), source, target, type, time]
        self.number_of_edges += 1
        return new_idx

    def convert_to_igraph(self):
        import igraph as ig
        pass

    def update(self, new_loc, edgelist, time):
        # add new node and all edges in edgelist

        self.add_node("source", new_loc, time)

        for edge in edgelist:
            self.add_edge(edge[0], edge[1], "line", time)

def node_proposition(new_loc, coordinator, providers, regulator):

    assert isinstance(coordinator, Coordinator)
    assert [isinstance(provider, Provider) for provider in providers]
    assert isinstance(regulator, Regulator)

    # coordinator proposes possible extensions to connect new node
    possible_extensions = coordinator.proposition(new_loc)

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


if __name__ == '__main__':
    raise RuntimeWarning("Nothing will happen...")