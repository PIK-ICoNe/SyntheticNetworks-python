__author__ = "Paul Schultz, Jobst Heitzig"
__date__ = "Dec 15, 2016"
__version__ = "v2.2"

# This file is based on the network creation algorithm published in:
#
# A Random Growth Model for Power Grids and Other Spatially Embedded Infrastructure Networks
# Paul Schultz, Jobst Heitzig, and Juergen Kurths
# Eur. Phys. J. Special Topics on "Resilient power grids and extreme events" (2014)
# DOI: 10.1140/epjst/e2014-02279-6
#


# TODOs:
# reparamterize ns


import numpy as np
from igraph import Graph, load

import scipy.spatial.distance as sp
from rtree.index import Index as rtree # install via conda install --channel https://conda.anaconda.org/IOOS rtree

class RpgAlgorithm(object):


    def __init__(self, k):

        # parameters for the algorithm
        self.k = k

        self.w = [.945,.05,.005] # JH: list of relative frequencies of nodes by level
        self.n = [100,100,99800] #[2000,2000,60]
        self.n0 = [100,100,100] #[1250,250,50]
        self.p = [0,.1,.3]
        self.q = [0,.075,.075]
        self.r = [0.,0.7,1.4]
        self.s = [.2,.05,.0]

        self.sampling = "clotty"
        self.clottiness = 0.95

        # counters
        self.levnodes = [[] for l in range(self.k)] # list of nodes in level
        self.cumnodes = [[] for l in range(self.k)] # list of nodes in level or higher
        self.added_nodes = [0 for l in range(self.k)]
        self.added_edges = [0 for l in range(self.k)]
        self.noffset = 0 # total no. nodes added so far

        # node coordinates
        self.lon = []
        self.lat = []
        self.lev = []  # level of node

        # CHANGE WITH CAUTION!
        self.distance_measure = "euclidean"
        self.debug = False

        
    def __str__(self):
        print("----------")
        #print self.graph.num_vertices(), "nodes and", self.graph.num_edges(), "edges"
        for attr in vars(self):
            if attr in ["identifier", "added_nodes", "n", "n0", "p", "q", "r", "s", "k", "w"]:
                print(attr, ":", str(getattr(self, attr)))
        return "----------"


    ###############################################################################
    # ##                       PUBLIC FUNCTIONS                                ## #
    ###############################################################################


    def set_params(self, **kwargs):
        for key in kwargs:
            if not hasattr(self, key):
                print("ERROR: There is no parameter called:", key)
                print("Possible choices: n,n0,p,q,r,s")
                continue
            else:
                if self._validation(key, kwargs[key]):
                    setattr(self, key, kwargs[key])
                else:
                    print("ERROR: invalid parameter value for", key, kwargs[key])

    def prepare(self):
        """this should be called after set_params"""
        self.totaln = np.sum(self.n)

        self.setup_locations(sampling=self.sampling)

        self.mst_edges = [None for l in range(self.k)]
        self.init_edges = [None for l in range(self.k)]

        #TODO: find way to dynamically add nodes without index problems

        self.levgraph = [Graph(self.totaln) for l in range(self.k)] # one graph per level, storing only that level's edges
        self.cumgraph = [Graph(self.totaln) for l in range(self.k)] # one graph per level, storing that and all higher levels' edges

        self.levrtree = [None for l in range(self.k)] # one RTree per level, storing coordinates of level's nodes
        self.cumrtree = [None for l in range(self.k)] # one RTree per level, storing coordinates of this and higher level's nodes

    def initialise(self, l): # JH: l = level to initialise
        assert self.n[l] >= self.n0[l]

        # step I1: draw random locations from rho and add nodes
        #######################################################
        self._get_locations(l, self.noffset, self.noffset + self.n0[l], init=True)

        # step I2: construct minimum spanning tree
        ##########################################
        edge_mask = self._initial_mst(l) # will store only nodes and edges on this level

        if self.debug:
            print("I2", edge_mask)

        # step I3: add redundant links
        ##############################
        # CAUTION: logic has changed from original version! now it's simply the same as step G4, i.e., a little less optimal!
        
        m = min(int(np.floor(self.n0[l] * (1 - self.s[l]) * (self.p[l] + self.q[l]))), self.n0[l] * (self.n0[l] - 1) / 2 - (self.n0[l] - 1))

        for dummy in range(m):
            self._G34(l, np.random.choice(self.levnodes[l]))


    #        assert self.added_edges[l] == (len(self.adjacency[l].keys()) / 2)
        assert self.added_edges[l] == self.levgraph[l].ecount()
        assert self.added_nodes[l] == len(self.levnodes[l])

        # label initial edges
        self.init_edges[l] = self.levgraph[l].get_edgelist()
        

    def grow(self, lmax):
        """adds total no. of n[lmax] nodes to levels 0 <= l <= lmax"""

        new_nodes = range(self.noffset, self.n[lmax] - self.n0[lmax] + self.noffset)

        for node in new_nodes:

            # draw level for node:
            l = np.random.choice(range(lmax + 1), p=self.w[:lmax + 1] / np.sum(self.w[:lmax + 1]))
            self.lev.append(l)

            # register new node
            self._update_graphs(l, nodes=[node])

            if self.debug:
                print("---------")
                print("adding node", node, "to level", l)

            # step G5: split random link at midpoint
            ########################################
            if (np.random.random() < self.s[l]) and self.levgraph[l].ecount() > 0:
                self._G5(l, node)

            else:
                # step G2: link to nearest
                ##########################
                self._G2(l, node)

                # step G3: add optimal redundant link to node
                #############################################
                if np.random.random() < self.p[l]:
                    self._G34(l, node)

                # step G4: add another optimal redundant link to random node
                ############################################################
                if np.random.random() < self.q[l]:
                    self._G34(l, np.random.choice(self.levnodes[l]))

    def cleanup(self):
        """ remove objects from memory"""
        del self.levrtree
        del self.cumrtree
        for level in range(self.k):
            del self.levgraph[level]
            del self.cumgraph[level]



    def setup_locations(self, sampling="uniform", locations = None, centre=None, boundaries=None):
        """
        setup function that returns locations, either randomly or from data
        :param sampling:
        :param locations:
        :param centre:
        :param boundaries:
        :return:
        """
        if locations is not None:
            assert len(locations) == np.sum(self.n[:self.k])
            self.locations = locations
            self.counter = 0

        self.sampling = sampling
        self.centre = centre
        self.boundaries = boundaries

        # JH: docstring says this returns locations, but it returns nothing??
        

    ###############################################################################
    # ##                       PRIVATE FUNCTIONS                               ## #
    ###############################################################################

    def _get_coords(self, sampling=None, centre=None, boundaries=None):

        if sampling is not None:
            # override default sampling method
            self.sampling = sampling

        if self.sampling == "uniform":
            return self._uniformunitsquare(centre, boundaries)
        elif self.sampling == "data":
            pos = self.locations[self.counter]
            self.counter += 1
            return pos
        elif self.sampling == "clotty":
            l = len(self.lat)
            if l==0: return (0,0)
            i = np.random.choice(range(l))
            pos0 = np.array([self.lat[i], self.lon[i]])
            pos1 = pos0 + np.random.normal(size=2)
            b = np.random.choice([self.clottiness, 0])
            return tuple(b*pos0 + (1-b)*pos1)
        else:
            print("ERROR: Not implemented yet.")
            exit(1)


    def _get_distances(self, sources, targets):
        """
        return array of distances from nodes "sources" to list of nodes "targets"
        """
        x = np.c_[np.array(self.lon)[sources], np.array(self.lat)[sources]]
        y = np.c_[np.array(self.lon)[targets], np.array(self.lat)[targets]]
        return sp.cdist(x, y, metric=self.distance_measure)



    def _uniformunitsquare(self, centre=None, boundaries=None):
        """
        return point drawn uniformly at random

        :param centre: centre distribution around this point
        :param boundaries: array containing [width, height]
        :return: coordinate tuple
        """

        if centre is None:
            centre = -1.
        if boundaries is None:
            boundaries = -1.

        return (.5 - np.random.uniform(size=2)) * np.array(boundaries) + np.array(centre)


    def _G2(self, l, node):
        # only now get one new location and nearest earlier node:
        target = self._get_locations(l, node, node+1)

        if target is not None:
            # update graphs:
            d = self._get_distances([node], [target])[0, 0]
            self._update_graphs(l, edges=[(target, node)], weights=[d])

            if self.debug:
                print("G2", (node, target))

    def _G34(self, l, node):
        targets = list(set(self.cumnodes[l]).difference(self.cumgraph[l].neighbors(node)).difference([node]))
        if len(targets):
            dists = self._get_distances([node], targets)[0,:]
            prices = dists / (dists + self.cumgraph[l].shortest_paths_dijkstra(node, targets))**self.r[l] if self.r[l]>0 else dists
            best = np.argmin(prices)
            a, b = self._s((targets[best], node))

            # update graphs:
            d = dists[best]
            self._update_graphs(l, edges=[(a, b)], weights=[d])

            if self.debug:
                print("G3/4", (a, b))

    def _G5(self, l, node):
        # choose link at random:
        elist = self.levgraph[l].get_edgelist()
        a, b = elist[np.random.choice(range(len(elist)))]

        # NOTE: CHANGED BEHAVIOUR: now split somewhere, not in middle:
        pos = np.random.random()  # 0:a, 1:b

        # add node at midpoint and calc distances:
        lat = (1 - pos) * self.lat[a] + pos * self.lat[b]
        lon = (1 - pos) * self.lon[a] + pos * self.lon[b]
        self.lat.append(lat)
        self.lon.append(lon)

        # update graphs and rtrees:

        eid = self.levgraph[l].get_eid(a, b)
        d = self.levgraph[l].es["weight"][eid]

        self.levrtree[l].insert(node, (lat, lon, lat, lon))
        for l2 in range(l + 1):
            self.cumrtree[l2].insert(node, (lat, lon, lat, lon))

        self._update_graphs(l, edges=[(a, b)], delete_edges=True)
        self._update_graphs(l, edges=[(a, node), (b, node)], weights=[pos * d, (1 - pos) * d])

        if self.debug:
            print("G5", (int(a), int(b)))

    def _validation(self, attr, value):
        value = np.array(value)
        if attr == "n0" or attr == "n":
            if any(value < 1):
                return False
            else:
                return True
        elif attr == "r":
            if any(value < 0):
                return False
            else:
                return True
        elif attr in ["p", "q", "s", "w"]:
            if any(value < 0) or any(value > 1):
                return False
            else:
                return True
        elif attr == "k":
            if value < 1:
                return False
            else:
                return True
        elif attr == "clottiness":
            if value < 0 or value > 1:
                return False
            else:
                return True


    def _initial_mst(self, l):

        self.lev += [l for i in range(self.n0[l])]
        nodes = list(range(self.noffset, self.noffset+self.n0[l]))
        self.mst_edges[l] = elist = self._get_mst(l)
        self._update_graphs(l, nodes=nodes, edges=elist)

        return elist


    def _get_mst(self, l):
        nodes = range(self.noffset, self.noffset + self.n0[l])
        distmatrix = self._get_distances(nodes, nodes)
        full_graph = Graph.Full(self.n0[l])
        factor = 1e5  # since small weights lead to MST problems
        weights = [factor * distmatrix[i,j] for (i,j) in full_graph.get_edgelist()]
        G = full_graph.spanning_tree(weights).as_undirected()
        return [self._s((i+self.noffset,j+self.noffset)) for (i,j) in G.get_edgelist()]


    def _get_locations(self, l, offset, _m, init=False):
        m = int(_m)
        poss = np.zeros((m,2))
        for i in range(offset, m):
            poss[i,:] = pos = self._get_coords(self.sampling, self.centre, self.boundaries)
            self.lat.append(pos[0])
            self.lon.append(pos[1])
            # update earlier rtree spatial indices:
            for l2 in range(l):
                self.cumrtree[l2].insert(i, (pos[0],pos[1],pos[0],pos[1]))
            if not init: # otherwise en bulk (below)
                nearest = list(self.cumrtree[l].nearest((pos[0],pos[1],pos[0],pos[1]),1))[0] if m > 0 else None # query before adding!
                self.levrtree[l].insert(i, (pos[0],pos[1],pos[0],pos[1]))
                self.cumrtree[l].insert(i, (pos[0],pos[1],pos[0],pos[1]))                
#        self._update_distance(offset, m, m)
        if init: # bulk insert: # TODO: CAUTION: must only be used at initialization of level!
            # set up additional rtree spatial indices:
            def f():
                for i in range(offset, m):
                    yield (i, (poss[i,0],poss[i,1],poss[i,0],poss[i,1]),None)
            self.levrtree[l] = lrt = rtree(f())
            self.cumrtree[l] = crt = rtree(f()) # sadly, rtrees cannot be cloned yet
        else:
            return nearest


    def _update_counters(self, level, nodes=0, edges=0):
        self.added_nodes[level] += nodes
        self.noffset += nodes
        self.added_edges[level] += edges

    def _update_graphs(self, level, nodes=[], edges=[], weights=[], delete_edges=False):
        if delete_edges:
            eid = self.levgraph[level].get_eids(edges)
            self.levgraph[level].delete_edges(eid)

            for l in range(level + 1):
                eid = self.cumgraph[l].get_eids(edges)
                self.cumgraph[l].delete_edges(eid)

            self._update_counters(level, edges=-1)
        else:
            if nodes:
                # PS: this is not necessary, as all Graphs are created with size totaln.
                # otherwise, the difference between index and name is going to cause many problems
                # self.levgraph[level].add_vertices(nodes)
                # for l in range(level + 1):
                #     self.cumgraph[l].add_vertices(nodes)

                self.levnodes[level].extend(nodes)
                for l in range(level + 1):
                    self.cumnodes[l].extend(nodes)

            if edges:
                if not weights:
                    weights = [self._get_distances([i],[j])[0,0] for (i,j) in edges]

                for idx, (i, j) in enumerate(edges):
                    # level graphs do not contain links between levels,
                    #if self.lev[i] == self.lev[j]:
                    self.levgraph[level].add_edge(i, j, weight=weights[idx])

                    for l in range(level + 1):
                        self.cumgraph[l].add_edge(i, j, weight=weights[idx])

            self._update_counters(level, nodes=len(nodes), edges=len(edges))


    def _s(self, tuple):
        if tuple[0] < tuple[1]:
            return tuple
        else:
            return (tuple[1], tuple[0])


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


def calc(name="test"):

    np.random.seed(0)

    # initialise algorithm
    g = RpgAlgorithm(k=3)
    assert(isinstance(g, RpgAlgorithm))

    # for detailed output set 
    g.debug = True

    branching = np.array([6084.,84.,2.])

    # set desired parameters and perform algorithm
    g.set_params(k=1,
                 n=[1000],
                 n0=[10],
                 w=[1.],
                 p=[0.2],
                 q=[0.75],
                 r=[0.6],
                 s=[.2],
                 clottiness=0.95
                 )

    # use predefined locations ...
    # g.setup_locations(sampling="data", locations=np.random.random([g.n, 2]))

    g.prepare()
    for l in range(g.k):
        g.initialise(l)
        g.grow(l)


    print g

    G = g.cumgraph[0].copy()

    elist = np.array(G.get_edgelist())
    G.es['level'] = map(lambda (a,b): min(g.lev[a], g.lev[b]), elist)
    G.vs["level"] = g.lev
    G.vs["lat"] = g.lat
    G.vs["lon"] = g.lon

    G.write_pickle(name)

    return G

def plot(G=None, name="test", groups=False):

    if G is None:
        G = Graph.Read_Pickle(name)

    cols = {0: "grey", 1: "blue", 2: "red"}
    weights = {0: 1, 1: 1.5, 2: 2}
    sizes = weights

    G.vs['color'] = map(lambda y: cols[y], G.vs["level"])
    G.es['color'] = map(lambda y: cols[y], G.es["level"])
    G.es['width'] = map(lambda y: 10. * weights[y], G.vs["level"])
    G.vs['size'] = map(lambda y: 20. * sizes[y], G.vs["level"])

    print "connected graph:", G.is_connected()

    from igraph import plot, Layout
    l = [(xy[0], xy[1]) for xy in np.array([G.vs["lat"], G.vs["lon"]]).T]

    w = 2000 * np.sqrt(G.vcount())
    if groups:
        comp = G.clusters()
        sort = np.argsort([len(c) for c in comp])[::-1]
        comp = [comp[i] for i in sort]
        cmap = np.tile(["grey", "blue", "red", "yellow"], 3)
        group_markers = []

        print "components:", len(comp)
        for i, c in enumerate(comp):
            if i >= len(cmap):
                break
            print i, len(c), cmap[i]
            group_markers.append((c, cmap[i]))

        plot(G, "output.pdf",
             bbox=(w, w),
             layout=Layout(coords=l),
             vertex_order=np.argsort(G.vs["level"]),
             mark_groups=group_markers
             )
    else:
        plot(G, "output.pdf",
             bbox=(w, w),
             layout=Layout(coords=l),
             vertex_order=np.argsort(G.vs["level"])
             )

def hist(G=None, name="test"):

    if G is None:
        G = Graph.Read_Pickle(name)
        assert isinstance(G, Graph)

    import pandas as pd
    import matplotlib.pyplot as plt
    import scipy.stats as st

    df = pd.DataFrame({"length": G.es["weight"], "loglength":np.log10(G.es["weight"]), "level": G.es["level"]})

    df.pivot(columns="level", values="loglength").plot(kind="hist", bins=40, stacked=True, log=True, grid=True)
    plt.xlabel(r"$\log_{10}$ length")
    plt.savefig(name + "_loglength_dist.png")


    # print "aspl", G.average_path_length(), "transitivity", G.transitivity_undirected()

    df_nodes = pd.DataFrame({"level": G.vs["level"],
                             "degree": G.vs.degree(),
                             "clust": G.transitivity_local_undirected(),
                             "betw": 2. * np.array(G.betweenness()) / (G.vcount() * (G.vcount() - 1.))
                             })

    df_nodes = df_nodes.pivot(columns="level")

    df_nodes.degree.plot(kind="hist", bins=40, stacked=True, log=True, grid=True)
    plt.xlabel("degree")
    plt.savefig(name + "_degree_dist.png")

    df_nodes.clust.plot(kind="hist", bins=40, stacked=True, log=True, grid=True)
    plt.xlabel("local transitivity")
    plt.savefig(name + "_clust_dist.png")

    df_nodes.betw.plot(kind="hist", bins=40, stacked=True, log=True, grid=True)
    plt.xlabel("shortest path betweenness")
    plt.savefig(name + "_betw_dist.png")

if __name__ == "__main__":
    G = calc()
    plot(G, groups=False)
    hist(G)


