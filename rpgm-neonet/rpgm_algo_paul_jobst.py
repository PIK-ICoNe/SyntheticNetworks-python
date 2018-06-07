__author__ = "Paul Schultz, Jobst Heitzig"
__date__ = "Dec 12, 2016"
__version__ = "v2.1b"

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
from scipy.sparse import dok_matrix
from igraph import Graph

from rtree.index import Index as rtree # install via conda install --channel https://conda.anaconda.org/IOOS rtree

xrange = range




class RpgAlgorithm(object):


    def __init__(self):

        # parameters for the algorithm
        self.k = 3
        self.w = [.945,.05,.005] # JH: list of relative frequencies of nodes by level
        self.n = [1,1,9998] #[2000,2000,60]
        self.n0 = [1,1,100] #[1250,250,50]
        self.p = [0,.1,.3]
        self.q = [0,.075,.075]
        self.r = [0.,.1,.4]
        self.s = [.2,.05,.0]
        
        # node coordinates
        self.lev = [] # level of node
        self.levnodes = [[] for l in range(self.k)] # list of nodes in level
        self.cumnodes = [[] for l in range(self.k)] # list of nodes in level or higher
        self.lon = []
        self.lat = []
        self.added_nodes = [0 for l in range(self.k)]
        self.added_edges = [0 for l in range(self.k)]
        self.noffset = 0 # total no. nodes added so far

        self.totaln = np.sum(self.n)
        # JH: this was the performance bottleneck, so I sped it up.
#        self.distance = {}
#        self.distmatrix = np.zeros((totaln,totaln))

        self.setup_locations()

        # CHANGE WITH CAUTION!
        self.distance_measure = "euclidean"
        self.low_memory = False
        self.debug = False
        
        self.adjacency = [None for l in range(self.k)]
        self.mst_edges = [None for l in range(self.k)]
        self.init_edges = [None for l in range(self.k)]
#        self.dGmatrix = [None for l in range(self.k)]

        self.levgraph = [None for l in range(self.k)] # one graph per level, storing only that level's edges
        self.cumgraph = [None for l in range(self.k)] # one graph per level, storing that and all higher levels' edges
        
        self.levrtree = [None for l in range(self.k)] # one RTree per level, storing coordinates of level's nodes
        self.cumrtree = [None for l in range(self.k)] # one RTree per level, storing coordinates of this and higher level's nodes

        
    def __str__(self):
        print("----------")
        #print self.graph.num_vertices(), "nodes and", self.graph.num_edges(), "edges"
        for attr in vars(self):
            if attr in ["identifier", "added_nodes", "n", "n0", "p", "q", "r", "s"]:
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


    def initialise(self, l): # JH: l = level to initialise
        assert self.n[l] >= self.n0[l]

        # keep track of added nodes
        self.added_nodes[l] = 0

        # step I1: draw random locations from rho and add nodes
        #######################################################
        self._get_locations(l, self.noffset, self.noffset + self.n0[l], init=True)
        self.added_nodes[l] += self.n0[l]

        # step I2: construct minimum spanning tree
        ##########################################
        edge_mask = self._initial_mst(l) # will store only nodes and edges on this level

        if self.debug:
            print("I2", edge_mask)

        # step I3: add redundant links
        ##############################
        # CAUTION: logic has changed from original version! now it's simply the same as step G4, i.e., a little less optimal!
        
        m = min(int(np.floor(self.n0[l] * (1 - self.s[l]) * (self.p[l] + self.q[l]))), self.n0[l] * (self.n0[l] - 1) / 2 - (self.n0[l] - 1))

        for dummy in xrange(m): 
            self._G34(l,np.random.choice(self.levnodes[l]))

        assert self.added_edges[l] == (len(self.adjacency[l].keys()) / 2)

        # label initial edges
        self.init_edges[l] = self.levgraph[l].get_edgelist()
        

    def grow(self, lmax):
        """adds total no. of n[lmax] nodes to levels 0 <= l <= lmax"""
        
        nnew = self.n[lmax] - self.n0[lmax]
        for l in range(lmax + 1):
            # extend existing adjacency matrices:
            self.adjacency[l]._shape = (self.n[lmax] + self.noffset, self.n[lmax] + self.noffset)
            # get distances:
#             if self.r[l] > 0:
#                 self.dGmatrix[l] = self._get_graph_distances(l,lmax) # JH: why is this too small and needs to be extended?
#                 self.dGmatrix[l] = np.concatenate([np.concatenate([self.dGmatrix[l], 
#                                                                    np.inf*np.ones((self.n[lmax] - self.n0[lmax], self.n0[lmax] + self.noffset))],axis=0),
#                                                    np.concatenate([np.inf*np.ones((self.n0[lmax] + self.noffset, self.n[lmax] - self.n0[lmax])),
#                                                                    np.inf*(np.ones((nnew,nnew)) - np.eye(nnew,nnew))],axis=0)],axis=1)
#                 np.fill_diagonal(self.dGmatrix[l], 0)
                
        # connect new vertices
        for i in xrange(self.n0[lmax] + self.noffset, self.n[lmax] + self.noffset):
            self._growth_step(i, lmax)

        self.noffset += self.n[l]


#     def edgelist(self, l):
#         return sorted(set([self._s(key) for key in self.adjacency[l].iterkeys()]))


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
        else:
            print("ERROR: Not implemented yet.")
            exit(1)


#     def _get_distance(self, u, v):
#         if self.distance_measure == "euclidean":
#             return self._euclidean(int(u), int(v))
#         else:
#             print("ERROR: Not implemented yet.")
#             exit(1)

    def _get_distances(self, sources, targets):
        """
        return array of distances from node "u" to list of nodes "targets"
        """
        assert self.distance_measure == "euclidean"
        x = np.array([self.lat, self.lon])
        return np.sqrt(np.sum((x[:,sources].reshape((2,-1,1)) - x[:,targets].reshape((2,1,-1)))**2,axis=0))
        

#     def _update_distance(self, v1, v2, vmax):
#         print("updating distances from",v1,"to",v2-1)
# #        for v in xrange(v1, v2):
# #            for u in xrange(v):
# #                self.distance[(u, v)] = self._get_distance(u, v)
# #            for u in xrange(v+1,vmax):
# #                self.distance[(v, u)] = self._get_distance(v, u)
#         assert self.distance_measure == "euclidean"
#         x = np.array([self.lat, self.lon])
#         self.distmatrix[v1:v2,:vmax] = d = np.sqrt(
#                 np.sum((x[:,v1:v2].reshape((2,-1,1)) - x[:,:vmax].reshape((2,1,-1)))**2,axis=0)
#                 )
#         self.distmatrix[:vmax,v1:v2] = d.T

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


#     def _euclidean(self, u, v):
#         """
#         return euclidean distance between x and y
#         """
#         x = np.array([self.lat[u], self.lon[u]])
#         y = np.array([self.lat[v], self.lon[v]])
#         return np.sqrt(sum((x-y)**2))


    def _growth_step(self, node, lmax):
        # draw level for node:
        l = np.random.choice(range(lmax + 1),p=self.w[:lmax+1]/np.sum(self.w[:lmax+1]))
        self.lev.append(l)
        self.levnodes[l].append(node)
        for l2 in range(l+1):
            self.cumnodes[l2].append(node)
        self.added_nodes[l] += 1

        if self.debug:
            print("---------")
            print("adding node", node, "to level", l)

        # step G5: split random link at midpoint
        ########################################
        if (np.random.random() < self.s[l]) and len(self.adjacency[l]) > 1:
            # choose link at random:
            elist = self.levgraph[l].get_edgelist()
            a,b = elist[np.random.choice(range(len(elist)))]
            
            # NOTE: CHANGED BEHAVIOUR: now split somewhere, not in middle:
            pos = np.random.random() # 0:a, 1:b

            if self.debug:
                print("G5", (int(a), int(b)))

            # add node at midpoint and calc distances:
            if self.debug:
                print(node,a,b,len(self.lat))
            lat = (1-pos)*self.lat[a] + pos*self.lat[b]
            lon = (1-pos)*self.lon[a] + pos*self.lon[b]
            self.lat.append(lat)
            self.lon.append(lon) 
#             if not self.low_memory:
#                 self._update_distance(node,node+1,self.noffset+self.n[lmax])
            self.adjacency[l][a, b] = self.adjacency[l][b, a] = 0
            self.adjacency[l][a, node] = self.adjacency[l][node, a] = 1
            self.adjacency[l][b, node] = self.adjacency[l][node, b] = 1
            
            # update graphs and rtrees:
            eid = self.levgraph[l].get_eid(a,b)
            d = self.levgraph[l].es["weight"][eid]
            self.levgraph[l].delete_edges(eid)
            self.levgraph[l].add_edge(a, node, weight = pos*d)
            self.levgraph[l].add_edge(b, node, weight = (1-pos)*d)
            self.levrtree[l].insert(node, (lat,lon,lat,lon))
            for l2 in range(l+1):
                eid = self.cumgraph[l2].get_eid(a,b)
                self.cumgraph[l2].delete_edges(eid)
                self.cumgraph[l2].add_edge(a, node, weight = pos*d)
                self.cumgraph[l2].add_edge(b, node, weight = (1-pos)*d)
                self.cumrtree[l2].insert(node, (lat,lon,lat,lon))

#             d = self.distmatrix[a,b]
#             for l2 in range(l+1):
#                 if self.r[l2] > 0:
# #                    self.dGmatrix[l2][:(node + 1), :(node + 1)] = self._get_graph_distances(l2,lmax)
#                     self.dGmatrix[l2][:node,node] = self.dGmatrix[l2][node,:node] = \
#                         np.minimum(pos*d + self.dGmatrix[l2][a,:node], (1-pos)*d + self.dGmatrix[l2][b,:node])

            self.added_edges[l] += 1

        else:
            # only now get one new location and nearest earlier node:
            target = self._get_locations(l, node, node+1) 
            if target is not None:
    
                # step G2: link to nearest
                ##########################
    #             if node == 1:
    #                 target = 0
    #             else:
    #                 target = self._get_closest_connected_node(node, node - 1, l)
                self.adjacency[l][node, target] = self.adjacency[l][target, node] = 1
                # update graphs:
                d = self._get_distances([node], [target])[0,0]
                self.levgraph[l].add_edge(target, node, weight = d)
                for l2 in range(l+1):
                    self.cumgraph[l2].add_edge(target, node, weight = d)
    
    #             for l2 in range(l+1):
    #                 if self.r[l2] > 0:
    #                     # adjust network distances:
    #                     #self.dGmatrix[l2][:(node + 1), :(node + 1)] = self._get_graph_distances(l2,lmax)
    # #                    self.dGmatrix[l2][node, :(node + 1)] = self.dGmatrix[l2][target, :(node + 1)] + 1
    # #                    self.dGmatrix[l2][:(node + 1), node] = self.dGmatrix[l2][:(node + 1), target] + 1
    #                     d = self.distmatrix[node,target]
    #                     self.dGmatrix[l2][node, :(node + 1)] = self.dGmatrix[l2][target, :(node + 1)] + d
    #                     self.dGmatrix[l2][:(node + 1), node] = self.dGmatrix[l2][:(node + 1), target] + d
    #                     self.dGmatrix[l2][node, node] = 0
    
                self.added_edges[l] += 1
    
                if self.debug:
                    print("G2", (node, target))
    
                # step G3: add optimal redundant link to node
                #############################################
                if np.random.random() < self.p[l]:
                    self._G34(l,node)
    
                # step G4: add another optimal redundant link to random node
                ############################################################
                if np.random.random() < self.q[l]: 
                    self._G34(l,np.random.choice(self.levnodes[l]))
    
    #         if self.debug and self.r[l] > 0:
    #             # check that update went well
    # #            assert np.all(self.dGmatrix[l][:(node + 1), :(node + 1)] == self._get_graph_distances(l,lmax))
    #             assert not np.any(abs(self.dGmatrix[l][:(node + 1), :(node + 1)]-self._get_graph_distances(l,lmax))>1e-10) # check that update went well


    def _G34(self,l,node):
#                 candidates = {}
#                 for v in xrange(node - 1):
#                     if self.adjacency[l][v, node] == 0 and self.lev[v] >= l: # JH: cannot connect to lower-level nodes
# #                        candidates[(v, node)] = self._get_distance(v, node) if self.low_memory else self.distance[(v, node)]
#                         candidates[(v, node)] = self._get_distance(v, node) if self.low_memory else self.distmatrix[v,node]
# 
#                 # there might be no candidates if n0 = 1
#                 if len(candidates) > 0:
#                     if self.r[l] > 0:
#                         #dGmatrix = self._get_graph_distances()
# 
#                         for (u, v) in candidates.iterkeys():
#                             candidates[(u, v)] /=  ( 1. + self.dGmatrix[l][u, v])**self.r[l]
# 
#                     a, b = min(candidates, key=candidates.get)

                targets = list(set(self.cumnodes[l]).difference(self.cumgraph[l].neighbors(node)).difference([node]))
                if len(targets):
                    dists = self._get_distances([node], targets)[0,:]
                    prices = dists / (dists + self.cumgraph[l].shortest_paths_dijkstra([node], targets))**self.r[l] if self.r[l]>0 else dists 
                    best = np.argmin(prices)
                    a,b = self._s((targets[best],node))
                    
                    self.adjacency[l][a, b] = self.adjacency[l][b, a] = 1
                    # update graphs:
                    d = dists[best]
                    self.levgraph[l].add_edge(a, b, weight = d)
                    for l2 in range(l+1):
                        self.cumgraph[l2].add_edge(a, b, weight = d)

#                     d = self.distmatrix[a,b] * np.ones((self.n[lmax] + self.noffset, self.n[lmax] + self.noffset))
#                     if self.r[l] > 0:
#                         if a == node:
#                             target = b
#                         else:
#                             target = a
#                         # adjust network distances:
#                         for l2 in range(l+1):
#                             if self.r[l2] > 0:
# #                                self.dGmatrix[l2][:(node + 1), :(node + 1)] = self._get_graph_distances(l2)
# #                                self.dGmatrix[l2] = np.minimum(
# #                                     np.minimum(self.dGmatrix[l2],
# #                                                self.dGmatrix[l2][:,[node]] +
# #                                                np.ones((self.n[lmax] + self.noffset, self.n[lmax] + self.noffset)) +
# #                                                self.dGmatrix[l2][[target],:]
# #                                                ),
# #                                     self.dGmatrix[l2][:,[target]] +
# #                                     np.ones((self.n[lmax] + self.noffset, self.n[lmax] + self.noffset)) +
# #                                     self.dGmatrix[l2][[node],:]
#                                 da = self.dGmatrix[l2][:,[a]]
#                                 db = self.dGmatrix[l2][:,[b]]
#                                 self.dGmatrix[l2] = np.minimum(np.minimum(self.dGmatrix[l2], da + d + db.T), db + d + da.T) 

                    self.added_edges[l] += 1

                    if self.debug:
                        print("G3/4", (a, b))

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


    def _initial_mst(self, l):
        adjacency = np.zeros([self.n0[l] + self.noffset, self.n0[l] + self.noffset])
        np.fill_diagonal(adjacency, 0)

        self.mst_edges[l] = elist = self._get_mst(l)
        for (i,j) in elist:
            adjacency[i,j] = adjacency[j,i] = 1
            self.added_edges[l] += 1

        self.adjacency[l] = dok_matrix(adjacency)

        self.lev += [l for i in range(self.n0[l])]
        self.levnodes[l] = nodes = list(range(self.noffset,self.noffset+self.n0[l]))
        self.levgraph[l] = lg = Graph(self.totaln)
        self.cumgraph[l] = cg = Graph(self.totaln) # will store nodes and edges on this and higher levels
        lg.add_edges(elist)     
        cg.add_edges(elist)     
        for l2 in range(l):
            self.cumnodes[l] += nodes
            self.cumgraph[l2].add_edges(elist)

        return elist


    def _get_mst(self, l):
        full_graph = Graph.Full(self.n0[l])
        factor = 1e5 # since small weights lead to MST problem
        nodes = range(self.noffset,self.noffset+self.n0[l])
        distmatrix = self._get_distances(nodes, nodes)
        weights = [factor * distmatrix[i,j] for (i,j) in full_graph.get_edgelist()]
        G = full_graph.spanning_tree(weights).as_undirected()
        return [self._s((i+self.noffset,j+self.noffset)) for (i,j) in G.get_edgelist()]


#     def _get_graph_distances(self, l0, lmax):
#         """distances using edges of levels l0...lmax"""
#         elist = []
#         for l in range(l0, lmax + 1):
#             elist += self.edgelist(l)
#         G = Graph(np.sum(self.added_nodes))
#         G.add_edges(elist)
# #        return np.array(G.shortest_paths())
#         return np.array(G.shortest_paths_dijkstra(weights=[self.distmatrix[i,j] for (i,j) in elist]))


    def _get_locations(self, l, offset, _m, init=False):
        m = int(_m)
        poss = np.zeros((m,2))
        for i in xrange(offset, m):
            poss[i,:] = pos = self._get_coords(self.sampling, self.centre, self.boundaries)
            self.lat.append(pos[0])
            self.lon.append(pos[1])
            # update earlier rtree spatial indices:
            for l2 in range(l):
                self.cumrtree[l2].insert(i, (pos[0],pos[1],pos[0],pos[1]))
            if not init: # otherwise en bulk (below)
                nearest = list(self.cumrtree[l].nearest((pos[0],pos[1],pos[0],pos[1]),1))[0] if m>0 else None # query before adding!
                self.levrtree[l].insert(i, (pos[0],pos[1],pos[0],pos[1]))
                self.cumrtree[l].insert(i, (pos[0],pos[1],pos[0],pos[1]))                
#        self._update_distance(offset, m, m)
        if init: # bulk insert: # TODO: CAUTION: must only be used at initialization of level!
            # set up additional rtree spatial indices:
            def f():
                for i in xrange(offset, m):
                    yield (i, (poss[i,0],poss[i,1],poss[i,0],poss[i,1]),None)
            self.levrtree[l] = lrt = rtree(f())
            self.cumrtree[l] = crt = rtree(f()) # sadly, rtrees cannot be cloned yet
        else:
            return nearest



    def _get_closest_connected_node(self, source, connected, l):
        """return closest connected node in level >= l"""
        # vertices need to be properly ordered for this to work, i.e. nodes in the connected component
        # should be labeled from 0 to connected-1
        min = np.inf
        target = source
        for node in xrange(connected):
            if source == node:
                # should actually not happen
                continue
            elif self.adjacency[l][source, node] == 0 and self.lev[node] >= l:
#                d = self._get_distance(node, source) if self.low_memory else self.distance[(node, source)]
                d = self._get_distance(node, source) if self.low_memory else self.distmatrix[node, source]
                if d < min:
                    min = d
                    target = node
        return target

    def _s(self, tuple):
        if tuple[0] < tuple[1]:
            return tuple
        else:
            return (tuple[1], tuple[0])


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


def main():

    # initialise algorithm
    g = RpgAlgorithm()
    assert(isinstance(g, RpgAlgorithm))

    # for detailed output set 
    g.debug = True

    # set desired parameters and perform algorithm
    g.set_params(k=3, n=[1000,1000,1000], n0=[1,1,10])

    # use predefined locations ...
    # g.setup_locations(sampling="data", locations=np.random.random([g.n, 2]))

    for l in range(g.k):
        g.initialise(l)
        g.grow(l)


    print(g)
    # print g.adjacency

    # create igraph object
    cols = {0:"grey",1:"blue",2:"red"}
    weights = {0:1,1:1.5,2:2}
    sizes = weights
    G = Graph(np.sum(g.n[:g.k]))
    for l in range(g.k):
        el = g.levgraph[l].get_edgelist()
        G.add_edges(el)
        for i,j in el:
            eid = G.get_eid(i,j)
            G.es[eid]['color'] = cols[l] 
            G.es[eid]['width'] = 10*weights[l] 
        # TODO: store color and thickness for edges!
    G.vs["lat"] = g.lat
    G.vs["lon"] = g.lon
    G.vs["color"] = [cols[l] for l in g.lev]
    G.vs["size"] = [20*sizes[l] for l in g.lev]

#    from igraph.drawing.graph import DefaultGraphDrawer
#    GD = DefaultGraphDrawer()
#    GD.draw(G)

    from igraph import plot, Layout
    l = [(xy[0],xy[1]) for xy in np.array([g.lat,g.lon]).T]
    
    G.write("test2.dot")
    
    plot(G, "output2.pdf", bbox=(8000, 8000), layout=Layout(coords=l))


    return G


if __name__ == "__main__":
    main()



