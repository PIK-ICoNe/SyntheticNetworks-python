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



class RpgAlgorithm(object):
    def __init__(self):

        # parameters for the algorithm
        self.n = 20
        self.n0 = 10
        self.p = 0.2
        self.q = 0.3
        self.r = 1. / 3.
        self.s = 0.1

        # node coordinates
        self.lon = []
        self.lat = []
        self.distance = {}
        self.added_nodes = 0
        self.added_edges = 0

        # CHANGE WITH CAUTION!
        self.distance_measure = "euclidean"
        self.sampling = "uniform"
        self.low_memory = True
        self.debug = False

    def __str__(self):
        print "----------"
        #print self.graph.num_vertices(), "nodes and", self.graph.num_edges(), "edges"
        for attr in vars(self):
            if attr in ["identifier", "added_nodes", "n", "n0", "p", "q", "r", "s"]:
                print attr, ":", str(getattr(self, attr))
        return "----------"

    ###############################################################################
    # ##                       PUBLIC FUNCTIONS                                ## #
    ###############################################################################

    def set_params(self, **kwargs):
        for key in kwargs:
            if not hasattr(self, key):
                print "ERROR: There is no parameter called:", key
                print "Possible choices: n,n0,p,q,r,s"
                continue
            else:
                if self._validation(key, kwargs[key]):
                    setattr(self, key, kwargs[key])
                else:
                    print "ERROR: invalid parameter value for", key

    def initialise(self):
        assert self.n >= self.n0

        # keep track of added nodes
        self.added_nodes = 0

        # step I1: draw random locations from rho and add nodes
        #######################################################
        self._add_random_locations(self.n0)
        self.added_nodes += self.n0

        # step I2: construct minimum spanning tree
        ##########################################
        self._initial_mst()

        edge_mask = sorted(set([self._s(key) for key in self.adjacency.iterkeys()]))

        if self.debug:
            print "I2", edge_mask


        # step I3: add redundant links
        ##############################
        m = min(int(np.floor(self.n0 * (1 - self.s) * (self.p + self.q))), self.n0 * (self.n0 - 1) / 2 - (self.n0 - 1))

        candidates = {}
        for (u, v) in self.distance.iterkeys():
            if not (u, v) in edge_mask:
                candidates[(u, v)] = self.distance[(u, v)]

        if self.r > 0:
            dGmatrix = self._get_graph_distances()
            onesquare = np.ones([self.n0, self.n0])

        for k in xrange(m):
            if self.r > 0:
                #dGmatrix = self._get_graph_distances()

                for (u, v) in candidates.iterkeys():
                    candidates[(u, v)] = self.distance[(u, v)] / ( 1. + dGmatrix[u, v])**self.r

            a, b = min(candidates, key=candidates.get)
            self.adjacency[a, b] = self.adjacency[b, a] = 1
            # make sure i,j is not selected again:
            candidates.pop((a, b), None)

            if self.r > 0:
                dGmatrix = np.minimum(np.minimum(dGmatrix, dGmatrix[:,[a]] + onesquare + dGmatrix[[b],:]),
                                      dGmatrix[:,[b]] + onesquare + dGmatrix[[a],:])

            if self.debug:
                print "I3", (a, b)

        self.added_edges += m

        assert self.added_edges == (len(self.adjacency.keys()) / 2)

        # label initial edges
        self.init_edges = sorted(set([self._s(key) for key in self.adjacency.iterkeys()]))

        if self.debug and self.r > 0:
            assert ((dGmatrix - self._get_graph_distances())**2).sum() == 0 # check that update went well


    def grow(self):
        self._add_random_locations(self.n - self.n0)
        self.adjacency._shape = (self.n, self.n)

        if self.r > 0:
            self.dGmatrix = self._get_graph_distances()
            self.dGmatrix = np.concatenate([np.concatenate([self.dGmatrix, np.zeros((self.n - self.n0, self.n0))],axis=0),
                                            np.zeros((self.n, self.n - self.n0))],axis=1)
        # connect new vertices
        for i in xrange(self.n0, self.n):
            self.added_nodes += 1
            self._growth_step(i)


        # TODO: this is probably redundant
        assert self.added_nodes == self.n

        if self.debug and self.r > 0:
            assert ((self.dGmatrix - self._get_graph_distances())**2).sum() == 0 # check that update went well

        if self.r > 0:
            delattr(self, "dGmatrix")



    ###############################################################################
    # ##                       PRIVATE FUNCTIONS                               ## #
    ###############################################################################

    def _get_coords(self):
        if self.sampling == "uniform":
            return self._uniformunitsquare()
        else:
            print "ERROR: Not implemented yet."
            exit(1)

    def _get_distance(self, u, v):
        if self.distance_measure == "euclidean":
            return self._euclidean(int(u), int(v))
        else:
            print "ERROR: Not implemented yet."
            exit(1)

    def _update_distance(self):
        N = len(self.lat)
        for v in xrange(N):
            for u in xrange(v):
                self.distance[(u, v)] = self._get_distance(u, v)

    def _uniformunitsquare(self):
        """
        return point drawn uniformly at random
        from unit square -> 2D coordinates
        """
        return np.random.uniform(size=2)

    def _euclidean(self, u, v):
        """
        return euclidean distance between x and y
        """
        x = np.array([self.lat[u], self.lon[u]])
        y = np.array([self.lat[v], self.lon[v]])
        return np.sqrt(sum((x-y)**2))

    def _growth_step(self, node):
        if self.debug:
            print "---------"
            print "adding node", node




        # step G5: split random link at midpoint
        ########################################
        if (np.random.random() < self.s) and len(self.adjacency) > 1:
            # choose link at random:
            elist = sorted(set([self._s(key) for key in self.adjacency.iterkeys()]))

            ei = np.random.randint(len(elist))
            e = elist[ei]
            a, b = e[0], e[1]


            # add node at midpoint and calc distances:
            self.lat[node] = (self.lat[a] + self.lat[b]) / 2.
            self.lon[node] = (self.lon[a] + self.lon[b]) / 2.
            if not self.low_memory:
                self._update_distance()

            self.adjacency[a, b] = self.adjacency[b, a] = 0
            self.adjacency[a, node] = self.adjacency[node, a] = 1
            self.adjacency[b, node] = self.adjacency[node, b] = 1

            if self.r > 0:
                self.dGmatrix[:(node + 1), :(node + 1)] = self._get_graph_distances()

            self.added_edges += 1

            if self.debug:
                print "G5", (int(a), int(b))

            #TODO: make shure (a, node) and (b, node) are not selected again?
        else:
            # step G2: link to nearest
            ##########################
            if node == 1:
                target = 0
            else:
                target = self._get_closest_connected_node(node, node - 1)
            self.adjacency[node, target] = self.adjacency[target, node] = 1


            if self.r > 0:
                # adjust network distances:
                #self.dGmatrix[:(node + 1), :(node + 1)] = self._get_graph_distances()
                self.dGmatrix[node, :self.added_nodes] = self.dGmatrix[target, :self.added_nodes] + 1
                self.dGmatrix[:self.added_nodes, node] = self.dGmatrix[:self.added_nodes, target] + 1
                self.dGmatrix[node, node] = 0

            self.added_edges += 1

            if self.debug:
                print "G2", (node, target)

            # step G3: add optimal redundant link to node
            #############################################
            if np.random.random() < self.p:
                candidates = {}
                for v in xrange(node - 1):
                    if self.adjacency[v, node] == 0:
                        candidates[(v, node)] = self._get_distance(v, node) if self.low_memory else self.distance[(v, node)]

                # there might be no candidates if n0 = 1
                if len(candidates) > 0:
                    if self.r > 0:
                        #dGmatrix = self._get_graph_distances()

                        for (u, v) in candidates.iterkeys():
                            candidates[(u, v)] /=  ( 1. + self.dGmatrix[u, v])**self.r

                    a, b = min(candidates, key=candidates.get)

                    self.adjacency[a, b] = self.adjacency[b, a] = 1

                    if self.r > 0:
                        if a == node:
                            target = b
                        else:
                            target = a
                        # adjust network distances:
                        self.dGmatrix = np.minimum(
                            np.minimum(self.dGmatrix,
                                       self.dGmatrix[:, [node]] +
                                       np.ones([self.n, self.n]) +
                                       self.dGmatrix[[target],:]
                            ),
                            self.dGmatrix[:,[target]] +
                            np.ones([self.n, self.n]) +
                            self.dGmatrix[[node], :]
                        )

                    self.added_edges += 1

                    if self.debug:
                        print "G3", (a, b)


            # step G4: add another optimal redundant link to random node
            ############################################################
            if np.random.random() < self.q:
                i2 = np.random.randint(node)

                candidates = {}
                for v in xrange(node):
                    if v == i2:
                        continue
                    if self.adjacency[v, i2] == 0:
                        candidates[self._s((v, i2))] = self._get_distance(v, i2) if self.low_memory else self.distance[self._s((v, i2))]


                # there might be no candidates if n0 = 1
                if len(candidates) > 0:
                    if self.r > 0:
                        #dGmatrix = self._get_graph_distances()

                        for (u, v) in candidates.iterkeys():
                            candidates[(u, v)] /= ( 1. + self.dGmatrix[u, v])**self.r

                    a, b = min(candidates, key=candidates.get)
                    self.adjacency[a, b] = self.adjacency[b, a] = 1


                    if self.r > 0:
                        if a == i2:
                            target = b
                        else:
                            target = a
                        # adjust network distances:
                        self.dGmatrix = np.minimum(
                            np.minimum(self.dGmatrix,
                                       self.dGmatrix[:, [i2]] +
                                       np.ones([self.n, self.n]) +
                                       self.dGmatrix[[target], :]
                            ),
                            self.dGmatrix[:, [target]] +
                            np.ones([self.n, self.n]) +
                            self.dGmatrix[[i2], :]
                        )

                    self.added_edges += 1

                    if self.debug:
                        print "G4", i2, (a, b)

        if self.debug and self.r > 0:
            # check that update went well
            assert ((self.dGmatrix[:(node + 1), :(node + 1)] - self._get_graph_distances())**2).sum() == 0


    def _validation(self, attr, value):
        if attr == "n0" or attr == "n":
            if value < 1:
                return False
            else:
                return True
        elif attr == "r":
            if value < 0:
                return False
            else:
                return True
        else:
            if value < 0 or value > 1:
                return False
            else:
                return True

    def _initial_mst(self):
        adjacency = np.zeros([self.n0, self.n0])
        np.fill_diagonal(adjacency, 0)

        self.mst_edges = self._get_mst()
        for edge in self.mst_edges:
            adjacency[edge[0], edge[1]] = adjacency[edge[1], edge[0]] = 1
            self.added_edges += 1

        self.adjacency = dok_matrix(adjacency)


    def _get_mst(self):
        full_graph = Graph.Full(self.n0)
        factor = 1e5 # since small weights lead to MST problem
        weights = [factor * self.distance[self._s((edge.source,edge.target))] for edge in full_graph.es]
        G = full_graph.spanning_tree(weights).as_undirected()
        return G.get_edgelist()


    def _get_graph_distances(self):
        elist = sorted(set([self._s(key) for key in self.adjacency.iterkeys()]))
        G = Graph(self.added_nodes)
        G.add_edges(elist)
        return np.array(G.shortest_paths())


    def _add_random_locations(self, _m):
        m = int(_m)
        if m < 1:
            print "ERROR: You have to add a positive integer number of nodes."
        else:
            for i in xrange(m):
                pos = self._get_coords()
                self.lat.append(pos[0])
                self.lon.append(pos[1])
            self._update_distance()

    def _get_closest_connected_node(self, source, connected):
        # vertices need to be properly ordered for this to work, i.e. nodes in the connected component
        # should be labeled from 0 to connected-1
        min = np.inf
        target = source
        for node in xrange(connected):
            if source == node:
                # should actually not happen
                continue
            elif self.adjacency[source, node] == 0:
                d = self._get_distance(node, source) if self.low_memory else self.distance[(node, source)]
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

class RPG(RpgAlgorithm):
    def __init__(self):

        super(RPG, self).__init__()


        # where am I?
        self.basepath = os.getcwd()
        self.figdir = os.path.join(self.basepath, "figures/")
        self.netdir = os.path.join(self.basepath, "networks/")
        if not os.path.exists(self.figdir):
            os.makedirs(self.figdir)
        if not os.path.exists(self.netdir):
            os.makedirs(self.netdir)

        i = 0
        self.identifier = "random_net" + str(i)
        while os.path.exists(self.netdir + self.identifier + ".gml"):
            i += 1
            self.identifier = "random_net" + str(i)



    ###############################################################################
    # ##                       PUBLIC FUNCTIONS                                ## #
    ###############################################################################

    def save_graph(self, info_file=True):
        elist = sorted(set([self._s(key) for key in self.adjacency.iterkeys()]))

        G = Graph(self.added_nodes)
        G.add_edges(elist)
        G.vs['lat'] = self.lat
        G.vs['lon'] = self.lon

        G.write_graphml(self.netdir + self.identifier + ".gml")

        if info_file:
            if not hasattr(self, "_stat"):
                self._stat = self.stats
            with open(self.netdir + self.identifier + ".gml.info", "w") as f:
                f.write(self._stat + "\n")

        return self.netdir + self.identifier + ".gml"


    @property
    def stats(self):
        from scipy.linalg import eigvals

        elist = sorted(set([self._s(key) for key in self.adjacency.iterkeys()]))

        G = Graph(self.added_nodes)
        G.add_edges(elist)

        STR = str()
        STR += "connected:" + str(G.is_connected()) + "\n"
        STR += "undirected:" + str(not G.is_directed()) + "\n"
        STR += "\n"
        STR += "mean degree kbar:" + str(np.mean(G.degree())) + "\n"
        STR += "ratio r{k>kbar}:" + str(sum(G.degree() > np.mean(G.degree())) * 1. / self.added_nodes) + "\n"
        STR += "avg. neighbour's degree distribution:" + str(np.min(G.knn()[0])) + \
            "..." + "<k>=" + str(np.mean(G.knn()[0])) + \
            "..." + str(np.max(G.knn()[0])) + "\n"
        STR += "degree - degree frequency - avg. neighbour's degree:"
        fd, _ = np.histogram(G.degree(), bins=np.max(G.degree()))
        for i, val in enumerate(G.knn()[1]):
            STR += "    " + str(i + 1) + "-" + str(fd[i]) + "-" + str(val) + "\n"
        STR += "s.p. betweenness distribution:" + str(np.min(G.betweenness())) + \
            "..." + "<b>=" + str(np.mean(G.betweenness())) + \
            "..." + str(np.max(G.betweenness())) + "\n"
        STR += "\n"
        STR += "average shortest path length:" + str(G.average_path_length()) + "\n"
        STR += "network transistivity:" + str(G.transitivity_undirected()) + "\n"
        STR += "degree assortativity:" + str(G.assortativity_degree()) + "\n"
        STR += "Fiedler eigenvalue lambda2:" + str(sorted(eigvals(G.laplacian()))[1]) + "\n"
        STR += "\n"
        STR += "number of dead ends:" + str(G.degree().count(1)) + "\n"
        n_a_dt = G.betweenness().count(self.added_nodes - 2) + \
            G.betweenness().count(2 * self.added_nodes - 5) + \
            G.betweenness().count(2 * self.added_nodes - 6) + \
            G.betweenness().count(3 * self.added_nodes - 10) + \
            G.betweenness().count(4 * self.added_nodes - 17) + \
            G.betweenness().count(5 * self.added_nodes - 26)
        STR += "number of nodes adjacent to dead trees (approx.):" + str(n_a_dt) + "\n"
        STR += "number of detour nodes:" + str(G.betweenness().count(0) - G.degree().count(1)) + "\n"
        self._stat = STR

        return STR


    def plot_net(self, name="random_network", labels=False, crop=False):

        elist = sorted(set([self._s(key) for key in self.adjacency.iterkeys()]))

        filename = self.figdir + name + "_" + self.identifier + ".pdf"

        x = (np.max(self.lat) - np.min(self.lat))
        y = np.cos(np.mean(self.lat)) * (np.max(self.lon) - np.min(self.lon))
        #print x, y

        G = Graph(self.added_nodes)
        G.add_edges(elist)

        visual_style = {}

        def edgecolor(i,j):
            if hasattr(self, "mst_edges"):
                if self._s((i,j)) in self.mst_edges:
                    return "red"
            else:
                return "black"

        scale = 2
        visual_style["layout"] = zip(self.lat, self.lon)
        visual_style["bbox"] = (x / y * 1024, 1024)
        visual_style["margin"] = 10 * scale
        visual_style["palette"] = palettes["heat"]



        visual_style["edge_color"] = [edgecolor(edge.source,edge.target) for edge in G.es]
        visual_style["edge_width"] = [2 * scale for edge in G.es]

        visual_style["vertex_size"] = [10 * scale for i in range(G.vcount())]
        visual_style["vertex_color"] = [int(x) for x in rescale(1. - self.bs_predictor(),
                                                                out_range=(0, len(visual_style["palette"]) - 1))]
        if labels:
            visual_style["vertex_label"] = [str(i) for i in range(G.vcount())]

        plot(G, filename, **visual_style)

    def bs_predictor(self):
        try:
            from pyunicorn.core import resistive_network as rn
        except:
            print "ERROR: pyunicorn not availiable."
            return np.zeros(self.added_nodes)

        def min_clust(W, N):
            ''' W is the admittance matrix'''
            C = np.zeros(N)
            for i in xrange(N):
                norm = 0
                for j in xrange(N):
                    for k in xrange(j):
                        if j!=k:
                            norm += W[i,k]*W[i,j]
                            if W[i,j]*W[i,k]*W[j,k]>0:
                                C[i] += min( W[i,k],W[j,k] ) * min( W[i,j],W[j,k] )
                                #print i,norm[i],C[i]
                C[i] /=  max(norm,1)
                #print C[i]
            return C

        weights = self.adjacency
        for key in weights.keys():
            weights[key] = self.distance[self._s(key)]

        net = rn.ResNetwork(resistances=weights.todense())

        ''' explanatory variables '''
        ED = net.admittive_degree()
        #print 'ED'
        ANED = np.divide(net.average_neighbors_admittive_degree(), ED)
        #print 'ANED'
        minC = min_clust(net.get_admittance(), net.N)
        #print 'minC'
        VCFB = np.zeros(net.N)
        for a in xrange(net.N):
            VCFB[a] = 1.0*net.vertex_current_flow_betweenness(a) *((net.N*(net.N-1)) / 2) /net.N
        #print 'VCFB'
        ERCC = np.zeros(net.N)
        for a in xrange(net.N):
            ERCC[a] = net.effective_resistance_closeness_centrality(a)
        #print 'ERCC'
        dead = np.zeros(net.N)
        for i, b in enumerate(net.betweenness()):
            if b in [net.N - 2, 2 * net.N - 5, 2 * net.N - 6, 3 * net.N - 10, 4 * net.N - 17, 5 * net.N - 26]:
                dead[i] = 1
            else:
                dead[i] = 0
        #print 'dead'

        '''predicted probability to have poor basin stability'''
        # %Coefficients:
        # %              Estimate Std. Error z value Pr(>|z|)
        # %(Intercept)   3.325612   0.143404   23.19   <2e-16 ***
        # %ED            0.088743   0.003832   23.16   <2e-16 ***
        # %ANED         -0.348762   0.011514  -30.29   <2e-16 ***
        # %minC        -10.389121   0.271047  -38.33   <2e-16 ***
        # %VCFB         -0.107782   0.007994  -13.48   <2e-16 ***
        # %ERCC         -1.515226   0.033785  -44.85   <2e-16 ***
        # %dead          4.925139   0.084272   58.44   <2e-16 ***


        g = 3.325612 + 0.088743*ED -0.348762*ANED -10.389121*minC -0.107782*VCFB -1.515226*ERCC + 4.925139*dead
        prob = 1.0/(1.0+np.exp(-g))
        t = 0.15
        poor_bs = np.zeros(net.N)

        for i in xrange(net.N):
            if prob[i]>t:
                poor_bs[i]=1
        return prob



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

    print g.stats


if __name__ == "__main__":
    main()






