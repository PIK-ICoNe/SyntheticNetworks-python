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

        return os.path.join(self.netdir, self.identifier + ".gml")


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
        #visual_style["palette"] = palettes["heat"]

        visual_style["edge_color"] = [edgecolor(edge.source,edge.target) for edge in G.es]
        visual_style["edge_width"] = [2 * scale for edge in G.es]

        visual_style["vertex_size"] = [10 * scale for i in range(G.vcount())]
        visual_style["vertex_color"] = "black"
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

    g.plot_net()


if __name__ == "__main__":
    main()






