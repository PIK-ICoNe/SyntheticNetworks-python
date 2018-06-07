
# coding: utf-8

from rpgm_algo import RpgAlgorithm
import numpy as np
import awesomeplot.core as ap

# In[1]:

class NeoNetPlot(ap.Plot):
    def add_network(self, adjacency, styles={}, height=False):
        from matplotlib import pyplot
        import seaborn
        
        if height:
            from mpl_toolkits.mplot3d import Axes3D
        from scipy.sparse import issparse, isspmatrix_dok

        if issparse(adjacency):
            assert isspmatrix_dok(adjacency)
            # print "Build network from sparse dok matrix."
            N = adjacency.shape[0]
            edgelist = sorted(set([tuple(np.sort(key)) for key in adjacency.iterkeys()]))
        else:
            N = len(adjacency)
            edgelist = np.vstack(np.where(adjacency > 0)).transpose()
            edgelist = sorted(set([tuple(np.sort(edgelist[i])) for i in range(len(edgelist))]))

        cmap = self.dfcmp

        visual_style = dict(
            edge_color=np.repeat('#8e908f', len(edgelist)),
            edge_width=seaborn.axes_style()["axes.linewidth"],
            vertex_size=100,
            vertex_label=range(N)
        )

        if styles:
            visual_style.update(styles)

        if visual_style.has_key("edge_color_dict"):
            f = lambda x: visual_style["edge_color_dict"][x]
            for i, e in enumerate(edgelist):
                visual_style["edge_color"][i] = f(e)

        if height:
            fig = pyplot.figure()
            ax = fig.gca(projection='3d')
            x, y, z = zip(*visual_style["layout"])
            args = (x, y, z)
        else:
            fig, ax = pyplot.subplots(nrows=1, ncols=1)
            fig.tight_layout()
            x, y = zip(*visual_style["layout"])
            args = (x, y)

        ax.set_axis_off()

        for i, e in enumerate(edgelist):
            edge = np.vstack((visual_style["layout"][e[0]], visual_style["layout"][e[1]]))
            if height:
                xyz = edge[:, 0], edge[:, 1], edge[:, 2]
            else:
                xyz = edge[:, 0], edge[:, 1]
            ax.plot(*xyz,
                    color=visual_style["edge_color"][i],
                    linestyle='-',
                    lw=visual_style["edge_width"],
                    alpha=0.4,
                    zorder=1)

        margin = max(0.05 * (np.max(x) - np.min(x)), 0.05 * (np.max(y) - np.min(y)))
        ax.set_xlim([np.min(x) - margin, np.max(x) + margin])
        ax.set_ylim([np.min(y) - margin, np.max(y) + margin])

        if not visual_style.has_key("vertex_color"):
            nodes = ax.scatter(*args,
                               c='#8e908f',
                               s=visual_style["vertex_size"],
                               cmap=cmap,
                               edgecolor='w',
                               zorder=2, 
                               alpha=1)
        else:
            nodes = ax.scatter(*args,
                               c=visual_style["vertex_color"],
                               s=visual_style["vertex_size"],
                               cmap=cmap,
                               vmin=np.floor(np.min(visual_style["vertex_color"])),
                               vmax=np.ceil(np.max(visual_style["vertex_color"])),
                               edgecolor='w',
                               zorder=2, 
                               alpha=1)
            #cb = fig.colorbar(nodes, orientation='horizontal', shrink=0.66, format=r"%.1f")

        # we may adjust the background colour to make light nodes more visible
        #ax.set_axis_bgcolor((.9, .9, .9))
        
        fig.tight_layout()
        
        self.figures.append(fig)

        return fig, ax
    
    def add_hist(self, data, label='x', nbins=None, log=False, legend=True):
        from matplotlib import pyplot
        import seaborn

        if not isinstance(data, dict):
            if hasattr(data, "__iter__"):
                data = {0: data}
                legend = False
            else:
                ValueError("Input data should be dict or list/array.")

        x_range = np.unique(np.concatenate(data.values()))
        xmin = np.min(x_range)
        xmax = np.max(x_range)

        fig, ax = pyplot.subplots(nrows=1, ncols=1)


        fig.tight_layout()

        ymax = 0.
        counter = 0

        if nbins is None:
            incr = np.diff(x_range)
            bin_edges = np.unique(np.append(x_range[:-1] - incr * .5, x_range[1:] + incr * .5))
            width = np.min(np.diff(bin_edges))
        else:
            width = 1. * (xmax - xmin) / nbins
            bin_edges = np.linspace(xmin - width, xmax + width, nbins + 1, endpoint=True)

        bottom = np.zeros(len(bin_edges) - 1)

        for i, k in enumerate(data.keys()):
            h, b = np.histogram(data[k], bins=bin_edges, density=True)
            ax.bar(.5 * (b[1:] + b[:-1]), h, bottom=bottom, color=self.dfcmp.colors[i], edgecolor="none", align='center', zorder=1, width=width, label=k,
                   alpha=1, log=log)
            bottom += h
            counter += 1
            ymax += h.max()

        ax.set_xlim([bin_edges[0], bin_edges[-1]])
        ax.set_ylim([0., ymax * 1.1])
        ax.set_xlabel(label)
        ax.set_ylabel(r'density')

        if legend:
            pyplot.legend()

        ax.yaxis.grid(color='w', linestyle='-', zorder=2)

        self.figures.append(fig)

        return fig




# In[ ]:

class Parameters(object):
    def __init__(self, **kwargs):
        self.centre = [0, 0]
        self.boundaries = [1, 1]
        self.n = int(10)
        self.n0 = int(1)
        self.p = .2
        self.q = .3
        self.r = 1./3.
        self.s = .1
        
        for k, v in kwargs.iteritems():
            if hasattr(self, k):
                setattr(self, k, v)
                
        self.m = 1
        self.grow = True
    
    @classmethod
    def mv(cls, centre, n):
        return cls(centre=centre, boundaries=[.3, .3], q=0, p=0.1, n=n, n0=int(n-1), grow=False)
    
    @classmethod
    def lv(cls, centre, n):
        return cls(centre=centre, boundaries=[.1, .1], q=0, p=0, n=n, n0=int(n-1), grow=False)
    
    


# In[ ]:

def shift(l, m):
    if isinstance(l, tuple):
        return (l[0]+m, l[1]+m)
    else:
        return [(v[0]+m, v[1]+m) for v in l]


# In[ ]:

def centered_coords(m, centre, boundaries):
    coords = np.zeros([m, 2])
    coords[0] = centre
    for i in xrange(1, m):
        coords[i] = (.5 - np.random.uniform(size=2)) * np.array(boundaries) + np.array(centre)
    return coords


# In[ ]:

def network(pars, loc):
    g = RpgAlgorithm()
    g.debug = False
    g.set_params(n=pars.n, n0=pars.n0, r=pars.r, p=pars.p, q=pars.q, s=pars.s)

    # first node will be connected via a transformer
    g.setup_locations(sampling="data", locations=loc)

    g.initialise()
    
    if pars.grow:
        g.grow()
        
    g.added_transformers = 0
    
    return g


# In[ ]:

def append_net(a, b, connections):
    assert isinstance(a, RpgAlgorithm)
    assert isinstance(b, RpgAlgorithm)
    a.adjacency.resize([a.n + b.n, a.n + b.n])
    for e in b.edgelist():
        (s, t) = shift(e, a.n)
        a.adjacency[s, t] = a.adjacency[t, s] = 1
    # add connections between a and b
    for edge in connections:
        (s, t) = (edge[0], edge[1] + a.n)
        a.adjacency[s, t] = a.adjacency[t, s] = 1
    a.added_transformers += len(connections)
    # update locations
    for prop in ["lat", "lon", "locations"]:
        xa = np.array(getattr(a, prop))
        xb = np.array(getattr(b, prop))
        setattr(a, prop, np.concatenate([xa, xb]))
    # update distance
    a.distance.update({shift(key, a.n): value for (key, value) in  b.distance.items()})
    a.n += b.n
    a.added_nodes += b.added_nodes
    a.added_edges += b.added_edges
    #return a
    


# In[ ]:

def independent_layers(n_high, n_middle, n_low, high_to_middle, middle_to_low, debug=0, seed=42):
    np.random.seed(seed)

    pars = Parameters(n=n_high, n0=int(.8 * n_high))
    loc = centered_coords(pars.n, pars.centre, pars.boundaries)
    g = network(pars, loc)


    if debug:
        print "HV"
        print g

    # for bookkeeping ...
    v_type = 3 * np.ones(g.n)

    # select transformers to mv level
    transformers = np.random.choice(range(g.n), size=int(high_to_middle * g.n), replace=False)

    # attach mv grid to each transformer
    for i, t in enumerate(transformers):

        pars = Parameters.mv(centre=g.locations[t], n=n_middle)
        loc = centered_coords(pars.n, pars.centre, pars.boundaries)
        g_mv = network(pars, loc)

        if debug:
            print "MV @", t, g.locations[t]
            print g_mv

        v_type = np.append(v_type, 2 * np.ones(g_mv.n))

        # list all nodes that are not connected to the hv layer (i.e. all beside the first one ...)
        mv_nodes = range(g.n + 1, g.n + g_mv.n)

        append_net(g, g_mv, [(t, 0)])

        # select transformers to lv layer, attach lv grid
        for ii, tt in enumerate(np.random.choice(mv_nodes, int(middle_to_low * g_mv.n))):

            pars = Parameters.lv(centre=g.locations[tt], n=n_low)
            loc = centered_coords(pars.n, pars.centre, pars.boundaries)
            g_lv = network(pars, loc)

            if debug:
                print "LV @", tt, g.locations[tt]
                print g_lv

            v_type = np.append(v_type, 1 * np.ones(g_lv.n))

            append_net(g, g_lv,  [(tt, 0)])


    e_type = {}
    for edge in g.edgelist():
        if v_type[edge[0]] == v_type[edge[1]]:
            e_type[edge] = v_type[edge[0]]
        else:
            e_type[edge] = 0
            
    return g, v_type, e_type

def main():
    names = ["high", "middle", "low"]
    degrees = {names[i]: [] for i in range(3)}
    link_length = {names[i]: [] for i in range(3)}

    for sample in range(2):
        g, v_type, e_type = independent_layers(10, 15, 12, 6. / 10., 14. / 15., seed=sample)
        print sample,

        for layer in range(3):
            degrees[names[layer]].extend(
                [g.adjacency[node, :].sum() for node in np.where(v_type == (layer + 1))[0].tolist()]
            )

        for layer in range(3):
            link_length[names[layer]].extend(
                np.log([g.distance[k] for (k, v) in e_type.items() if v == (layer + 1)])
            )

    canvas = NeoNetPlot.paper(font_scale=3)
    fig2 = canvas.add_hist(link_length, label="k", nbins=20)

    canvas.show()

if __name__ == '__main__':
    main()
