import awesomeplot.core as plt
import pandas as pd
from scipy.sparse import dok_matrix
import numpy as np

canvas = plt.AwesomePlot(output="paper")

nodes = pd.read_pickle("test.network.nodes")
edges = pd.read_pickle("test.network.edges")

N = len(nodes)
m = len(edges)

adj = dok_matrix(np.zeros([N, N]))

for i in xrange(m):
    adj[int(edges.source[i]), int(edges.target[i])] = \
        adj[int(edges.target[i]), int(edges.source[i])] = edges["length"].values[i]


c = {}
for i, t in enumerate(edges["type"].values):
    c[(int(edges["source"][i]), int(edges["target"][i]))] = "red" if t == "rpg_edge" else "black"



canvas.add_network(adj,
                   height=True,
                   sym=False,
                   styles={"layout":  nodes[["lat", "lon", "time"]].values,
                           "edge_color_dict": c},
                   axis_labels=["lat", "lon", "time"])

canvas.add_network(adj,
                   height=False,
                   sym=False,
                   labels=True,
                   styles={"layout":  nodes[["lat", "lon"]].values,
                           "edge_color_dict": c},
                   axis_labels=["lat", "lon", "time"])

canvas.show()