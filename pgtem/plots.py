import awesomeplot.core as plt
import igraph as ig
import numpy as np

canvas = plt.AwesomePlot(output="paper")
g = ig.read("test.network", format="picklez")

adj = np.array(g.get_adjacency().data)

c = np.array(map(lambda x: canvas.pik_colours.colors[0] if x == "rpg_edge" else canvas.pik_colours.colors[1],
               g.es["type"]))

canvas.add_network(adj,
                   height=True,
                   sym=False,
                   styles={"layout": np.c_[g.vs["lat"], g.vs["lon"], g.vs["time"]], "edge_color": c},
                   axis_labels=["lat", "lon", "time"])

canvas.add_network(adj,
                   height=False,
                   sym=False,
                   labels=True,
                   styles={"layout": np.c_[g.vs["lat"], g.vs["lon"]], "edge_color": c},
                   axis_labels=["lat", "lon", "time"])

canvas.show()