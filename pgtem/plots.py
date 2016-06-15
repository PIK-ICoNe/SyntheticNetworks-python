import awesomeplot.core as plt
import igraph as ig
import numpy as np

canvas = plt.AwesomePlot(output="paper")
g = ig.read("test.network", format="picklez")

adj = np.array(g.get_adjacency().data)

canvas.add_network(adj, height=True, styles={"layout": np.c_[g.vs["lat"], g.vs["lon"], g.vs["time"]]},
                   axis_labels=["lat", "lon", "time"])

canvas.show()