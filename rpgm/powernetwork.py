#!/usr/bin/python
# -*- coding: utf-8 -*-

__version__ = "1.0"
__author__ = "Paul Schultz, pschultz@pik-potsdam.de"
__copyright__ = "Copyright (C) 2016 Paul Schultz, GNU GPL 3"


"""
Module contains class PowerNetwork.

Provides a data structure power grid topologies with additional data
assigned to vertices and edges. It is subclassed from igraph.Graph
and inherits all methods from Graph objects.

The class has the no instance variables.

Overriden inherited methods::

    (str)       __str__         : extended description
    (NoneType)  save            : include file identifier and fixed file format to graphml
"""


# Import NumPy for the array object and fast numerics.
import numpy as np
# Import os to set output file paths.
import os
# Import time to assign timestamps to PowerNetwork objects.
import time
# Import the class we inherit from.
from igraph import Graph

# TODO: future interoperability with PyPSA, problem with v-attribute name already used in igraph
# Import component specifications from PyPSA
import grid_components as com

# use pint to handle physical units
# Conventions:
# - W for real power
# - var for reactive power
# - VA for apparent power

from pint import UnitRegistry
ureg = UnitRegistry()
ureg.define("var = W")


class PowerNetwork(Graph):
    """
    PowerNetwork class

    The class PowerNetwork builds upon an igraph.Graph to provide a convenient way to
    handle power grid topologies with a fixed set of vertex and edge properties needed
    for visualisation and simulations.

    **Examples:**

    >>> print Graph.Lattice([3, 3], 1)
    IGRAPH U--- 9 18 --
    + edges:
    0 -- 1 2 3 6   2 -- 0 1 5 8   4 -- 1 3 5 7   6 -- 0 3 7 8   8 -- 2 5 6 7
    1 -- 0 2 4 7   3 -- 0 4 5 6   5 -- 2 3 4 8   7 -- 1 4 6 8

    """

    # set default output directory for figures / graphml
    loc = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # position of the script
    output_dir = os.path.join(loc, "figures/")

    # implemented components from grid_components.py
    implemented_components = ["generator", "load", "prosumer", "passive", "line", "transformer"]

    # TODO: ensure per unit system.

    ###############################################################################
    # ##                       MAGIC FUNCTIONS                                 ## #
    ###############################################################################

    def __init__(self, name="", n=0, edgelist=None, v_types=None, e_types=None):
        """
        Initialise an instance of PowerNetwork. This class is derived from an igraph.Graph class.

        Parameters
        ----------
        name: string
            label used to identify a PowerNetwork object
        n: int, optional
            number of vertices
        edgelist: list, optional
            list of pairs '(source, target)' defining edges
        v_types: dict, optional
            dictionary containing attribute values for each vertex
            allowed attributes are set in Components:
            generator, passive, load
        e_dict : dict, optional
            dictionary containing attribute values for each edge
            allowed attributes are set in Components:
            line,transformer
        """

        if n is 0:
            if type(edgelist) == list and len(edgelist) > 0:
                n = np.max(edgelist) + 1
        else:
            assert type(n) == int and n > 0
            if type(edgelist) == list and len(edgelist) > 0:
                assert n == np.max(edgelist) + 1

        if type(edgelist) is not list:
            print type(edgelist)
            raise ValueError("Invalid type of edgelist!")

        # system base in MVA
        self.pbase = 3 * ureg("MW")
        self.vbase = 20 * ureg("kV")
        # system reference frequency in Hertz
        self.rfreq = 50 * ureg("Hz")

        super(self.__class__, self).__init__(
            n=n,
            edges=edgelist,
            directed=False,
            graph_attrs={"Identifier": name,
                         "Timestamp": time.asctime(time.localtime(time.time())),
                         "Frequency": self.rfreq.magnitude,
                         "PBase": self.pbase.magnitude}
        )

        # the following sets up descriptor instances for each vertex and link
        # TODO: possibility to add external data

        if v_types and e_types:
            self._set_default_components(v_types, e_types)



    def __str__(self):
        """
        Output of passing a PowerNetwork instance to the print command.

        Returns
        -------
        text: string

        **Examples:**

        >>> print PowerNetwork.SmallTestPowerGrid() # doctest: +ELLIPSIS
        IGRAPH U--- 5 6 --
        + attr: Frequency (g), Identifier (g), PBase (g), Timestamp (g), Bus_ID (v),
          _T (v), lat (v), lon (v), _T (e)
        + edges:
        0--1 0--1 1--2 1--3 2--3 3--4
        + SmallTestPowerGrid
        ...
        + Average degree: 2.4

        """

        text = Graph.__str__(self)
        text += "\n+ " + self["Identifier"] + "\n+ " + self["Timestamp"]
        # TODO: include some basic statistics here
        if self.vcount() > 0:
            text += "\n+ Average degree: " + str(np.mean(self.degree()))
        return text

    ###############################################################################
    # ##                       PUBLIC FUNCTIONS                                ## #
    ###############################################################################
    def save(self, filename, use_identifier=True):
        """
        Saves a PowerNetwork instance to a GraphML file.

        Element descriptors are converted to bus and branch attributes.

        Parameters
        ----------
        filename: str
            name of output graphML appended to self["Identifier"]

        """

        _graph = self.copy()

        _graph._move_attributes(_graph)
        #print _graph.vs.attribute_names()
        _graph._separate_units(_graph)
        #print _graph.vs.attribute_names()

        identifier = ""
        if use_identifier:
            identifier = self["Identifier"] + "_"

        i = 0
        path = os.path.join(self.output_dir, identifier + filename + "_" + str(i) + ".graphml")
        while os.path.exists(path):
            i += 1
            path = os.path.join(self.output_dir, identifier + filename + "_" + str(i) + ".graphml")

        _graph.write(path, format="graphml")

        return path

    def visualise(self, filename):
        """
        Basic Network visualisation.

        Parameters
        ----------
        filename: str
            name of output PDF appended to self["Identifier"]

        """
        if self.vcount() == 0:
            print "Nothing to visualise, no vertices in the graph!"
            return

        from igraph import plot, palettes

        _graph = self.copy()
        assert isinstance(_graph, PowerNetwork)

        self._move_attributes(_graph)

        # dummy values for plotting if no further data is given
        if not hasattr(_graph.vs, "_T"):
            _graph.vs["volt"] = 220
            _graph.es["volt"] = 220
            _graph.vs["Type"] = None
            _graph.es["Type"] = None

        levels = np.unique(_graph.vs["volt"])

        visual_style = dict(
            margin=30,
            mark_groups=None if len(levels) == 1 \
                else {_graph.vs.select(volt=l): int(i * 256. / len(levels)) for i, l in enumerate(levels)},
            edge_color="black",
            # edge_curved=0.1, # does not work with multiple edges
            edge_width=2,
            edge_order_by="volt",  # draw low-voltage grids first
            layout=zip(_graph.vs["lat"], _graph.vs["lon"]),
            palette=palettes["heat"],
            vertex_size=30,
            vertex_shape="circle",
            vertex_color="gray",
            vertex_frame_color="white",
            vertex_frame_width=3,
            vertex_label=_graph.vs["Bus_ID"],
            vertex_label_dist=0,
            vertex_label_color="black",
            vertex_label_size=12,
            vertex_label_angle=0,
            vertex_order_by="volt"  # draw low-voltage grids first
        )

        def vertex_shape(val):
            if val == "generator":
                return "circle"
            elif val == "load":
                return "square"
            elif val == "prosumer":
                return "circle"
            elif val == "passive":
                return "diamond"
            else:
                return "circle"

        def edge_width(val):
            if val == "line":
                return 2
            elif val == "transformer":
                return 6
            else:
                return 2

        visual_style["vertex_shape"] = map(vertex_shape, _graph.vs["Type"])
        visual_style["edge_width"] = map(edge_width, _graph.es["Type"])

        try:
            visual_style["vertex_label"] = _graph.vs["Bus_ID"]
        except:
            visual_style["vertex_label"] = range(self.vcount())

        if None in _graph.vs["lat"] or None in _graph.vs["lon"]:
            visual_style["layout"] = _graph.layout_auto()
            print "Assign random layout for plotting."

        plot(_graph, self.output_dir + self["Identifier"] + "_" + filename + ".pdf", **visual_style)

    def vertex_filter(self, att, val, op, return_subgraph=False):
        """
        Filters a PowerNetwork graph by a given vertex attribute.

        Parameters
        ----------
        att: string
            a valid vertex attribute name
        val: int/float/string
            value of vertex attribute used for filtering
        op: string
            filter operation, one of (_eq, _ne, _lt, _gt, _le, _ge, _in, _notin)
        return_subgraph: bool, optional
            if true, a new Graph object containing the filtered subgraph is returned,
            else, vertices not fulfilling the op condition are deleted (default)

        Returns
        -------
        subgraph: PowerNetwork object
            returned only if return_subgraph=True

        """

        # TODO: check what happens to edges adjecent to deleted vertices

        assert op in ["_eq", "_ne", "_lt", "_gt", "_le", "_ge", "_in", "_notin"]
        assert att in self.vs.attribute_names()
        if return_subgraph:
            subgraph = self.vs.select(**{att + op: val}).subgraph()
            return subgraph
        else:
            self.vs.select(**{att + op: val}).delete()
            pass

    def edge_filter(self, att, val, op, return_subgraph=False):
        """
        Filters a PowerNetwork graph by a given edge attribute.

        Parameters
        ----------
        att: string
            a valid vertex attribute name
            ignored if op is one of (_source, _target, _within, _between)
        val: int/float/string
            value of edge attribute used for filtering
        op: string
            filter operation, one of (_eq, _ne, _lt, _gt, _le, _ge, _in, _notin, _source, _target, _within, _between)
        return_subgraph: bool, optional
            if true, a new Graph object containing the filtered subgraph is returned,
            else, edges not fulfilling the op condition are deleted (default)

        Returns
        -------
        subgraph: PowerNetwork object
            returned only if return_subgraph=True

        """

        assert op in ["_eq", "_ne", "_lt", "_gt", "_le", "_ge", "_in", "_notin", "_source", "_target", "_within",
                      "_between"]
        if op in ["_source", "_target", "_within", "_between"]:
            att = ""
        else:
            assert att in self.es.attribute_names()
        if return_subgraph:
            return self.es.select(**{att + op: val}).subgraph()
        else:
            self.es.select(**{att + op: val}).delete()
            pass

    @classmethod
    def convert_igraph_object(cls, graph, name):
        """
        Builds a PowerNetwork instance directly upon an existing igraph object.

        Parameters
        ----------
        graph: igraph.Graph
            a Graph object
        name: string
            label used to identify the created PowerNetwork onject

        """

        assert isinstance(graph, Graph)
        return cls(name=name, n=graph.vcount(), edgelist=graph.get_edgelist())

    @classmethod
    def RandomPowerGrid(cls, number_of_nodes, **params):
        """
        Creates a spatially embedded synthetic power grid topology.

        Reference:
            A Random Growth Model for Power Grids and Other Spatially Embedded Infrastructure Networks
            Paul Schultz, Jobst Heitzig, and Juergen Kurths
            Eur. Phys. J. Special Topics on "Resilient power grids and extreme events" (2014)
            DOI: 10.1140/epjst/e2014-02279-6

        Parameters
        ----------
        number_of_nodes: int
            desired number of vertices
        params: kwargs
            custom parameters overwriting the default values of the RPG algorithm

        """

        from rpgm_core import RpgAlgorithm

        # crete random power grid
        rpg = RpgAlgorithm()
        assert isinstance(rpg, RpgAlgorithm)
        rpg.set_params(n=number_of_nodes, n0=max(1, int(number_of_nodes / 10.)), r=1. / 3.)
        rpg.set_params(**params)  # overwrites default parameter values
        rpg.initialise()
        rpg.grow()
        edgelist = sorted(set([rpg._s(key) for key in rpg.adjacency.iterkeys()]))

        graph = cls(name="rpg (p,q,r,s,n0) = " + str(rpg.p) + "," + str(rpg.q) + "," + \
                             str(rpg.r) + "," + str(rpg.s) + "," + str(rpg.n0),
                    n=number_of_nodes,
                    edgelist=edgelist
                    )

        # take care of attributes
        graph.vs["Bus_ID"] = range(number_of_nodes)
        graph.vs["lat"] = rpg.lat
        graph.vs["lon"] = rpg.lon
        graph.vs["volt"] = np.ones(number_of_nodes) * 110. * ureg("kV")
        graph.es["Branch_ID"] = range(len(edgelist))
        graph.es["source"] = zip(*edgelist)[0]
        graph.es["target"] = zip(*edgelist)[1]

        return graph

    @classmethod
    def SmallTestPowerGrid(cls):
        """
        Small network for testing purposes.

        Creates a PowerNetwork instance with a small test network comprising
        basic properties of the following topology.

            0------1--------3------4
                    \      /
                     \    /
                      \  /
                       \/
                       2


        **Examples:**

        >>> stpg = PowerNetwork.SmallTestPowerGrid()
        >>> isinstance(stpg, PowerNetwork)
        True
        >>> print stpg.attributes()
        ['Timestamp', 'PBase', 'Frequency', 'Identifier']

        """

        edgelist = [(0, 1), (0, 1), (1, 2), (1, 3), (2, 3), (3, 4)]
        n = 5

        graph = cls(name="SmallTestPowerGrid", n=n, edgelist=edgelist)

        # random bus locations
        graph.vs["Bus_ID"] = range(5)
        graph.vs["lat"] = np.random.random(5)
        graph.vs["lon"] = np.random.random(5)

        graph._set_default_components(["generator", "load", "generator", "passive", "generator"],
                                      ["line", "line", "line", "line", "line", "transformer"])
        for j in [3, 4]:
            graph.vs[j]["_T"].volt = 110 * ureg("kV")

        return graph

    def add_buses(self, bus_list, v_types=None):
        """

        :param bus_list: vertex sequence of Bus_IDs
        :return:
        """

        #check for vertices already present
        new_vertices = set(bus_list) - set(self.vs['Bus_ID'])
        assert len(new_vertices) == len(bus_list)

        for i, bus in enumerate(bus_list):

            if v_types[i] == "generator":
                T = com.Generator.default(Bus_ID=bus, Type=v_types[i])
            elif v_types[i] == "load":
                T = com.Load.default(Bus_ID=bus, Type=v_types[i])
            elif v_types[i] == "passive":
                T = com.Passive.default(Bus_ID=bus, Type=v_types[i])
            else:
                raise ValueError("Not an implemented bus type!")

            self.add_vertex(bus, _T=T)


    def add_branches(self, branch_list):
        """
        :param branch_list: list of vertex tuples
        :return:
        """

        n_branches = self.ecount()

        for i, branch in enumerate(branch_list):

            if self.vs[branch[0]]['_T'].volt == self.vs[branch[1]]['_T'].volt:
                T = com.Line.default(Branch_ID=n_branches+i, Type="line", Bus1=branch[0], Bus2=branch[1])
            else:
                T = com.Transformer.default(Branch_ID=n_branches+i, Type="transformer", Bus1=branch[0], Bus2=branch[1])

            self.add_edge(branch[0], branch[1], _T=T)


    ###############################################################################
    # ##                       PRIVATE FUNCTIONS                               ## #
    ###############################################################################

    def _set_default_components(self, v_types=None, e_types=None):
        assert all(t in self.implemented_components for t in v_types)

        for i, bus in enumerate(self.vs):
            if v_types[i] == "generator":
                bus["_T"] = com.Generator.default(Bus_ID=i, Type=v_types[i])
                assert isinstance(bus["_T"], com.Generator)
            elif v_types[i] == "load":
                bus["_T"] = com.Load.default(Bus_ID=i, Type=v_types[i])
                assert isinstance(bus["_T"], com.Load)
            elif v_types[i] == "prosumer":
                bus["_T"] = com.Prosumer.default(Bus_ID=i, Type=v_types[i])
                assert isinstance(bus["_T"], com.Prosumer)
            elif v_types[i] == "passive":
                bus["_T"] = com.Passive.default(Bus_ID=i, Type=v_types[i])
                assert isinstance(bus["_T"], com.Passive)
            else:
                raise ValueError("Not an implemented bus type!")

        for j, branch in enumerate(self.es):
            f, t = branch.tuple
            if e_types[j] == "line":
                branch["_T"] = com.Line.default(Branch_ID=j, Type=e_types[j], Bus1=f, Bus2=t)
                assert isinstance(branch["_T"], com.Line)
            elif e_types[j] == "transformer":
                branch["_T"] = com.Transformer.default(Branch_ID=j, Type=e_types[j], Bus1=f, Bus2=t)
                assert isinstance(branch["_T"], com.Transformer)

    @staticmethod
    def _move_attributes(_graph):
        """ move attributes from "_T" to new dict, so that they are saved in gml-file, if "_T" exists
            :param _graph a igraph object
            :return:no return
        """
        assert isinstance(_graph, PowerNetwork)
        if "_T" in _graph.vs.attribute_names():

            for bus in _graph.vs:
                _graph._add_dict(bus, bus["_T"].__dict__)

            for branch in _graph.es:
                _graph._add_dict(branch, branch["_T"].__dict__)

            del _graph.vs["_T"]
            del _graph.es["_T"]

            # print "Bus after deleting _T", bus
            # print "Bus3 in graph.vs:" ,_graph.vs[3]
            # print "BUS_ID bus 3 in graph.vs:",_graph.vs[3]["Bus_ID"]
            # print "PG bus 3 in graph.vs:",_graph.vs[3]["PG"]

        else:
            pass

    @staticmethod
    def _separate_units(_graph):
        """ seperate the magnitude and the unit to save them in the gml file seperately
            :param _graph a igraph object
            :return:no return
        """
        assert isinstance(_graph, PowerNetwork)
        for bus in _graph.vs:

            for k in bus.attribute_names():
                if isinstance(bus[k], ureg.Quantity):
                    bus[k + "_u"] = str(bus[k].units)
                    bus[k] = bus[k].magnitude

        for branch in _graph.es:
            for k in branch.attribute_names():
                if isinstance(branch[k], ureg.Quantity):
                    # print "seperate Units for", branch[k]
                    branch[k + "_u"] = str(branch[k].units)
                    branch[k] = branch[k].magnitude


    @staticmethod
    def _add_dict(obj, d):
        """
        Adding attributes to a vertex or link from a dictionary.

        Parameters
        ----------
        obj: VertexSeq or EdgeSeq
            vertex or edge to add the attributes to
        d: dict
            attribute dictionary

        """

        for k, v in d.iteritems():
            obj[k] = v


    @staticmethod
    def _sort(t):
        """
        Returns ordered tuple.

        **Examples:**

        >>> print PowerNetwork._sort((2, 1))
        (1, 2)
        >>> print PowerNetwork._sort(("b", "a"))
        ('a', 'b')
        >>> print PowerNetwork._sort((2 * 1j, 1))
        Traceback (most recent call last):
        ...
        TypeError: no ordering relation is defined for complex numbers

        """

        if t[0] < t[1]:
            return t
        else:
            return t[1], t[0]


    ###############################################################################
    # ##                       FUNCTION ATTIC                                  ## #
    ###############################################################################

    @staticmethod
    def vertex_name_index(graph):
        raise NotImplementedError("Not supported yet.")
        # return dict((v, k) for k, v in enumerate(graph.vs["name"]))

    @staticmethod
    def graph_union_by_name(g1, g2):
        raise NotImplementedError("Not supported yet.")
        # result_names = sorted(set(g1.vs["name"]) + set(g2.vs["name"]))
        # result = Graph(len(result_names))
        # result.vs["name"] = "names"
        # result_name_index = "vertex_name_index(graph)""
        # edges = []
        # for edge in g1:
        #     src = result_name_index[g1.vs[edge.source]["name"]]
        #     dest = result_name_index[g1.vs[edge.target]["name"]]
        #     edges.append((src, dest))
        # for edge in g2:
        #     src = result_name_index[g2.vs[edge.source]["name"]]
        #     dest = result_name_index[g2.vs[edge.target]["name"]]
        #     edges.append((src, dest))
        # result.add_edges(edges)
        # return result


###############################################################################
###############################################################################


def main():
    # g = PowerNetwork("test", edgelist, v_dict, e_dict)
    # g = PowerNetwork("test", n=n, edgelist=edgelist, v_dict=v_dict, e_dict=e_dict)

    # import igraph as ig
    # f = ig.Graph.Erdos_Renyi(10, 0.5)
    # h = ig.Graph.Erdos_Renyi(10, 0.5)
    #
    # f.vs["name"] = range(10)
    # h.vs["name"] = range(5, 15)
    #
    # g = PowerNetwork.convert_igraph_object(f + h, "ER")
    # g = PowerNetwork("test", n=100)
    # g.vertex_filter("lat", 0.5, "_ge")
    #g.visualise("plot")

    g = PowerNetwork.SmallTestPowerGrid()
    print g

    # g.add_buses([10, 21], ["generator", "generator"])

    # g.add_branches([(0, 5), ])


    g.output_dir = "./"
    g.visualise("plot")
    g.save("data")


if __name__ == "__main__":
    main()
    #import doctest
    #doctest.testmod()
