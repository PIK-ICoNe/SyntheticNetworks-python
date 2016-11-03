#!/usr/bin/python
# -*- coding: utf-8 -*-

from pgtem.components import *
from os.path import join, dirname


def main():

    #from scipy.spatial import kdtree
    # TODO: use kdtree to get (approximate) distances in a point cloud quickly

    ## init actors
    producer = Creator.source()
    consumer = Creator.sink()
    coord = Coordinator()
    goHV = Provider.HV()
    regulator = Regulator()

    ## init empty power grid with desired number of nodes
    net = PowerGrid(100)
    assert isinstance(net, PowerGrid)

    ####### initially, we might again have a MST ###########

    net.import_from_rpg(n=100, n0=10)

    print net.number_of_nodes, net.number_of_edges


    ####### growth mechanism #######

    for time in xrange(0, 20, 1):

        # propose new plant
        node = producer.new_node(net.number_of_nodes + 1, time)

        if net.number_of_nodes == 0:
            net.update(node, edgelist=[])
        else:
            # coordinator proposes possible extensions to connect new node
            possible_extensions = coord.proposition(node, net)

            # coordinator contacts providers and the regulator
            allowed, ranking = coord.negotiation(possible_extensions, [goHV,], regulator)


            if allowed:
                net.update(node, edgelist=ranking[0])

    print net.edges

    print net.nodes

    net.save(join(dirname(__file__), "networks", "test.network"))




if __name__ == '__main__':
    main()