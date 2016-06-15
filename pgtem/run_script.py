from components import *

def main():

    ## init acttors
    producer = Creator.source()
    consumer = Creator.sink()
    coord = Coordinator()
    goHV = Provider.HV()
    regulator = Regulator()

    ## iniit empty power grid
    net = PowerGrid(10)

    for time in xrange(0, 10, 1):

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

    print net.nodes
    print net.edges


    net.save("test.network")




if __name__ == '__main__':
    main()