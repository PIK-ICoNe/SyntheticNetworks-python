from actors import *
from evolution import *

def main():

    ## init acttors
    producer = Creator.source()
    consumer = Creator.sink()
    coord = Coordinator()
    goLV = Provider.LV()
    goMV = Provider.MV()
    goHV = Provider.HV()
    regulator = Regulator()

    ## iniit empty power grid
    net = PowerGrid()

    for time in xrange(0, 10, 1):

        # propose new plant
        loc = producer.new_node()

        # coordinate construction
        allowed, ranking = node_proposition(loc, coord, [goLV, goMV, goHV], regulator)


        if allowed:
            net.update(loc, ranking[0], time)


    net.save("test.network")




if __name__ == '__main__':
    main()