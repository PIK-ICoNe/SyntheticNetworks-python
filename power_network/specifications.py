#!/usr/bin/python
# -*- coding: utf-8 -*-

__version__ = "1.0"
__author__ = "Paul Schultz, pschultz@pik-potsdam.de"
__copyright__ = "Copyright (C) 2016 Paul Schultz, GNU GPL 3"


"""
Library of PowerNetwork.

Provides data structures for standard bus/branch types
and options to initialise default values.

"""

# use pint to handle physical units
# Conventions:
# - W for real power
# - var for reactive power
# - VA for apparent power

from pint import UnitRegistry
ureg = UnitRegistry()
ureg.define("var = W")

class Bus(object):
    """
    Generic bus class containing attributes common to all vertices.

    BusID <int>     : Numerical identifier
    Name <str>      : name
    Type <str>      : type of node
                        1) generator
                        2) motor
                        3) slack bus
    source <str>    : energy source
    operator <str>  : grid operator
    lon <float>     : longitude
    lat <float>     : latitude
    comments <str>  : open for comments
    """

    def __init__(self):
        self.Bus_ID = 0
        self.Name = "NaN"
        self.Type = "NaN"
        self.source = "NaN"
        self.comments = "NaN"
        self.operator = "NaN"
        self.lat = 0
        self.lon = 0
        self.volt = 0 * ureg("kV")

    def __setattr__(self, name, value):
        if hasattr(self, name):
            if isinstance(getattr(self, name), ureg.Quantity):
                assert value.dimensionality == getattr(self, name).dimensionality
            else:
                assert not isinstance(value, ureg.Quantity)
        super(Bus, self).__setattr__(name, value)

    def __str__(self):
        return "BUS: " + str(self.Bus_ID) + " " + str(self.Type)+" " + str(self.volt)

class Generator(Bus):
    """
    Specification for synchronous machine generators.

    For the definition of D and H, see the following publication:

    Auer, S., Kleis, K., Schultz, P., Kurths, J., & Hellmann, F. (2016).
    The impact of model detail on power grid resilience measures.
    The European Physical Journal Special Topics, 225(3), 609–625.
    http://doi.org/10.1140/epjst/e2015-50265-9


    PG <float>      : real power output
    QG <float>      : reactive power output
    QMax <float>    : max reactive power output
    QMin <float>    : min reactive power output
    MBase <float>   : total MVA base (apparent power)
    Status <int>    : machine status
                        > 0 in-service
                        < 0 off
    Pmax <float>    : max real power output [MW]
    Pmin <float>    : min real power output [MW]
    H <float>       : inertia time constant
    D <float>       : speed damping time constant
    """
    def __init__(self, **kwargs):
        super(self.__class__, self).__init__()

        self.PG = "NaN" * ureg("MW")
        self.QG = "NaN" * ureg("Mvar")
        self.QMax = "NaN" * ureg("Mvar")
        self.QMin = "NaN" * ureg("Mvar")
        self.MBase = "NaN" * ureg("MVA")
        self.Status = 0
        self.PMax = "NaN" * ureg("MW")
        self.PMin = "NaN" * ureg("MW")
        self.H = "NaN" * ureg("s")
        self.D = "NaN" * ureg("s")

        for key, val in kwargs.iteritems():
            if hasattr(self, key):
                setattr(self, key, val)
            else:
                raise AttributeError(key + " is not a proper generator attribute!")




    @classmethod
    def default(cls, **kwargs):
        new_cls = cls(**kwargs)
        new_cls.PG = 1. * ureg("MW")
        new_cls.QG = 1. * ureg("Mvar")
        new_cls.QMax = 2. * ureg("Mvar")
        new_cls.QMin = 0. * ureg("Mvar")
        new_cls.MBase = 1. * ureg("MVA")
        new_cls.Status = 1.
        new_cls.PMax = 2. * ureg("MW")
        new_cls.PMin = 0. * ureg("MW")
        new_cls.H = 1. * ureg("s")
        new_cls.D = 0.1 * ureg("s")

        return new_cls

class Load(Bus):
    """
    Specification for loads.

    PD <float>      : real power demand [MW]
    QD <float>      : reactive power demand
    GS <float>      : shunt conductance
    BS <float>      : shunt susceptance
    VM <float>      : voltage magnitude [kV]
    VA <float>      : voltage angle
    VBase <float>   : base voltage
    VMax <float>    : max voltage magnitude
    Vmin <float>    : min voltage magnitude
    """
    def __init__(self, **kwargs):
        super(self.__class__, self).__init__()
        self.PD = "NaN" * ureg("MW")
        self.QD = "NaN" * ureg("Mvar")
        self.GS = "NaN" * ureg("nS")
        self.BS = "NaN" * ureg("nS")
        self.VM = "NaN" * ureg("kV")
        self.VA = "NaN" * ureg("radian")
        self.VBase = "NaN" * ureg("kV")
        self.VMax = "NaN" * ureg("kV")
        self.VMin = "NaN" * ureg("kV")

        for key, val in kwargs.iteritems():
            if hasattr(self, key):
                setattr(self, key, val)
            else:
                raise AttributeError(key + " is not a proper load attribute!")

    @classmethod
    def default(cls, **kwargs):
        new_cls = cls(**kwargs)
        new_cls.PD = 1. * ureg("MW")
        new_cls.QD = 1. * ureg("Mvar")
        new_cls.GS = 1 * ureg("nS")
        new_cls.BS = 1 * ureg("nS")
        new_cls.VM = 220 * ureg("kV")
        new_cls.VA = 0.2 * ureg("radian")
        new_cls.VBase = 220 * ureg("kV")
        new_cls.VMax = 221 * ureg("kV")
        new_cls.VMin = 219 * ureg("kV")

        new_cls.volt = new_cls.VBase
        return new_cls

class Prosumer(Bus):
    """
    PG <float>      : real power output
    QG <float>      : reactive power output
    QMax <float>    : max reactive power output
    QMin <float>    : min reactive power output
    MBase <float>   : total MVA base (apparent power)
    Status <int>    : machine status
                        > 0 in-service
                        < 0 off
    Pmax <float>    : max real power output [MW]
    Pmin <float>    : min real power output [MW]
    H <float>       : inertia time constant
    D <float>       : speed damping time constant


    PD <float>      : real power demand [MW]
    QD <float>      : reactive power demand
    GS <float>      : shunt conductance
    BS <float>      : shunt susceptance
    VM <float>      : voltage magnitude [kV]
    VA <float>      : voltage angle
    VBase <float>   : base voltage
    VMax <float>    : max voltage magnitude
    Vmin <float>    : min voltage magnitude
    """
    def __init__(self, **kwargs):
        super(self.__class__, self).__init__()

        self.PG = "NaN" * ureg("MW")
        self.QG = "NaN" * ureg("Mvar")
        self.QMax = "NaN" * ureg("Mvar")
        self.QMin = "NaN" * ureg("Mvar")
        self.MBase = "NaN" * ureg("MVA")
        self.Status = 0
        self.PMax = "NaN" * ureg("MW")
        self.PMin = "NaN" * ureg("MW")
        self.H = "NaN" * ureg("s")
        self.D = "NaN" * ureg("s")

        self.kP = 1. * ureg("1/(W*s)")
        self.M = 1. * ureg("s**2*W")

        self.PD = "NaN" * ureg("MW")
        self.QD = "NaN" * ureg("Mvar")
        self.GS = "NaN" * ureg("nS")
        self.BS = "NaN" * ureg("nS")
        self.VM = "NaN" * ureg("kV")
        self.VA = "NaN" * ureg("radian")
        self.VBase = "NaN" * ureg("kV")
        self.VMax = "NaN" * ureg("kV")
        self.VMin = "NaN" * ureg("kV")

        for key, val in kwargs.iteritems():
            if hasattr(self, key):
                setattr(self, key, val)
            else:
                raise AttributeError(key + " is not a proper prosumer attribute!")




    @classmethod
    def default(cls, **kwargs):
        new_cls = cls(**kwargs)
        new_cls.PG = 1. * ureg("MW")
        new_cls.QG = 1. * ureg("Mvar")
        new_cls.QMax = 2. * ureg("Mvar")
        new_cls.QMin = 0. * ureg("Mvar")
        new_cls.MBase = 1. * ureg("MVA")
        new_cls.Status = 1.
        new_cls.PMax = 2. * ureg("MW")
        new_cls.PMin = 0. * ureg("MW")
        new_cls.H = 1. * ureg("s")
        new_cls.D = 1/30 * ureg("s")

        new_cls.kP = 1. * ureg("1/(W*s)")
        new_cls.M = 1. * ureg("s**2*W")


        new_cls.PD = 1. * ureg("MW")
        new_cls.QD = 1. * ureg("Mvar")
        new_cls.GS = 1 * ureg("nS")
        new_cls.BS = 1 * ureg("nS")
        new_cls.VM = 220 * ureg("kV")
        new_cls.VA = 0.2 * ureg("radian")
        new_cls.VBase = 20 * ureg("kV")
        new_cls.VMax = 221 * ureg("kV")
        new_cls.VMin = 219 * ureg("kV")

        new_cls.volt = new_cls.VBase

        return new_cls

class Passive(Bus):
    """
    Specification of a passive bus where nothing happens, just power flows through.
    """
    def __init__(self, **kwargs):
        #raise NotImplementedError("Not supported yet.")
        super(self.__class__, self).__init__()

        for key, val in kwargs.iteritems():
            if hasattr(self, key):
                setattr(self, key, val)
            else:
                raise AttributeError(key + " is not a proper Passive attribute!")


    @classmethod
    def default(cls, **kwargs):
        new_cls = cls(**kwargs)
        return new_cls

class Branch(object):
    """
    Generic branch class containing attributes common to all links.

    BranchID <int>     : Numerical identifier
    Name <str>      : name
    Type <str>      : type of branch
                        1) line
                        2) transformer
    operator <str>  : grid operator
    volt <float>    : operating voltage
    Bus1 <str>      : source bus name
    Bus2 <str>      : to bus  name
    source <int>    : source bus ID
    target <int>    : target bus  ID
    comments <str>  : open for comments
    """
    def __init__(self):
        self.Branch_ID = 0
        self.Name = "NaN"
        self.Type = "NaN"
        self.operator = "NaN"
        self.volt = 0 * ureg("kV")
        self.source = 0
        self.target = 0
        self.Bus1 = "NaN"
        self.Bus2 = "NaN"
        self.comments = "NaN"

    def __setattr__(self, name, value):
        if hasattr(self, name):
            if isinstance(getattr(self, name), ureg.Quantity):
                assert value.dimensionality == getattr(self, name).dimensionality
            else:
                assert not isinstance(value, ureg.Quantity)
        super(Branch, self).__setattr__(name, value)



class Line(Branch):
    """
    Transmission line.

    cables <int>    : number of electrically separated power-carrying conductors in a power line
    wires <int>     : number of wires per power cable
    R <float>       : resistance [ohm/km]
    X <float>       : reactance [ohm/km]
    C <float>       : capacity [nF/km]
    G <float>       : conductance [µS/km]
    B <float>       : susceptance [µS/km]
    I <float>       : max. threshold current [A]
    L <float>       : line length [km]

    """
    def __init__(self, **kwargs):
        super(self.__class__, self).__init__()
        self.cables = 0
        self.wires = 0
        self.R = 0 * ureg("ohm") / ureg("km")
        self.X = 0 * ureg("ohm") / ureg("km")
        self.C = 0 * ureg("nF") / ureg("km")
        self.G = 0 * ureg("microS") / ureg("km")
        self.B = 0 * ureg("microS") / ureg("km")
        self.I = 0 * ureg("A")
        self.L = 0 * ureg("km")

        for key, val in kwargs.iteritems():
            if hasattr(self, key):
                setattr(self, key, val)
            else:
                raise AttributeError(key + " is not a proper line attribute!")

    @classmethod
    def default(cls, **kwargs):
        new_cls = cls(**kwargs)
        new_cls.cables = 3
        new_cls.wires = 1
        new_cls.R = 0.4 * ureg("ohm") / ureg("km")
        new_cls.X = 0.3 * ureg("ohm") / ureg("km")
        new_cls.C = 0 * ureg("nF") / ureg("km")
        new_cls.G = 1 * ureg("microS") / ureg("km")
        new_cls.B = 1 * ureg("microS") / ureg("km")
        new_cls.I = 10 * ureg("A")
        new_cls.L = 27 * ureg("km")
        return new_cls

class Transformer(Branch):
    """
    Transformer substation.

    Tap <float>     : transformer off nominal turns ratio
    Shift <float>   : transformer phase shift angle
    Amin <float>    : minimum angle difference
    Amax <float>    : maximum angle difference

    """

    def __init__(self, **kwargs):
        super(self.__class__, self).__init__()
        self.Tap = 0
        self.Shift = 0 * ureg("radian")
        self.AMin = 0 * ureg("radian")
        self.AMax = 0 * ureg("radian")

        for key, val in kwargs.iteritems():
            if hasattr(self, key):
                setattr(self, key, val)
            else:
                raise AttributeError(key + " is not a proper transformer attribute!")

    @classmethod
    def default(cls, **kwargs):
        new_cls = cls(**kwargs)
        new_cls.Tap = 2
        new_cls.Shift = 0 * ureg("radian")
        new_cls.AMin = 0 * ureg("radian")
        new_cls.AMax = 0 * ureg("radian")
        return new_cls

