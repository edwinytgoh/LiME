"""
This example demonstrates how Cantera can be used with the 'multiprocessing'
module.

Because Cantera Python objects are built on top of C++ objects which cannot be
passed between Python processes, it is necessary to set up the computation so
that each process has its own copy of the relevant Cantera objects. One way to
do this is by storing the objects in (module) global variables, which are
initialized once per worker process.
"""

import multiprocessing
import numpy as np
import cantera as ct
import itertools
from time import time
import os

# Global storage for Cantera Solution objects
gases = {}
nets = {}
reactors = {}

def init_process(mech):
    """
    This function is called once for each process in the Pool. We use it to
    initialize any Cantera objects we need to use.
    """
    id = os.getpid()
    gases[id] = ct.Solution(mech)
    reactors[id] = ct.ConstPressureReactor(gases[id])
    nets[id] = ct.ReactorNet([reactors[id]])
    # gases[mech] = ct.Solution(mech)
    gases[id].transport_model = 'Multi'
    # gases[mech].transport_model = 'Multi'

def get_thermal_conductivity(args):
    # Pool.imap only permits a single argument, so we pack all of the needed
    # arguments into the tuple 'args'
    id = os.getpid()
    mech, T, P, X = args
    gas = gases[id]
    gas.TPX = T, P, X
    reactor = reactors[id]
    net = nets[id]
    reactor.syncState() # this will automatically call net.reinitialize()
    net.advance(10.000)  
    print(f"Hello from {id}. T_init = {T:.2f} K; T_out = {reactor.thermo.T:.2f} K\n")  
    return gas.thermal_conductivity

def get_viscosity(args):
    # Pool.imap only permits a single argument, so we pack all of the needed
    # arguments into the tuple 'args'
    id = os.getpid()
    mech, T, P, X = args
    gas = gases[id]
    gas.TPX = T, P, X
    reactor = reactors[id]
    net = nets[id]

    return gas.viscosity

def parallel(mech, predicate, nProcs, nTemps):
    """
    Call the function ``predicate`` on ``nProcs`` processors for ``nTemps``
    different temperatures.
    """
    P = 25*ct.one_atm
    X = 'CH4:1.0, O2:2.0'
    pool = multiprocessing.Pool(processes=nProcs,
                                initializer=init_process,
                                initargs=(mech,))

    y = pool.map(predicate,
                 zip(itertools.repeat(mech),
                     np.linspace(800, 1200, nTemps),
                     itertools.repeat(P),
                     itertools.repeat(X)))
    return y

def serial(mech, predicate, nTemps):
    P = ct.one_atm
    X = 'CH4:1.0, O2:2.0'
    init_process(mech)
    y = list(map(predicate,
                 zip(itertools.repeat(mech),
                     np.linspace(300, 900, nTemps),
                     itertools.repeat(P),
                     itertools.repeat(X))))
    return y

if __name__ == '__main__':
    multiprocessing.set_start_method("spawn")
    nPoints = 100
    nProcs = 8

    # For functions where the work done in each subprocess is substantial,
    # significant speedup can be obtained using the multiprocessing module.
    print('Thermal conductivity')
    t1 = time()
    parallel('gri30.xml', get_thermal_conductivity, nProcs, nPoints)
    t2 = time()
    print('Parallel: {0:.3f} seconds'.format(t2-t1))

    t1 = time()
    serial('gri30.xml', get_thermal_conductivity, nPoints)
    t2 = time()
    print('Serial: {0:.3f} seconds'.format(t2-t1))

    # On the other hand, if the work done per call to the predicate function is
    # small, there may be no advantage to using multiprocessing.
    print('\nViscosity')
    t1 = time()
    parallel('gri30.xml', get_viscosity, nProcs, nPoints)
    t2 = time()
    print('Parallel: {0:.3f} seconds'.format(t2-t1))

    t1 = time()
    serial('gri30.xml', get_viscosity, nPoints)
    t2 = time()
    print('Serial: {0:.3f} seconds'.format(t2-t1))