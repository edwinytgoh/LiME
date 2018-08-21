# For the Glory of God
import cantera as ct 
import numpy as np 
import pandas as pd
import time 
import matplotlib.pylab as plt 
import pdb 
# from argparse import ArgumentParser
# import os.path
# import pyarrow.parquet as pq
# import pyarrow as pa
import multiprocessing
# from pathos.multiprocessing import ProcessingPool

class Particle(object):
    """Class for particle in BatchPaSR.
    """
    
    @classmethod
    def fromGas(cls, gas, particle_mass = 1.0):
        """Initialize particle object with thermochemical state.
        
        Parameters
        ----------
        gas : `cantera.Solution`
            Initial thermochemical state of particle
        
        particle_mass : `float`
            Particle mass
        
        Returns
        -------
        cls : `Particle`
            Instance of the Particle class
        """    
#         pdb.set_trace()
        return cls(np.hstack([gas.T, gas.P, gas.X]), particle_mass, gas.name + ".xml")
    
    def __init__(self, state, particle_mass = 1.0, mech='gri30.xml'):
        """Initialize particle object with thermochemical state.

        Parameters
        ----------
        
        state : `numpy array or list`
            Initial thermochemical state of particle. Needs T, P, X. 
        
        gas : `cantera.Solution`
            Initial thermochemical state of particle
        
        particle_mass : `float`
            Particle mass
        
        Returns
        -------
        None

        """
        self.mech = mech;
        self.P = state[1]
        self.gas = ct.Solution(mech); 
        self.gas.TPX = [state[0], state[1], state[2:]]
        self.column_names = ['age', 'T', 'MW', 'h'] + ["Y_" + sn for sn in self.gas.species_names] + ["X_" + sn for sn in self.gas.species_names]        
        self.mass = particle_mass
        self.age = 0
        self.state = np.hstack((self.gas.enthalpy_mass, self.gas.Y))
#         self.time_history = pd.DataFrame(columns = self.column_names, data = np.hstack((0, gas.T, )))
        self.timeHistory_list = [[self.age, self.gas.T, self.gas.mean_molecular_weight, self.gas.enthalpy_mass] + self.gas.Y.tolist() + self.gas.X.tolist()]
        self.timeHistory_array = None; 
    
    def __call__(self, comp=None):
        """Return or set composition.
        Parameters
        ----------
        comp : Optional[cantera.Solution]

        Returns
        -------
        comp : numpy.array
            Thermochemical composition of particle (enthalpy + mass fractions).

        """
        if comp is not None:
            if isinstance(comp, Particle):
                h = comp.gas.enthalpy_mass
                Y = comp.gas.Y
            elif isinstance(comp, np.ndarray):
                h = comp[0]
                Y = comp[1:]
            elif isinstance(comp, ct.Solution):
                h = comp.enthalpy_mass
                Y = comp.Y
            else:
                return NotImplemented
            self.gas.HPY = h, self.gas.P, Y
        self.state = np.hstack((self.gas.enthalpy_mass, self.gas.Y))
        return self.state

    def __add__(self, other):
        """Add values to state of particle.

        Parameters
        ----------
        other : `Particle`, `numpy.array`, `int`, `float`
            Thermochemical state (enthalpy + mass fractions) to add to current state.

        Returns
        -------
        comp : numpy.array
            Thermochemical composition of particle (enthalpy + mass fractions).

        """
        if isinstance(other, Particle):
            return self.state + other.state
        elif isinstance(other, np.ndarray):
            return self.state + other
        elif isinstance(other, (int, float)):
            return self.state + other
        else:
            return NotImplemented

    def __radd__(self, other):
        """Add values to state of particle.

        Parameters
        ----------
        other : `Particle`, `numpy.array`, `int`, `float`
            Thermochemical state (enthalpy + mass fractions) to add to current state.

        Returns
        -------
        comp : numpy.array
            Thermochemical composition of particle (enthalpy + mass fractions).

        """
        if isinstance(other, Particle):
            return self.state + other.state
        elif isinstance(other, np.ndarray):
            return self.state + other
        elif isinstance(other, (int, float)):
            return self.state + other
        else:
            return NotImplemented

    def __sub__(self, other):
        """Subtract values from state of particle.

        Parameters
        ----------
        other : `Particle`, `numpy.array`, `int`, `float`
            Thermochemical state (enthalpy + mass fractions) to subtract from current state.

        Returns
        -------
        comp : numpy.array
            Thermochemical composition of particle (enthalpy + mass fractions).

        """
        if isinstance(other, Particle):
            return self.state - other.state
        elif isinstance(other, np.ndarray):
            return self.state - other
        elif isinstance(other, (int, float)):
            return self.state - other
        else:
            return NotImplemented

    def __rsub__(self, other):
        """Subtract state of particle from input state.

        Parameters
        ----------
        other : `Particle`, `numpy.array`, `int`, `float`
            Thermochemical state from which to subract Particle state.

        Returns
        -------
        comp : numpy.array
            Thermochemical composition of particle (enthalpy + mass fractions).

        """
        if isinstance(other, Particle):
            return other.state - self.state
        elif isinstance(other, np.ndarray):
            return other - self.state
        elif isinstance(other, (int, float)):
            h = other - self.gas.enthalpy_mass
            Y = other - self.gas.Y
            return np.hstack((h, Y))
        else:
            return NotImplemented

    def __mul__(self, other):
        """Multiply state of particle by value.

        Parameters
        ----------
        other : `int`, `float`
            Value to multiply `Particle` state by.

        Returns
        -------
        comp : numpy.array
            Thermochemical composition of particle (enthalpy + mass fractions).

        """
        if isinstance(other, (int, float)):
            return self.state * other
        else:
            return NotImplemented

    def __rmul__(self, other):
        """Multiply state of particle by value.

        Parameters
        ----------
        other : `int`, `float`
            Value to multiply `Particle` state by.

        Returns
        -------
        comp : numpy.array
            Thermochemical composition of particle (enthalpy + mass fractions).

        """
        if isinstance(other, (int, float)):
            return self.state * other
        else:
            return NotImplemented

    def __iadd__(self, other):
        """Add values to state of particle.

        Parameters
        ----------
        other : `Particle`, `numpy.array`, `int`, `float`
            Thermochemical state (enthalpy + mass fractions) to add to current state.

        Returns
        -------
        comp : numpy.array
            Thermochemical composition of particle (enthalpy + mass fractions).

        """
        if isinstance(other, Particle):
            h = self.gas.enthalpy_mass + other.gas.enthalpy_mass
            Y = self.gas.Y + other.gas.Y
        elif isinstance(other, np.ndarray):
            h = self.gas.enthalpy_mass + other[0]
            Y = self.gas.Y + other[1:]
        elif isinstance(other, (int, float)):
            h = self.gas.enthalpy_mass + other
            Y = self.gas.Y + other
        else:
            return NotImplemented
        self.gas.HPY = h, self.gas.P, Y
        return self

    def __isub__(self, other):
        """Subtract values from state of particle.

        Parameters
        ----------
        other : `Particle`, `numpy.array`, `int`, `float`
            Thermochemical state (enthalpy + mass fractions) to subtract from current state.

        Returns
        -------
        comp : numpy.array
            Thermochemical composition of particle (enthalpy + mass fractions).

        """
        if isinstance(other, Particle):
            h = self.gas.enthalpy_mass - other.gas.enthalpy_mass
            Y = self.gas.Y - other.gas.Y
        elif isinstance(other, np.ndarray):
            h = self.gas.enthalpy_mass - other[0]
            Y = self.gas.Y - other[1:]
        elif isinstance(other, (int, float)):
            h = self.gas.enthalpy_mass - other
            Y = self.gas.Y - other
        else:
            return NotImplemented
        self.gas.HPY = h, self.gas.P, Y
        return self

    def __imul__(self, other):
        """Multiply state of particle by value.

        Parameters
        ----------
        other : `int`, `float`
            Value to multiply `Particle` state by.

        Returns
        -------
        comp : numpy.array
            Thermochemical composition of particle (enthalpy + mass fractions).

        """
        if isinstance(other, (int, float)):
            h = self.gas.enthalpy_mass * other
            Y = self.gas.Y * other
        else:
            return NotImplemented
        self.gas.HPY = h, self.gas.P, Y
        return self

    def react(self, dt):
        """Perform reaction timestep by advancing network.

        Parameters
        ----------
        dt : float
            Reaction timestep [seconds]

        Returns
        -------
        None

        """
        reac = ct.ConstPressureReactor(self.gas,
            volume= self.mass/self.gas.density)
        netw = ct.ReactorNet([reac])
        netw.advance(dt)
        self.age += dt
        self.timeHistory_list.append([[self.age, self.gas.T, self.gas.mean_molecular_weight, self.gas.enthalpy_mass] + self.gas.Y.tolist() + self.gas.X.tolist()])
        self.state = np.hstack((self.gas.enthalpy_mass, self.gas.Y))

    def get_timeHistory(self, dataFrame=False):
        """Obtain particle's history. 
        
        Parameters
        ----------
        None
        
        Returns
        -------
        timeHistory_array : `numpy.array` 
            The array containing particle property time traces 
        """
        self.timeHistory_array = np.vstack(self.timeHistory_list)
        if dataFrame == True:
            df = pd.DataFrame(columns = self.column_names, data = self.timeHistory_array)
            df.set_index(['age']); 
            return df
        
        return self.timeHistory_array

    def __getstate__(self):
        # how to get the state data out of a Particle instance
        state = self.__dict__.copy()
        del state['gas']
        return state

    def __setstate__(self, state):
        # rebuild a Particle instance from state
        self.__dict__.update(state)
        self.gas = ct.Solution(self.mech)
        self.gas.HPY = self.state[0], self.P, self.state[1:]

class BatchPaSR(object):
    def __init__(self, particle_list, N_MAX=10, dt = 0.01e-3):
        """Initialize BatchPaSR from a list of Particles
        
        Parameters
        ----------
        particle_list : `List/numpy array of Particles`
            A python list/array of Particle objects
        
        Returns
        -------
        None 
        
        """ 
        
        self.column_names = ['age', 'T', 'MW', 'h'] + ["Y_" + sn for sn in particle_list[0].gas.species_names] + ["X_" + sn for sn in particle_list[0].gas.species_names]        
        self.P = particle_list[0].gas.P # NOTE: Assume all particles have same temp
        self.particle_list = particle_list
        self.N_MAX = N_MAX
        self.dt = 0.01e-3
        self.ParticleFlowController = None
        self.time = 0
        self.mass = 0
        self.state = None
        self.N = len(self.particle_list)
        self.mean_gas = ct.Solution(particle_list[0].gas.name + ".xml")
        self.mean_gas.name = "BatchPaSR Mean Gas"
        self.mean_gas.TPX = particle_list[0].gas.T, particle_list[0].gas.P, particle_list[0].gas.X
        self.updateState()
        self.timeHistory_list = [[self.time, self.mass, self.mean_gas.T, self.mean_gas.mean_molecular_weight, self.mean_gas.enthalpy_mass] + self.mean_gas.Y.tolist() + self.mean_gas.X.tolist()]        
        self.massList = [self.mass]

    def __call__(self):
        self.mean_gas()
        return self.state
    
    def updateState(self):
        """Update mass, number of particles, and BatchPaSR state vector
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
        self.mass = sum([p.mass for p in self.particle_list]) # Note: can remove this if we're not updating particle_list if we don't usually update particle_list before we call this method
        self.N = len(self.particle_list)  
        self.time += self.dt
        self.state = sum([p.mass * p() for p in self.particle_list])/self.mass
        self.mean_gas.HPY = self.state[0], self.mean_gas.P, self.state[1:]
        assert all([round(particle.gas.P) == round(self.P) for particle in self.particle_list]), "BatchPaSR does not support particles with different pressures (yet)"
        assert self.N <= self.N_MAX , "N > N_MAX; too many particles"
        
        
    def insert(self, particle):
        self.particle_list.append(particle); 
        self.updateState() # Is the cost of assigning N to a new value (in updateState()) higher than doing N += 1?
    
    @classmethod
    def reactHelper(cls, particle_dt):
        particle, dt = particle_dt
        particle.react(dt)
        return particle.__getstate__()
        
    
    def react(self, parallel=True):
        if parallel:
            pool = multiprocessing.Pool(processes=4)
            jobList = [(p, self.dt) for p in self.particle_list]
            # states = ProcessingPool().map(BatchPaSR.reactHelper, jobList)
            states = pool.map(BatchPaSR.reactHelper, jobList)
            pool.close()
            # pool.join()
            [self.particle_list[i].__setstate__(states[i]) for i in range(0, len(self.particle_list))]
        else:
            # [p.react(self.dt) for p in self.particle_list]
            for p in self.particle_list:
                p.react(self.dt)
        self.updateState()
        self.timeHistory_list.append([self.time, self.mass, self.mean_gas.T, self.mean_gas.mean_molecular_weight, self.mean_gas.enthalpy_mass] + self.mean_gas.Y.tolist() + self.mean_gas.X.tolist())

# if __name__ == '__main__':
#     # for j in range(1,16):
#     timeList = []
#     for i in range(0,2):
#         g1 = ct.Solution('gri30.xml');
#         g2 = ct.Solution('gri30.xml');
#         g1.TPX = 1500, 25*ct.one_atm, g1.X;
#         g2.TPX = 1200, 25*ct.one_atm, g2.X;
#         g1.set_equivalence_ratio(0.8, {'CH4':1.0}, {'N2':0.79, 'O2':0.21})
#         g2.set_equivalence_ratio(0.8, {'CH4':1.0}, {'N2':0.79, 'O2':0.21})
#         p1 = Particle.fromGas(g1);
#         p2 = Particle.fromGas(g2);
#         particle_list = [Particle.fromGas(g1) for _ in range(0,4)]
#         t1 = time.time();
#         t = np.arange(0, 10*1e-3, 0.01*1e-3)
#         bp = BatchPaSR(particle_list, dt=t[1]-t[0], N_MAX = 20);
#         [bp.react(parallel=True) for _ in t]
#         t2 = time.time();
#         timeList.append(t2-t1)
#         print("Run {0:d}: {1:.2f} seconds".format(i, t2-t1))
#         # bp.particle_list[0].get_timeHistory(dataFrame=True).to_clipboard();
#         print("Average time taken to run BatchPaSR of size {0:02d} for {1:04d} timesteps in series: {2:06.3f} seconds".format(bp.N, len(t), np.mean(np.array(timeList))))
    # pdb.set_trace()
    # bp.insert(p2)
    # [bp.react() for _ in t]
    # df = bp.particle_list[0].get_timeHistory(dataFrame=True)
    # df2 = bp.particle_list[1].get_timeHistory(dataFrame=True)
    # df2.to_clipboard()
    # df.append(bp.particle_list[0].get_timeHistory(dataFrame=True))


if __name__ == '__main__':
    timeList = []
    for i in range(0,10):
        g1 = ct.Solution('gri30.xml');
        g2 = ct.Solution('gri30.xml');
        g1.TPX = 1500, 25*ct.one_atm, g1.X;
        g2.TPX = 1200, 25*ct.one_atm, g2.X;
        g1.set_equivalence_ratio(0.8, {'CH4':1.0}, {'N2':0.79, 'O2':0.21})
        g2.set_equivalence_ratio(0.8, {'CH4':1.0}, {'N2':0.79, 'O2':0.21})
        p1 = Particle.fromGas(g1);
        p2 = Particle.fromGas(g2);

        t1 = time.time();
        t = np.arange(0, 10*1e-3, 0.01*1e-3)
        bp = BatchPaSR([p1, p1, p1, p1, p1, p1, p1, p1], dt=t[1]-t[0]);
        [bp.react(parallel=True) for _ in t]
        t2 = time.time();
        timeList.append(t2-t1)
        print("Run {0:d}: {1:.2f} seconds".format(i, t2-t1))
        # bp.particle_list[0].get_timeHistory(dataFrame=True).to_clipboard();
    print("Average time taken to run BatchPaSR of size {0:02d} for {1:04d} timesteps in series: {2:06.3f} seconds".format(bp.N, len(t), np.mean(np.array(timeList))))
