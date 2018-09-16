import cantera as ct 
import numpy as np 
import pandas as pd
import time 
import matplotlib.pylab as plt 
import pdb 
from argparse import ArgumentParser
import os.path
import pyarrow.parquet as pq 
import pyarrow as pa
import multiprocessing
from pathos.multiprocessing import ProcessingPool
from CanteraTools import *

class Particle(object):
    """Class for particle in BatchPaSR.
    """
    gas_template = None
    
    @classmethod
    def fromGas(cls, gas, particle_mass = 1.0, chemistry = True):
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
        if Particle.gas_template == None:
            Particle.gas_template = ct.Solution(gas.name + ".xml")
        return cls(np.hstack([gas.T, gas.P, gas.X]), particle_mass, gas.name + ".xml", chemistry)
    
    def __init__(self, state, particle_mass = 1.0, mech='gri30.xml', chemistry = True):
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
        if Particle.gas_template == None:
            Particle.gas_template = ct.Solution(mech)        
        self.mech = mech
        self.P = state[1]
        Particle.gas_template.TPX = [state[0], state[1], state[2:]]
        self.column_names = ['age', 'T', 'MW', 'h'] + ["Y_" + sn for sn in Particle.gas_template.species_names] + ["X_" + sn for sn in Particle.gas_template.species_names]        
        self.mass = particle_mass
        self.age = 0
        self.state = np.hstack((Particle.gas_template.enthalpy_mass, Particle.gas_template.Y))
        self.timeHistory_list = [[self.age, Particle.gas_template.T, Particle.gas_template.mean_molecular_weight, Particle.gas_template.enthalpy_mass] + Particle.gas_template.Y.tolist() + Particle.gas_template.X.tolist()]
        self.timeHistory_array = None
        self.chemistry_enabled = chemistry
    
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
            self.state[0] = h
            self.state[1:] = Y
        else:
            return self.state
   
    def __add__(self, other):
        """Add values to state of particle without changing the state of either particle.

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
        """Add values to state of particle without changing the state of either particle.

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
            assert len(other) == len(self.state) , "Please ensure that input array has the same length as current particle state array ({0.d})".format(len(self.state))            
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
        """Subtract state of particle from input state without changing the state of either particle.

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
            h = other - self.state[0]
            Y = other - self.state[1:]
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
        """Add values to state of particle, modifying the state itself.

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
            h = (self.mass*self.state[0] + other.mass*other.state[0])/(self.mass + other.mass)
            Y = (self.mass*self.state[1:] + other.mass*other.state[1:])/(self.mass + other.mass)
            self.mass += other.mass
        elif isinstance(other, np.ndarray):
            assert len(other) == len(self.state) , "Please ensure that input array has the same length as current particle state array ({0.d})".format(len(self.state))
            h = self.state[0] + other[0]
            Y = self.state[1:] + other[1:]
        elif isinstance(other, (int, float)):
            h = self.state[0] + other
            Y = self.state[1:] + other
        else:
            return NotImplemented
        self.state[0] = h
        self.state[1:] = Y
        return self

    def __isub__(self, other):
        """Subtract values from state of particle, modifying the particle itself.

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
            h = self.state[0] - other.state[0]
            Y = self.state[1:] - other.state[1:]
        elif isinstance(other, np.ndarray):
            assert len(other) == len(self.state) , "Please ensure that input array has the same length as current particle state array ({0.d})".format(len(self.state))
            h = self.state[0] - other[0]
            Y = self.state[1:] - other[1:]
        elif isinstance(other, (int, float)):
            h = self.state[0] - other
            Y = self.state[1:] - other
        else:
            return NotImplemented
        self.state[0] = h
        self.state[1:] = Y
        return self

    def __imul__(self, other):
        """Multiply state of particle by value, modifying the particle itself.

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
            self.state *= other
            return self
        else:
            return NotImplemented
        
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
        Particle.gas_template.HPY = [self.state[0], self.P, self.state[1:]]
        self.timeHistory_list.append([self.age, Particle.gas_template.T, Particle.gas_template.mean_molecular_weight, Particle.gas_template.enthalpy_mass] + Particle.gas_template.Y.tolist() + Particle.gas_template.X.tolist())        
        reac = ct.ConstPressureReactor(Particle.gas_template,
            volume= self.mass/Particle.gas_template.density)
        reac.chemistry_enabled = self.chemistry_enabled
        netw = ct.ReactorNet([reac])
        netw.advance(dt)
        self.age += dt
        #         self.timeHistory_list = [[self.age, Particle.gas_template.T, Particle.gas_template.mean_molecular_weight, Particle.gas_template.enthalpy_mass] + Particle.gas_template.Y.tolist() + Particle.gas.X.tolist()]        
        self.timeHistory_list.append([self.age, Particle.gas_template.T, Particle.gas_template.mean_molecular_weight, Particle.gas_template.enthalpy_mass] + Particle.gas_template.Y.tolist() + Particle.gas_template.X.tolist())
        self.state = np.hstack((Particle.gas_template.enthalpy_mass, Particle.gas_template.Y))

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
            df.set_index(['age'])
            return df
        
        return self.timeHistory_array

    def __getstate__(self):
        # how to get the state data out of a Particle instance
        state = self.__dict__.copy()
        return state

    def __setstate__(self, state):
        # rebuild a Particle instance from state
        self.__dict__.update(state)


class PaSBR(object):
    def __init__(self, particle_list, N_MAX=10, dt = 0.01e-3, coalesce = True):
        """Initialize BatchPaSR from a list of Particles
        
        Parameters
        ----------
        particle_list : `List/numpy array of Particles`
            A python list/array of Particle objects
        
        Returns
        -------
        None 
        
        """ 
        
        self.column_names = ['age', 'mass', 'T', 'MW', 'h'] + ["Y_" + sn for sn in Particle.gas_template.species_names] + ["X_" + sn for sn in Particle.gas_template.species_names]        
        self.particle_list = particle_list
        self.N_MAX = N_MAX
        self.dt = dt # note: make sure dt is smaller than tau_mix!!! 
        self.ParticleFlowController = None
        self.time = 0
        self.timenext = 10*dt   # Next point in time to try and coalesce particles
        self.mass = 0
        self.state = 0
        self.N = len(self.particle_list)
        self.mean_gas = ct.Solution(Particle.gas_template.name + ".xml")
        self.mean_gas.name = "BatchPaSR Mean Gas"
        self.P = 0
        self.timeHistory_list = []
        if len(particle_list) > 0:
            self.P = particle_list[0].P # NOTE: Assume all particles have same temp
            self.mean_gas.HPY = particle_list[0].state[0], self.P, particle_list[0].state[1:]
            self.updateState()
            self.timeHistory_list = [[self.time, self.mass, self.mean_gas.T, self.mean_gas.mean_molecular_weight, self.mean_gas.enthalpy_mass] + self.mean_gas.Y.tolist() + self.mean_gas.X.tolist()]        
        # self.chemistry_enabled = chemistry

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
        if self.P == 0:
            self.P = self.particle_list[0].P
            self.mean_gas.TPX = self.mean_gas.T, self.P, self.mean_gas.X
        self.mass = sum([p.mass for p in self.particle_list]) # Note: can remove this if we're not updating particle_list if we don't usually update particle_list before we call this method
        self.N = len(self.particle_list)  
        self.state = sum([p.mass * p.state for p in self.particle_list])/self.mass # use p() or p.state?
        self.mean_gas.HPY = self.state[0], self.mean_gas.P, self.state[1:]
        assert all([round(particle.P) == round(self.P) for particle in self.particle_list]), "BatchPaSR does not support particles with different pressures (yet)"
        assert self.N <= self.N_MAX , f"N ({self.N}) > N_MAX ({self.N_MAX}); too many particles"
        self.timeHistory_list.append([self.time, self.mass, self.mean_gas.T, self.mean_gas.mean_molecular_weight, self.mean_gas.enthalpy_mass] + self.mean_gas.Y.tolist() + self.mean_gas.X.tolist())
        
        
    def insert(self, particle):
        self.particle_list.append(particle)
        self.updateState() # Is the cost of assigning N to a new value (in updateState()) higher than doing N += 1?
    
    @classmethod
    def reactHelper(cls, particle_dt):
        particle, dt = particle_dt
        particle.react(dt)
        return particle.__getstate__()
        
    def react(self, parallel=False, coalesce = True):
        if parallel:
            pool = multiprocessing.Pool(processes=4)
            jobList = [(p, self.dt) for p in self.particle_list]
#             states = ProcessingPool().map(BatchPaSR.reactHelper, jobList)
            states = pool.map(PaSBR.reactHelper, jobList)
            pool.close()
            [self.particle_list[i].__setstate__(states[i]) for i in range(0, len(self.particle_list))]
        else:
            [p.react(self.dt) for p in self.particle_list]

        # Optional coalescing of particles for performance
        if coalesce and self.time >= self.timenext:
            self.gobble()

        self.time += self.dt        
        self.updateState()
    
    def gobble(self):
        self.timenext += 50*self.dt
        tol = 1e-12
        ind = 1
        p1_ind = 0
        while p1_ind < len(self.particle_list):
            p1 = self.particle_list[p1_ind]
            # loop through rest of the particles in the list to gobble if possible
            particles_to_delete = []
            for i in range(1, len(self.particle_list)):
                p2 = self.particle_list[i]
                if self._canCombine(p1,p2):
                    p1 += p2
                    particles_to_delete.append(i)
            [del self.particle_list[ind] for ind in particles_to_delete] # delete particles 
            p1_ind += 1

        
        # while ind < len(self.particle_list):
        #     p0 = self.particle_list[0]
        #     p = self.particle_list[ind]
        #     diffH = p.state[0] - p0.state[0]
        #     diffY = p.state[1:] - p0.state[1:]
        #     if ((diffH/(p0.state[0] + np.finfo(np.float64).eps))**2 < tol and
        #         np.linalg.norm(np.divide(diffY, p0.state[1:] + np.finfo(np.float64).eps)) < tol):
        #         p0 += p # combine particles and update p0
        #         del self.particle_list[ind]
        #     else:
        #         ind += 1        
    
    def _canCombine(p1, p2, tol=1e-12):
        diffH = p.state[0] - p0.state[0]
        diffY = p.state[1:] - p0.state[1:]        
        return ((diffH/(p0.state[0] + np.finfo(np.float64).eps))**2 < tol and np.linalg.norm(np.divide(diffY, p0.state[1:] + np.finfo(np.float64).eps)) < tol)

#     def iem(cls, paticle_list)
    def _iem(self, particle_list, k):
        k_avg = k*self.state
        for p in particle_list:
            p.state = p * (k + 1) - k_avg
 
    def mix(self, tau_mix): # note: make sure dt < tau_mix!
        # Constant k:
        k = -self.dt/tau_mix # note: actually, k = -0.5*C_phi*omega*dt, but since C_phi is usually 2, i canceled it out.
        k_avg = k*self.state
        [p(p * (k + 1) - k_avg) for p in self.particle_list] # setting p.state_new = p.state_old + k*p.state_old - k*avg_state
        self.updateState()
    
    def prepEntrainment(self, added_gas, total_mass_added, tau_ent, numParticles=10, method='constant', time_interval = None):
        """Initialize particles to be entrained
        
        Parameters
        ----------
        added_gas : `cantera.Solution`
            Thermochemical state of gas to be entrained
        
        total_mass_added : `float`
            Total mass to be entrained 
            (Edwin note: could this be a list of masses which some up to some total mass instead?)
        
        tau_ent : `float`
            Total time for which entrainment occurs (in seconds)
        
        numParticles : `int`
            Number of particles to be entrained 

        time_interval : `numpy.array`
            Custom defined points in time to be used in entrainment as opposed to constant time steps
            Note: Entrainment mass is determined
        
        Returns
        -------
        None 
        
        """ 
                
        # Todo: add multiple entrainment methods here (defaulting to constant for now)
        if time_interval is None:
            self.entrain_timer = np.arange(0, tau_ent, tau_ent/numParticles)
        else:
            self.entrain_timer = time_interval
        entrainmentMass = []
        for i in range(len(self.entrain_timer)):
            prev = self.entrain_timer[i]
            next = tau_ent
            if i < len(self.entrain_timer)-1:
                next = self.entrain_timer[i+1]
            # Compare difference in time to last 
            mass = total_mass_added*(next-prev)/tau_ent
            entrainmentMass.append(mass)

        self.inactive_particles = []    # Reset for sanity

        # Initialize particles
        for i in range(0, numParticles):
            tempParticle = Particle.fromGas(added_gas, particle_mass = entrainmentMass[i], chemistry = self.chemistry_enabled)
            self.inactive_particles.append(tempParticle)

    def entrain(self, current_time):
        # Adds the next inactive particle if the time is correct
        if self.entrainInd < len(self.inactive_particles) and current_time >= self.entrain_timer[self.entrainInd]:
            # Note: The code right now will wait until the time is at or past the defined checkpoints
            # It will not entrain exactly at the defined time if the time steps do not synchronize
            # However, this should be less of an issue if more particles and smaller time steps are used
            current_particle = self.inactive_particles[self.entrainInd]
            current_particle.react(current_time)  # Continue reacting until entrainment. Note: assumes current_time starts from 0            
            self.particle_list.append(current_particle)
            self.entrainInd += 1
        

    def get_timeHistory(self, dataFrame=True):
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
            df.set_index(['age'])
            return df
        
        return self.timeHistory_array        

class ParticleFlowController(object):
    def __init__(self, reactor, gas, totalmass, timestep, method):
        """
        Parameters
        ----------
        reactor : `BatchPaSR`
            The reactor that this PFC acts as an inlet for
        
        gas : `cantera.Solution`
            Thermochemical state of gas to be entrained

        totalmass : `float`
            Total amount of mass that can be entrained in

        timestep : `float`
            Time between entrainment steps (continue to add particles if there's still mass left)

        method : `string` or `function`
            Function of time that describes entrainment rate (i.e. '0.5' or '0.12*t**2-0.6*t')
            The expression must be in terms of 't' in seconds and in Python syntax
            Alternatively, pass in a function to be called that only takes in a time and returns a mass flow rate
        
        Returns
        -------
        None
        """
        self.bp = reactor
        self.gas = gas
        self.rn = ct.ReactorNet([ct.ConstPressureReactor(self.gas)])
        self.mass = totalmass
        self.dt = timestep
        self.mdotexp = method
        # self.nextstep = timestep/2   # Which time to add the next particle (using trapezoidal sum)
        self.nextstep = 0

    def canEntrain(self, t):
        # Check if we can even entrain at this time
        return (self.mass > 0) and (self.nextstep <= t)
        
    def entrain(self, t):
        # Will add a particle to the reactor if entrainment is possible
        if self.canEntrain(t):
            # Calculate the mass of added particle
            if callable(self.mdotexp):
                mdot1 = self.mdotexp(self.nextstep - self.dt/2)
                mdot2 = self.mdotexp(self.nextstep + self.dt/2)
            else:
                mdot1 = eval(self.mdotexp.replace('t', str(self.nextstep - self.dt/2)))
                mdot2 = eval(self.mdotexp.replace('t', str(self.nextstep + self.dt/2)))
            self.nextstep += self.dt
            entrained_mass =  min(self.mass, (mdot1 + mdot2)/2 * self.dt) # trapezoidal rule; if remaining self.mass < mass to be entrained then entrain self.mass
            # Keep track of remaining mass
            
            self.mass -= entrained_mass
            self.rn.advance(t)
            particle = Particle.fromGas(self.gas, particle_mass=entrained_mass)
            # print(f"Entraining particle of mass {entrained_mass:.2f} kg into PaSBR");
            # pdb.set_trace()
            # Add it into the reactor
            self.bp.insert(particle)
            return self.mass*particle()
        else:
            return self.mass, np.hstack((gas.enthalpy_mass, gas.Y))