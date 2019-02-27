import itertools
import multiprocessing
import os
import os.path
import pdb
import time
from argparse import ArgumentParser

import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import pyarrow as pa
import pyarrow.parquet as pq
from CanteraTools import *
from Particle import Particle
from pathos.multiprocessing import ProcessingPool

multiprocessing.set_start_method("spawn")

NUM_PROCESS = 4
particles = {}
gases = {}
def init_process(mech, P):
    id = os.getpid()
    gases[id] = ct.Solution(mech, name=mech[0:-4])
    gases[id].TPX = 300, P, {"N2":1.0}
    particles[id] = Particle.fromGas(gases{id})
    pass

class PaSBR(object):
    def __init__(self, particle_list, N_MAX=10, dt = 0.01e-3, coalesce = True, mech = "gri30.xml", parallel = False):
        """Initialize BatchPaSR from a list of Particles
        
        Parameters
        ----------
        particle_list : `List/numpy array of Particles`
            A python list/array of Particle objects
        
        Returns
        -------
        None 
        
        """ 
        if Particle.gas_template == None:
            Particle.gas_template = ct.Solution(mech)
        self.column_names = ['age', 'mass', 'T', 'MW', 'h', 'phi'] + ["Y_" + sn for sn in Particle.gas_template.species_names] + ["X_" + sn for sn in Particle.gas_template.species_names]        
        self.particle_list = particle_list
        self.N_MAX = N_MAX
        self.dt = dt # note: make sure dt is smaller than tau_mix!!! 
        self.ParticleFlowController = None
        self.time = 0.0
        self.timenext = 10*dt   # Next point in time to try and coalesce particles
        self.mass = 0.0
        self.state = 0.0
        self.N = len(self.particle_list)
        self.mech = mech
        # assert Particle.gas_template.name == mech[0:-4] # TODO: makesure mechs match in all objects
        self.mean_gas = ct.Solution(self.mech)
        self.mean_gas.name = "BatchPaSR Mean Gas"
        self.P = 0.0
        self.timeHistory_list = []
        self.particle_timeHistory_list = []
        self.particle_timeHistory_info = []
        self.pool = multiprocessing.Pool(processes=NUM_PROCESS, 
                                        initializer=init_process,
                                        initargs=(self.mech, self.P))
        if len(particle_list) > 0:
            self.P = particle_list[0].P # NOTE: Assume all particles have same temp
            self.mean_gas.HPY = particle_list[0].state[0], self.P, particle_list[0].state[1:]
            self.updateState()
            self.timeHistory_list = [[self.time, self.mass, self.mean_gas.T, self.mean_gas.mean_molecular_weight, self.mean_gas.enthalpy_mass, self.mean_gas.get_equivalence_ratio()] + self.mean_gas.Y.tolist() + self.mean_gas.X.tolist()]        
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
        self.timeHistory_list.append([self.time, self.mass, self.mean_gas.T, self.mean_gas.mean_molecular_weight, self.mean_gas.enthalpy_mass, self.mean_gas.get_equivalence_ratio()] + self.mean_gas.Y.tolist() + self.mean_gas.X.tolist())
        
    def insert(self, particle):
        self.particle_list.append(particle)
        self.updateState() # Is the cost of assigning N to a new value (in updateState()) higher than doing N += 1?
    
    @classmethod
    def reactHelper(cls, inputs):
        state, dt = inputs
        particle = particles[os.getpid()]
        particle.__setstate__(state)
        particle.react(dt)
        return particle.__getstate__()


    def react(self, parallel=False, coalesce = True):
        if parallel:
            if self.pool == None:
                self.pool = multiprocessing.Pool(processes=NUM_PROCESS, 
                                                initializer=init_process,
                                                initargs=(self.mech, self.P))            
            current_states = [p.state for p in self.particle_list]
            helper_input = zip(current_states, itertools.repeat(self.dt))
            new_states = pool.map(PaSBR.reactHelper, helper_input)
            [p.__setstate__(new_states[i]) for i,p in enumerate(self.particle_list)]
            [self.particle_list[i].__setstate__(states[i]) for i in range(0, len(self.particle_list))]
        else:
            [p.react(self.dt) for p in self.particle_list]

        # Optional coalescing of particles for performance
        if coalesce:
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
            for i in range(p1_ind + 1, len(self.particle_list)):
                p2 = self.particle_list[i]
                if self._canCombine(p1,p2):
                    # p1 += p2 # this does not work; p1 is an entirely new object
                    self.particle_list[p1_ind] += p2
                    particles_to_delete.append(i)
            particles_to_delete = np.array(particles_to_delete)
            # do two loops because we don't want to affect the first loop
            # edwin note: can do one loop if the first loop is a while loop, i.e. while p2_ind < len(self.particle_list), and p2_ind doesn't change if particles_to_delete contains p2_ind, because the next p2 will be at the same position in the new particle_list
            # there's probably an algorithm for this in an algorithms textbook...
            for i in range(0,len(particles_to_delete)):
                ind = particles_to_delete[i]
                if ind < len(self.particle_list):
                    self.particle_timeHistory_list.append(self.particle_list[ind].get_timeHistory())
                    self.particle_timeHistory_info.append([len(self.particle_timeHistory_list[-1]), self.time, self.particle_list[ind].age])
                    del self.particle_list[ind] # delete particles 
                    particles_to_delete -= 1 # the list is getting shorter, so adjust the index of particles to be deleted accordingly
            self.updateState() # WE DIDN'T HAVE THIS BEFORE, SO THE BATCHPASR DIDN'T KNOW TO UPDATE IT'S MASS/STATE
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
    
    def _canCombine(self, p1, p2, tol=1e-6):
        H_1 = p1.state[0]
        Y_1 = p1.state[1:]
        diffH = p2.state[0] - p1.state[0]
        diffY = p2.state[1:] - p1.state[1:]
        machine_epsilon = np.finfo(np.float64).eps
        diffH_percent = (diffH/(H_1 + machine_epsilon))**2
        diffY_percent = np.linalg.norm(np.divide(diffY, Y_1 + machine_epsilon))
        diffH_is_small = diffH_percent < tol
        diffY_is_small = diffY_percent < tol        
        # diffY_is_small = np.linalg.norm(np.divide(diffY, p0.state[1:] + np.finfo(np.float64).eps)) < tol
        if (diffH_is_small and diffY_is_small):
            print(f"We are combining particles at t = {self.time/1e-3:.3f} ms! diffH_norm = {diffH_percent:.5E}, diffY_percent = {diffY_percent:.5E}, system mass = {self.mass:.3f} kg")		
        
        return ( diffH_is_small and diffY_is_small )

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

    def get_particleTimeHistory(self, dataFrame=True):
        self.updateState()
        for i in range(0, len(self.particle_list)):
            self.particle_timeHistory_list.append(self.particle_list[i].get_timeHistory())
            self.particle_timeHistory_info.append([len(self.particle_timeHistory_list[-1]), self.time, self.particle_list[i].age])
        particle_timeHistory_array = np.vstack(self.particle_timeHistory_list)
        if dataFrame == True:
            df = pd.DataFrame(columns = self.particle_list[-1].column_names, data=particle_timeHistory_array)
            # {key:value for key, value in zip({'particle_length', 'pasr_time', 'particle_age'}, [particle_timeHistory_array[:,i] for i in range(0,particle_timeHistory_array.shape[1])])}
            particle_timeHistory_info_df = pd.DataFrame(columns=['particle_length', 'pasr_time', 'particle_age'], data=np.vstack(self.particle_timeHistory_info))
            # particle_timeHistory_info_df.set_index('particle_age')
            # pdb.set_trace()
            particle_end_ind = particle_timeHistory_info_df['particle_length'].cumsum().values - 1
            particle_start_ind = np.hstack([0, particle_timeHistory_info_df['particle_length'].iloc[:-1].cumsum().values])
            particle_timeHistory_info_df['particle_end_ind'] = pd.Series(particle_end_ind, dtype='int64', index=particle_timeHistory_info_df.index)
            particle_timeHistory_info_df['particle_start_ind'] = pd.Series(particle_start_ind, dtype='int64', index=particle_timeHistory_info_df.index)

            # df.set_index(['age'])
            return df, particle_timeHistory_info_df
        
        return particle_timeHistory_array, particle_timeHistory_info

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
        return (self.mass > 0) and (self.nextstep <= t) # EDWIN COMMENT: Is nextstep <= t necessary?
        
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
            return self.mass, particle()
        else:
            return self.mass, np.hstack((self.gas.enthalpy_mass, self.gas.Y))
    
if __name__ == "__main__":
    multiprocessing.set_start_method("spawn")
    gas = ct.Solution("gri30.xml");
    p1 = Particle.fromGas(gas);
    gas()
