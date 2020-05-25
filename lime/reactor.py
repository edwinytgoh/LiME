import itertools
import multiprocessing
import os.path

from .particle import Particle
from .utils.cantera_tools import *

NUM_PROCESS = 8
particles = {}


def init_process(mech, P):
    id = os.getpid()
    print(f"Running init_process on PID = {id}")
    particles[mech] = Particle(mech, name=mech[0:-4])
    particles[mech].TPX = 300, P, {"N2":1.0}


class LiME(Particle):

    @classmethod
    def from_gas(cls, gas, mech="", particle_mass = 0.0001, N_MAX=10):
        """Initialize LiME object with thermochemical state.
        
        Parameters
        ----------
        gas : `cantera.Solution`
            Initial thermochemical state of particle. 
            Must have NAME property if using custom mechanism. Use ct.Solution("mech.xml", name="mech")
        
        particle_mass : `float`
            Particle mass
        
        Returns
        -------
        cls : `Particle`
            Instance of the Particle class
        """
        state_vec = gas.state
        if mech == "":
            if isinstance(gas, cls):
                mech = gas.mech
            else:
                mech = f"{gas.name}.xml"
        return cls(infile=mech, state_vec=state_vec, particle_mass=particle_mass, P = gas.P, N_MAX=N_MAX)

    def __init__(self, infile, particle_list = [], particle_mass=1.0, state_vec=[], P=101325, N_MAX=10, dt = 0.01e-3, coalesce = True, parallel = False, **kwargs):
        """Initialize LiME reactor from a list of Particles
        
        Parameters
        ----------
        infile
        particle_mass
        state_vec
        P
        N_MAX
        dt
        coalesce
        parallel
        kwargs
        particle_list : `List/numpy array of Particles`
            A python list/array of Particle objects
        
        Returns
        -------
        None 
        
        """
        super().__init__(infile=infile, particle_mass=particle_mass, state_vec=state_vec, P=P)
        self.particle_list = particle_list
        self.N_MAX = N_MAX
        self.dt = dt # ! note: make sure dt is smaller than tau_mix!!!
        self.ParticleFlowController = None
        self.time = 0.0
        self.timenext = 10*dt   # Next point in time to try and coalesce particles
        # assert Particle.gas_template.name == mech[0:-4] # TODO: makesure mechs match in all objects
        self.particle_timeHistory_list = []
        self.particle_timeHistory_info = []
        try:
            multiprocessing.set_start_method("spawn")
        except Exception:
            pass
        self.pool = multiprocessing.Pool(processes=NUM_PROCESS,
                                        initializer=init_process,
                                        initargs=(infile, P))
        if len(particle_list) > 0:
            self.HPY = particle_list[0].HPY
            self.update_state()
        # self.chemistry_enabled = chemistry
        self.entrain_ind = 0

    @property
    def N(self):
        return len(self.particle_list)

    # def __call__(self):
    #     self.mean_gas()
    #     return self.state

    def update_state(self):
        """Update mass, number of particles, and BatchPaSR state vector
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
        if self.P != self.particle_list[0].P:
            self.TPX = self.T, self.particle_list[0].P, self.X

        hy = np.zeros(self.HY.shape)
        m = 0
        for i,p in enumerate(self.particle_list):
            hy += p.mass * p.HY
            m += p.mass
            assert round(p.P) == round(self.P), "BatchPaSR does not support particles with different pressures (yet)"
        self.HY = hy
        self.mass = m
        assert self.N <= self.N_MAX , f"N ({self.N}) > N_MAX ({self.N_MAX}); too many particles"
        self.timeHistory_list.append(self.out_state)

    def insert(self, particle):
        self.particle_list.append(particle)
        self += particle
        # self.update_state() # Is the cost of assigning N to a new value (in update_state()) higher than doing N += 1?

    @classmethod
    def react_helper(cls, inputs):
        HPY, dt = inputs
        particle = particles["gri30.xml"]
        particle(HPY)
        t1 = time.time()
        particle.react(dt)
        t2= time.time()
        print(f"len(particles) = {len(particles)}\nIn react helper now. PID = {os.getpid()}. Time taken to react for {dt} = {(t2-t1)/1e-3:.3f} ms\n")
        return particle.HPY

    def react(self, parallel=False, coalesce = True):
        if parallel:
            if self.pool == None:
                self.pool = multiprocessing.Pool(processes=NUM_PROCESS,
                                                initializer=init_process,
                                                initargs=(self.mech, self.P))
            current_states = [p.HPY for p in self.particle_list]

            print(f"In parallel block; len(current_states) = {len(current_states)}")
            helper_input = zip(current_states, itertools.repeat(self.dt))
            new_states = self.pool.map(LiME.react_helper, helper_input)
            [p(new_states[i]) for i,p in enumerate(self.particle_list)]
        else:
            [p.react(self.dt) for p in self.particle_list]

        # Optional coalescing of particles for performance
        if coalesce:
            self.gobble()

        self.time += self.dt
        self.update_state()

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
                if self._can_combine(p1, p2):
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
                    self.particle_timeHistory_list.append(self.particle_list[ind].get_time_history())
                    self.particle_timeHistory_info.append([len(self.particle_timeHistory_list[-1]), self.time, self.particle_list[ind].age])
                    del self.particle_list[ind] # delete particles
                    particles_to_delete -= 1 # the list is getting shorter, so adjust the index of particles to be deleted accordingly
            self.update_state() # WE DIDN'T HAVE THIS BEFORE, SO THE BATCHPASR DIDN'T KNOW TO UPDATE IT'S MASS/STATE
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

    def _can_combine(self, p1, p2, tol=1e-6):
        dHY = p2 - p1
        diffH = dHY[0]
        diffY = dHY[1:]
        machine_epsilon = np.finfo(np.float64).eps
        diffH_percent = (diffH/(p1.enthalpy_mass + machine_epsilon))**2
        diffY_percent = np.linalg.norm(np.divide(diffY, p1.Y + machine_epsilon))
        diffH_is_small = diffH_percent < tol
        diffY_is_small = diffY_percent < tol
        # diffY_is_small = np.linalg.norm(np.divide(diffY, p0.state[1:] + np.finfo(np.float64).eps)) < tol
        if (diffH_is_small and diffY_is_small):
            print(f"We are combining particles at t = {self.time/1e-3:.3f} ms! diffH_norm = {diffH_percent:.5E}, diffY_percent = {diffY_percent:.5E}, system mass = {self.mass:.3f} kg")

        return ( diffH_is_small and diffY_is_small )

#     def iem(cls, paticle_list)
    def _iem(self, particle_list, k):
        k_avg = k*self
        for p in particle_list:
            p(p * (k + 1) - k_avg)

    def mix(self, tau_mix): # note: make sure dt < tau_mix!
        # Constant k:
        k = -self.dt/tau_mix # note: actually, k = -0.5*C_phi*omega*dt, but since C_phi is usually 2, i canceled it out.
        k_avg = k*self
        [p(p * (k + 1) - k_avg) for p in self.particle_list] # setting p.state_new = p.state_old + k*p.state_old - k*avg_state
        self.update_state()

    def prep_entrainment(self, added_gas, total_mass_added, tau_ent, numParticles=10, method='constant', time_interval = None):
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
            tempParticle = Particle.from_gas(added_gas, particle_mass = entrainmentMass[i])
            self.inactive_particles.append(tempParticle)

    def entrain(self, current_time):
        # Adds the next inactive particle if the time is correct
        if self.entrainInd < len(self.inactive_particles) and current_time >= self.entrain_timer[self.entrainInd]:
            # Note: The code right now will wait until the time is at or past the defined checkpoints
            # It will not entrain exactly at the defined time if the time steps do not synchronize
            # However, this should be less of an issue if more particles and smaller time steps are used
            current_particle = self.inactive_particles[self.entrainInd]
            current_particle.react(current_time)  # Continue reacting until entrainment. Note: assumes current_time starts from 0
            self.insert(current_particle)
            self.entrainInd += 1


    def get_time_history(self, dataframe=True):
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
        if dataframe == True:
            df = pd.DataFrame(columns = self.column_names, data = self.timeHistory_array)
            df.set_index(['age'])
            return df

        return self.timeHistory_array

    def get_particleTimeHistory(self, dataframe=True):
        self.update_state()
        for i in range(0, len(self.particle_list)):
            self.particle_timeHistory_list.append(self.particle_list[i].get_time_history())
            self.particle_timeHistory_info.append([len(self.particle_timeHistory_list[-1]), self.time, self.particle_list[i].age])
        particle_timeHistory_array = np.vstack(self.particle_timeHistory_list)
        if dataframe == True:
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

        return particle_timeHistory_array, self.particle_timeHistory_info



if __name__ == "__main__":
    multiprocessing.set_start_method("spawn")
    gas = ct.Solution("gri30.xml");
    p1 = Particle.from_gas(gas);
    gas()
