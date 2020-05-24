import cantera as ct 
import numpy as np 
import pandas as pd 
import pdb
class Particle(ct.Solution):
    """Class for particle in BatchPaSR.
    """
    gas_template = None

    @classmethod
    def fromGas(cls, gas, mech="", particle_mass = 1.0, rates = True):
        """Initialize particle object with thermochemical state.
        
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
        return cls(infile=mech, phaseid='', source=None, thermo=None, species=(), kinetics=None, reactions=(), state_vec=state_vec, particle_mass=particle_mass, P=gas.P, rates=rates)

    @classmethod
    def fromReactor(cls, reactor, mech="", particle_mass = 0, rates = True):
        """Initialize particle object with reactor reference.
        
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
        reactor.syncState()
        state_vec = reactor.thermo.state
        if mech == "":
            if isinstance(reactor.thermo, cls):
                mech = reactor.thermo.mech
            else:
                mech = f"{reactor.thermo.name}.xml"
        if particle_mass == 0: 
            particle_mass = reactor.mass
        else: 
            particle_mass = particle_mass
            reactor.volume = particle_mass/reactor.thermo.density
        c = cls(infile=mech, phaseid='', source=None, thermo=None, species=(), kinetics=None, reactions=(), state_vec=state_vec, particle_mass=particle_mass, P=reactor.thermo.P, rates=rates)                
        c.reactor = reactor
        return c
    

    # def __init__(self, state, particle_mass = 1.0, mech='gri30.xml', chemistry = True):
    def __init__(self, infile='', phaseid='', source=None, thermo=None, species=(), kinetics=None, reactions=(), particle_mass=1.0, state_vec=[], P=101325, rates=True, **kwargs):
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
        super().__init__()
        if len(state_vec) > 0:
            self.state = state_vec
        # if Particle.gas_template == None:
        #     Particle.gas_template = ct.Solution(mech, name=mech[0:-4])        
        self.mech = infile
        self.TPX = self.T, P, self.X
        # self.P = self.state[0]*self.state[1]/self.mean_molecular_weight*ct.gas_constant # P = rho * R * T; 

        self.mass = particle_mass
        self.age = 0
        self.reactor = ct.ConstPressureReactor(self, volume=self.mass/self.density)
        self.timeHistory_list = [self.outState]
        self.rate_list = [np.hstack([self.age, self.net_production_rates])]
        self.net = ct.ReactorNet([self.reactor])
    
    # @property
    # def mass(self):
    #     return self.volume * self.density_mass
    # @mass.setter
    # def mass(self, new_mass):
    #     new_vol = new_mass/self.density_mass
    #     self.volume
    

    @property
    def column_names(self):
        return np.hstack(['age', 'T', 'MW', 'h', 'phi', 'mass', 'density_mole', 'density_mass', 'NO_production_rate'] + ["Y_" + sn for sn in self.species_names] + ["X_" + sn for sn in self.species_names])             

    @property
    def outState(self):
        self.reactor.syncState()
        self.state = self.reactor.thermo.state # only useful if particle was created using Particle.fromReactor
        return np.vstack([[self.age, self.T, self.mean_molecular_weight, self.enthalpy_mass, self.get_equivalence_ratio(), self.mass, self.density_mole, self.density_mass, self.net_production_rates[self.species_index('NO')]] 
        + self.Y.tolist() + self.X.tolist()])

    @property
    def reaction_rates(self):
        return np.hstack([self.age, self.net_production_rates])

    @property
    def HY(self):
        return np.insert(self.HPY[2], 0, self.HPY[0])
        
    @HY.setter
    def HY(self, vec):
        if (len(vec) == len(self.HY)):
            self.HPY = vec[0], self.P, vec[1:]
        elif (len(vec) == len(self.HPY)):
            self.HPY = vec

    @property
    def volume(self):
        return self.mass/self.density
    
    @volume.setter
    def volume(self, new_vol):
        self.reactor.volume = new_vol

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
            if isinstance(comp, ct.Solution):
                # h = comp.enthalpy_mass
                # Y = comp.Y
                self.HPY = comp.HPY
            elif isinstance(comp, np.ndarray):
                if len(comp) == len(self.HY):
                    self.HY = comp
            elif isinstance(comp, tuple) and len(comp) == len(self.HPY):
                self.HPY = comp # ignore 
            else:
                return NotImplemented
            return self.HPY
        else:
            print(self.report())
            pass
   
    def __add__(self, other):
        """Add values to state of particle without changing the state of either particle.

        Parameters
        ----------
        other : `Particle`, `numpy.array`, `int`, `float`
            Thermochemical state (ct.Solution.HPY) to add to current state.

        Returns
        -------
        comp : numpy.array
            Thermochemical composition of particle (enthalpy + mass fractions).

        """
        if isinstance(other, Particle):
            return self.HY + other.HY
        elif isinstance(other, np.ndarray):
            if (len(other) == len(self.HY)): 
                return self.HY + other
        elif isinstance(other,tuple) and (len(other) == len(self.HPY)):
            return (self.HPY[i] + other[i] for i in range(0, len(self.HPY)))
        # elif isinstance(other, (int, float)):
        #     return self.state + other
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
        return self.__add__(other)

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
            return self.HY - other.HY
        elif isinstance(other, np.ndarray):
            if (len(other) == len(self.HY)): 
                return self.HY - other
        elif isinstance(other, tuple) and (len(other) == len(self.HPY)):
                return (self.HPY[i] - other[i] for i in range(0, len(self.HPY)))
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
        return -1 * (self.__sub__(other))

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
            return self.HY * other
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
            return self.__mul__(other)
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
            # TODO: CHECK THIS!!!!
            # h = (self.mass*self.enthalpy_mass + other.mass*other.enthalpy_mass)/(self.mass + other.mass)
            # Y = (self.mass*self.Y + other.mass*other.Y)/(self.mass + other.mass)
            total_mass = self.mass + other.mass
            self.HY = (self.mass/total_mass)*self.HY + (other.mass/total_mass)*other.HY
            self.mass += other.mass
        elif isinstance(other, np.ndarray):
            assert len(other) == len(self.HY) , "Please ensure that input array has the same length as current particle state array ({0.d})".format(len(self.state))
            # h = self.state[0] + other[0]
            # Y = self.state[1:] + other[1:]
            self.HY = other
            # NOTE: DOES IT MAKE SENSE TO NOT HAVE MASS? 
        else:
            return NotImplemented
        # self.timeHistory_list.append(self.outState)
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
            # TODO: CHECK THIS!!!!
            # h = (self.mass*self.enthalpy_mass + other.mass*other.enthalpy_mass)/(self.mass + other.mass)
            # Y = (self.mass*self.Y + other.mass*other.Y)/(self.mass + other.mass)
            total_mass = self.mass - other.mass
            self.HY = (self.mass/total_mass)*self.HY - (other.mass/total_mass)*other.HY
            self.mass -= other.mass
        elif isinstance(other, np.ndarray):
            assert len(other) == len(self.HY) , "Please ensure that input array has the same length as current particle state array ({0.d})".format(len(self.state))
            # h = self.state[0] + other[0]
            # Y = self.state[1:] + other[1:]
            self.HY = other
            # !NOTE: DOES IT MAKE SENSE TO NOT HAVE MASS? 
        else:
            return NotImplemented
        # self.timeHistory_list.append(self.outState)            
        return self        
        # if isinstance(other, Particle):
        #     h = self.state[0] - other.state[0]
        #     Y = self.state[1:] - other.state[1:]
        # elif isinstance(other, np.ndarray):
        #     assert len(other) == len(self.state) , "Please ensure that input array has the same length as current particle state array ({0.d})".format(len(self.state))
        #     h = self.state[0] - other[0]
        #     Y = self.state[1:] - other[1:]
        # elif isinstance(other, (int, float)):
        #     h = self.state[0] - other
        #     Y = self.state[1:] - other
        # else:
        #     return NotImplemented
        # self.state[0] = h
        # self.state[1:] = Y
        # return self

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
            self.HY *= other
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
        # pdb.set_trace()
        # self.net.reinitialize()
        # self.reactor.syncState()
        self.net.set_initial_time(self.age)
        self.age += dt
        self.net.advance(self.age)
        if self.reactor.mass >= self.mass:
            self.mass = self.reactor.mass
        self.reactor.volume = self.mass/self.density
        self.reactor.syncState() # this will automatically call self.net.reinitialize()
        self.timeHistory_list.append(self.outState)
        self.rate_list.append(self.reaction_rates)
        # TODO: Make sure Particle updates state 
        
    def get_timeHistory(self, timeOffset=0, dataFrame=False, deleteFirstElem=True):
        """Obtain particle's history. 
        
        Parameters
        ----------
        None
        
        Returns
        -------
        timeHistory_array : `numpy.array` 
            The array containing particle property time traces 
        """
        if deleteFirstElem:
            timeHistory_array = np.vstack(self.timeHistory_list[1:])
        else:
            timeHistory_array = np.vstack(self.timeHistory_list)
        if timeOffset > 0:
            timeHistory_array[:,0] += timeOffset
        if dataFrame == True:
            df = pd.DataFrame(columns = self.column_names, data = timeHistory_array)
            df.set_index(['age'])
            return df
        
        return timeHistory_array

    def get_rateHistory(self, dataFrame=False, deleteFirstElem=True):
        """ Obtain particle's net production rate history. 

        Parameters 
        -----------
        dataFrame : bool 
            Boolean value indicating whether or not to output pandas DataFrame 
        
        deleteFirstElem : bool
            Boolean value indicating whether or not to delete first element at t = 0
        """

        if deleteFirstElem:
            rate_array = np.vstack(self.rate_list[1:])
        else:
            rate_array = np.vstack(self.rate_list)
        if dataFrame == True: 
            df = pd.DataFrame(columns = ['age'] + [f"{sp}_productionRate" for sp in self.species_names], data=rate_array)
            return df 
        return rate_array

    # Custom Equivalence Ratio test
    def get_global_equivalence_ratio(self, oxidizers = [], ignore = []):
        if not oxidizers:  # Default behavior, find all possible oxidizers
            oxidizers = [sn for sn in self.species_names] # if
                        # all(y not in s.composition for y in ['C', 'H', 'S'])]
            alpha = 0
            mol_O = 0
            for k, s in enumerate(self.gas_template.species):
                if s.name in ignore:
                    continue
                else: # elif s.name in oxidizers:
                    mol_O += s.composition.get('O', 0) * self.gas_template.X[k]
                # else:
                    nC = s.composition.get('C', 0)
                    nH = s.composition.get('H', 0)
                    # nO = s.composition.get('O', 0)
                    nS = s.composition.get('S', 0)

                    alpha += (2*nC + nH/2 + 2*nS) * self.gas_template.X[k]

            if mol_O == 0:
                return float('inf')
            else:
                return alpha / mol_O

    # def __getstate__(self):
    #         # how to get the state data out of a Particle instance
    #         state = self.__dict__.copy()
    #         return state

    # def __setstate__(self, state):
    #         # rebuild a Particle instance from state
    #         self.__dict__.update(state)

