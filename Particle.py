import cantera as ct 
import numpy as np 
class Particle(ct.Solution):
    """Class for particle in BatchPaSR.
    """
    gas_template = None
    
    @classmethod
    def fromGas(cls, gas, particle_mass = 1.0, chemistry = True):
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
        # if self.gas_template == None:
        #     self.gas_template = ct.Solution(gas.name + ".xml")
        # return cls(gas.state, particle_mass, gas.name + ".xml", chemistry)
        # pdb.set_trace()
        state_vec = gas.state
        return cls(infile=f"{gas.name}.xml", phaseid='', source=None, thermo=None, species=(), kinetics=None, reactions=(), state_vec=state_vec, particle_mass=particle_mass)
    
    # def __init__(self, state, particle_mass = 1.0, mech='gri30.xml', chemistry = True):
    def __init__(self, infile='', phaseid='', source=None, thermo=None, species=(), kinetics=None, reactions=(), particle_mass=1.0, state_vec=[], P=101325, **kwargs):
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

        self.column_names = ['age', 'T', 'MW', 'h', 'phi', 'mass'] + ["Y_" + sn for sn in self.species_names] + ["X_" + sn for sn in self.species_names]        
        self.mass = particle_mass
        # self.age = 0
        # # self.state = np.hstack((self.gas_template.enthalpy_mass, self.gas_template.Y))
        # self.timeHistory_list = [[self.age, self.T, self.mean_molecular_weight, self.enthalpy_mass, self.get_equivalence_ratio(), self.mass] + self.Y.tolist() + self.X.tolist()]
        # self.timeHistory_array = None
        # self.reac = ct.ConstPressureReactor(self)
        # self.reac.chemistry_enabled = True
        # self.net = ct.ReactorNet([reac])
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
                self.HPY = comp
            else:
                return NotImplemented
            return self.HPY
        else:
            print(self.report())
   
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
            return self.HPY + other.HPY
        elif isinstance(other, np.ndarray):
            return self.HPY + other
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
            # TODO: CHECK THIS!!!!
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
        self.gas_template.HPY = [self.state[0], self.P, self.state[1:]]
        # self.timeHistory_list.append([self.age, Particle.gas_template.T, Particle.gas_template.mean_molecular_weight, Particle.gas_template.enthalpy_mass, Particle.gas_template.get_equivalence_ratio(), self.mass] + Particle.gas_template.Y.tolist() + Particle.gas_template.X.tolist())        
        self.reac.syncState() # this will automatically call self.net.reinitialize()
        self.reac.volume = self.mass/Particle.gas_template.density
        self.net.advance(dt)
        self.age += dt
        #         self.timeHistory_list = [[self.age, Particle.gas_template.T, Particle.gas_template.mean_molecular_weight, Particle.gas_template.enthalpy_mass, Particle.gas_template.get_equivalence_ratio()] + Particle.gas_template.Y.tolist() + Particle.gas.X.tolist()]        

        self.timeHistory_list.append([self.age, self.gas_template.T, self.gas_template.mean_molecular_weight, self.gas_template.enthalpy_mass, self.gas_template.get_equivalence_ratio(), self.mass] + self.gas_template.Y.tolist() + self.gas_template.X.tolist())
        self.state = np.hstack((self.gas_template.enthalpy_mass, self.gas_template.Y))
        # self.timeHistory_list.append([self.age, Particle.gas_template.T, Particle.gas_template.mean_molecular_weight, Particle.gas_template.enthalpy_mass, Particle.gas_template.get_equivalence_ratio(), self.mass] + Particle.gas_template.Y.tolist() + Particle.gas_template.X.tolist())
        # self.state = np.hstack((Particle.gas_template.enthalpy_mass, Particle.gas_template.Y))

    def get_timeHistory(self, timeOffset=0, dataFrame=False):
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
        if timeOffset > 0:
            self.timeHistory_array[:,0] += timeOffset
        if dataFrame == True:
            df = pd.DataFrame(columns = self.column_names, data = self.timeHistory_array)
            df.set_index(['age'])
            return df
        
        return self.timeHistory_array

    # Custom Equivalence Ratio test
    def get_global_equivalence_ratio(self, oxidizers = [], ignore = []):
        self.gas_template.HPY = [self.state[0], self.P, self.state[1:]]
        if not oxidizers:  # Default behavior, find all possible oxidizers
            oxidizers = [s.name for s in self.gas_template.species()] # if
                        # all(y not in s.composition for y in ['C', 'H', 'S'])]
            alpha = 0
            mol_O = 0
            for k, s in enumerate(self.gas_template.species()):
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

    def __getstate__(self):
            # how to get the state data out of a Particle instance
            state = self.__dict__.copy()
            return state

    def __setstate__(self, state):
            # rebuild a Particle instance from state
            self.__dict__.update(state)

