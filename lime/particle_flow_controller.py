import cantera as ct
import numpy as np

from lime.particle import Particle


class ParticleFlowController(object):
    def __init__(self, reactor, gas, total_mass, timestep, method):
        """
        Parameters
        ----------
        reactor : `BatchPaSR`
            The reactor that this PFC acts as an inlet for

        gas : `cantera.Solution`
            Thermochemical state of gas to be entrained

        total_mass : `float`
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
        self.mass = total_mass
        self.dt = timestep
        self.mdotexp = method
        # self.nextstep = timestep/2   # Which time to add the next particle (using trapezoidal sum)
        self.nextstep = 0

    def can_entrain(self, t):
        # Check if we can even entrain at this time
        return (self.mass > 0) and (self.nextstep <= t)  # EDWIN COMMENT: Is nextstep <= t necessary?

    def entrain(self, t):
        # Will add a particle to the reactor if entrainment is possible
        if self.can_entrain(t):
            # Calculate the mass of added particle
            if callable(self.mdotexp):
                m_dot1 = self.mdotexp(self.nextstep - self.dt / 2)
                m_dot2 = self.mdotexp(self.nextstep + self.dt / 2)
            else:
                m_dot1 = eval(self.mdotexp.replace('t', str(self.nextstep - self.dt / 2)))
                m_dot2 = eval(self.mdotexp.replace('t', str(self.nextstep + self.dt / 2)))
            self.nextstep += self.dt
            # trapezoidal rule; if remaining self.mass < mass to be entrained then entrain self.mass
            entrained_mass = min(self.mass, (
                    m_dot1 + m_dot2) / 2 * self.dt)
            # Keep track of remaining mass

            self.mass -= entrained_mass
            self.rn.advance(t)
            particle = Particle.from_gas(self.gas, particle_mass=entrained_mass)
            # print(f"Entraining particle of mass {entrained_mass:.2f} kg into LiME");
            # Add it into the reactor
            self.bp.insert(particle)
            return self.mass, particle()
        else:
            return self.mass, np.hstack((self.gas.enthalpy_mass, self.gas.Y))