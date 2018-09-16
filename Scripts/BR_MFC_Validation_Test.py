from BatchPaSR import *
import matplotlib.pyplot as plt
# from BatchPaSR_3_Performance_Assessment import BatchPaSR

milliseconds = 0.001

[mfm, mam, mfs, mas] = solvePhi_airSplit(0.638, 0.4, 100, 1)
[vit_reactor, main_burner_DF] = runMainBurner(0.4, 19*milliseconds)

# vit_particle = Particle.fromGas(vit_reactor.thermo, particle_mass = entrainmentMass)
# g1 = ct.Solution('gri30.xml')

secondary_gas = ct.Solution('gri30.xml')
secondary_gas.TPX = 300, 25*ct.one_atm, {'CH4':1}
sec_reactor = ct.ConstPressureReactor(secondary_gas, volume = (mfs+mas)/secondary_gas.density_mass)
sec_reactor.chemistry_enabled = False
# secondary_part = Particle.fromGas(secondary_gas, particle_mass=mfs+mas)

totalTime = 1*milliseconds
mfc = ct.MassFlowController(vit_reactor, sec_reactor, mdot = (mfm+mam)/totalTime)

reactorNet = ct.ReactorNet([vit_reactor, sec_reactor])
states = ct.SolutionArray(secondary_gas, extra=['t'])
enthalpy = []
mass = []

tnow = 0
dt = totalTime/100
for i in range(100):
    reactorNet.advance(tnow)
    states.append(sec_reactor.thermo.state, t=tnow)
    enthalpy.append(sec_reactor.thermo.enthalpy_mass)
    mass.append(sec_reactor.mass)
    tnow += dt

fig1, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
ax1.set_title('Particle Temperature over Time in Second Stage')
ax2.set_title('Total Enthalpy of Reactor over Time')
fig2, (ax3, ax4) = plt.subplots(nrows=2, ncols=1)
ax3.set_title('NO Concentration over time (15% O2 corr ppm)')
ax4.set_title('CO Concentration over time (15% O2 corr ppm)')
fig3, ax5 = plt.subplots()
fig3.suptitle('Mass of System')

NO_corr = correctNOx(states('NO').X, states('H2O').X, states('O2').X)
CO_corr = correctNOx(states('CO').X, states('H2O').X, states('O2').X)
ax1.plot(states.t, states.T)
ax2.plot(states.t, enthalpy)
ax3.plot(states.t, NO_corr)
ax4.plot(states.t, CO_corr)
ax5.plot(states.t, mass)

for ax in [ax1, ax2, ax3, ax4, ax5]:
    ax.grid(True)
    ax.label_outer()
#fig1.savefig('BR_MFC_Temperature_Plot.png')
#fig2.savefig('BR_MFC_NO_CO_Concentrations_Plot.png')
plt.show()