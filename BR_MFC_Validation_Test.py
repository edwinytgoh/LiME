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
# secondary_part = Particle.fromGas(secondary_gas, particle_mass=mfs+mas)

totalTime = 1*milliseconds
mfc = ct.MassFlowController(vit_reactor, sec_reactor, mdot = (mfm+mam)/totalTime)

reactorNet = ct.ReactorNet([vit_reactor, sec_reactor])
states = ct.SolutionArray(secondary_gas, extra=['t'])

tnow = 0
dt = totalTime/100
for i in range(100):
    reactorNet.advance(tnow)
    states.append(sec_reactor.thermo.state, t=tnow)
    tnow += dt

fig1, ax1 = plt.subplots()
fig1.suptitle('Particle Temperature over Time in Second Stage')
fig2, (ax2, ax3) = plt.subplots(nrows=2, ncols=1)
ax2.set_title('NO Concentration over time (15% O2 corr ppm)')
ax3.set_title('CO Concentration over time (15% O2 corr ppm)')

NO_corr = correctNOx(states('NO').X, states('H2O').X, states('O2').X)
CO_corr = correctNOx(states('CO').X, states('H2O').X, states('O2').X)
ax1.plot(states.t, states.T)
ax2.plot(states.t, NO_corr)
ax3.plot(states.t, CO_corr)

for ax in [ax1, ax2, ax3]:
    ax.grid(True)
    ax.label_outer()
fig1.savefig('BR_MFC_Temperature_Plot.png')
fig2.savefig('BR_MFC_NO_CO_Concentrations_Plot.png')