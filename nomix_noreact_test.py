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
secondary_part = Particle.fromGas(secondary_gas, particle_mass=mfs+mas, chemistry = False)

# Initial PaSBR 
bp = PaSBR(particle_list=[secondary_part], dt=0.001*milliseconds, N_MAX=50, chemistry = False)
totalTime = 1*milliseconds
times = np.arange(0, totalTime, bp.dt)

# Constant entrainment rate for now (decide on more realistic with reference)
customInterval = np.array([0., 1., 1., 1., 2., 2., 2., 3., 3., 4.])
customInterval *= 0.95*totalTime/np.sum(customInterval) # 0.95 so that last particle enters before end
customInterval = customInterval.cumsum()
bp.prepEntrainment(added_gas = vit_reactor.thermo, total_mass_added = (mam+mfm), tau_ent = totalTime, numParticles=10, method='constant')   #, time_interval = customInterval)

states = ct.SolutionArray(bp.mean_gas, extra=['t'])
enthalpy = []
for ti in times:
    states.append(bp.mean_gas.state, t=ti)
    enthalpy.append(bp.mean_gas.enthalpy_mass)
    
    bp.entrain(ti)
    bp.react()
#    bp.mix(tau_mix=0.1*milliseconds)


# Separate BR + MFC run to plot on top of each other

[vit_reactor2, main_burner_DF] = runMainBurner(0.4, 19*milliseconds)

# vit_particle = Particle.fromGas(vit_reactor.thermo, particle_mass = entrainmentMass)
# g1 = ct.Solution('gri30.xml')

secondary_gas2 = ct.Solution('gri30.xml')
secondary_gas2.TPX = 300, 25*ct.one_atm, {'CH4':1}
sec_reactor = ct.ConstPressureReactor(secondary_gas2, volume = (mfs+mas)/secondary_gas2.density_mass)
sec_reactor.chemistry_enabled = False
# secondary_part = Particle.fromGas(secondary_gas, particle_mass=mfs+mas)

totalTime = 1*milliseconds
mfc = ct.MassFlowController(vit_reactor, sec_reactor, mdot = (mfm+mam)/totalTime)

reactorNet = ct.ReactorNet([vit_reactor, sec_reactor])
states2 = ct.SolutionArray(secondary_gas2, extra=['t'])
enthalpy2 = []
mass = []

tnow = 0
dt = totalTime/100
for i in range(100):
    reactorNet.advance(tnow)
    states2.append(sec_reactor.thermo.state, t=tnow)
    enthalpy2.append(sec_reactor.thermo.enthalpy_mass)
    mass.append(sec_reactor.mass)
    tnow += dt

fig1, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
ax1.set_title('Temperature over Time in Second Stage')
ax2.set_title('Enthalpy of Second Stage over Time')
fig2, (ax3, ax4) = plt.subplots(nrows=2, ncols=1)
ax3.set_title('NO Concentration over time (15% O2 corr ppm)')
ax4.set_title('CO Concentration over time (15% O2 corr ppm)')
fig3, ax5 = plt.subplots()
fig3.suptitle('Mass of System')

NO_corr = correctNOx(states('NO').X, states('H2O').X, states('O2').X)
CO_corr = correctNOx(states('CO').X, states('H2O').X, states('O2').X)
ax1.plot(states.t, states.T, label = 'Entrained Case')
ax2.plot(states.t, enthalpy, label = 'Entrained Case')
ax3.plot(states.t, NO_corr, label = 'Entrained Case')
ax4.plot(states.t, CO_corr, label = 'Entrained Case')
df = bp.get_timeHistory()
ax5.plot(df['age'], df['mass'], label = 'Entrained Case')

NO_corr2 = correctNOx(states2('NO').X, states2('H2O').X, states2('O2').X)
CO_corr2 = correctNOx(states2('CO').X, states2('H2O').X, states2('O2').X)
ax1.plot(states2.t, states2.T, label = 'BR + MFC')
ax2.plot(states2.t, enthalpy2, label = 'BR + MFC')
ax3.plot(states2.t, NO_corr2, label = 'BR + MFC')
ax4.plot(states2.t, CO_corr2, label = 'BR + MFC')
ax5.plot(states2.t, mass, label = 'BR + MFC')

for ax in [ax1, ax2, ax3, ax4, ax5]:
    ax.grid(True)
    ax.label_outer()
    ax.legend()
#fig1.savefig('BR_MFC_Temperature_Plot.png')
#fig2.savefig('BR_MFC_NO_CO_Concentrations_Plot.png')
plt.show()