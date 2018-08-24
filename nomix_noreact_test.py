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
bp.prepEntrainment(added_gas = vit_reactor.thermo, total_mass_added = (mam+mfm), tau_ent = totalTime, numParticles=10, method='constant', time_interval = customInterval)

states = ct.SolutionArray(bp.mean_gas, extra=['t'])
enthalpy = []
for ti in times:
    states.append(bp.mean_gas.state, t=ti)
    enthalpy.append(bp.mean_gas.enthalpy_mass)
    
    bp.entrain(ti)
    bp.react()
#    bp.mix(tau_mix=0.1*milliseconds)

fig1, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
ax1.set_title('Temperature over Time in Second Stage')
ax2.set_title('Enthalpy of Second Stage over Time')
fig2, (ax3, ax4) = plt.subplots(nrows=2, ncols=1)
ax3.set_title('NO Concentration over time (15% O2 corr ppm)')
ax4.set_title('CO Concentration over time (15% O2 corr ppm)')
fig3, ax5 = plt.subplots()
fig3.suptitle('Mass of System')

for i in range(1):
    # df = bp.particle_list[i].get_timeHistory(dataFrame=True)
    # labelstr = 'Secondary Stream Particle' if i == 0 else 'Entrained Particle ' + str(i)
    NO_corr = correctNOx(states('NO').X, states('H2O').X, states('O2').X)
    CO_corr = correctNOx(states('CO').X, states('H2O').X, states('O2').X)
    ax1.plot(states.t, states.T)
    ax2.plot(states.t, enthalpy)
    ax3.plot(states.t, NO_corr)
    ax4.plot(states.t, CO_corr)
    df = bp.get_timeHistory()
    ax5.plot(df['age'], df['mass'])

for ax in [ax1, ax2, ax3, ax4, ax5]:
    ax.grid(True)
    ax.label_outer()
#    ax.legend()
# fig1.savefig('Entrained_Temperature_Plot.png')
# fig2.savefig('Entrained_NO_CO_Concentrations_Plot.png')
plt.show()