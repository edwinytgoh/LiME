from BatchPaSR import *
import matplotlib.pyplot as plt
# from BatchPaSR_3_Performance_Assessment import BatchPaSR

milliseconds = 0.001

[mfm, mam, mfs, mas] = solvePhi_airSplit(0.638, 0.4, 100, 1)
[vit_reactor, main_burner_DF] = runMainBurner(0.4, 19*milliseconds)

# vit_particle = Particle.from_gas(vit_reactor.thermo, particle_mass = entrainmentMass)
# g1 = ct.Solution('gri30.xml')

secondary_gas = ct.Solution('gri30.xml')
secondary_gas.TPX = 300, 25*ct.one_atm, {'CH4':1}
secondary_part = Particle.from_gas(secondary_gas, particle_mass=mfs + mas)

# Initial LiME
bp = PaSBR(particle_list=[secondary_part], dt=0.001*milliseconds, N_MAX=100)
states = ct.SolutionArray(bp.mean_gas, extra=['t'])
enthalpy = []
totalTime = 1*milliseconds
t = np.arange(0, totalTime, bp.dt)

# Constant entrainment rate for now (decide on more realistic with reference)
bp.prep_entrainment(added_gas = vit_reactor.thermo, total_mass_added = (mam + mfm), tau_ent = totalTime, numParticles=50, method='constant')
for i in range(0,t.size):
    bp.entrain(t[i])
    bp.react()
    bp.mix(tau_mix=0.01*milliseconds)
    states.append(bp.mean_gas.state, t=t[i])
    enthalpy.append(bp.mean_gas.enthalpy_mass)

# BR + MFC Code
[vit_reactor2, main_burner_DF] = runMainBurner(0.4, 19*milliseconds)

# vit_particle = Particle.from_gas(vit_reactor.thermo, particle_mass = entrainmentMass)
# g1 = ct.Solution('gri30.xml')

secondary_gas2 = ct.Solution('gri30.xml')
secondary_gas2.TPX = 300, 25*ct.one_atm, {'CH4':1}
sec_reactor = ct.ConstPressureReactor(secondary_gas2, volume = (mfs+mas)/secondary_gas2.density_mass)
sec_reactor.chemistry_enabled = True
# secondary_part = Particle.from_gas(secondary_gas, particle_mass=mfs+mas)

totalTime = 1*milliseconds
mfc = ct.MassFlowController(vit_reactor, sec_reactor, mdot = (mfm+mam)/totalTime)

reactorNet = ct.ReactorNet([vit_reactor, sec_reactor])
states2 = ct.SolutionArray(secondary_gas2, extra=['t'])
enthalpy2 = []

tnow = 0
dt = 0.001*milliseconds
for i in range(1000):
    reactorNet.advance(tnow)
    states2.append(sec_reactor.thermo.state, t=tnow)
    enthalpy2.append(sec_reactor.thermo.enthalpy_mass)
    tnow += dt

fig1, (ax1, ax4) = plt.subplots(nrows=2, ncols=1)
ax1.set_title('Average Temperature over Time')
ax4.set_title('Specific Enthalpy over Time')
plt.xlabel('Time (seconds)')
fig2, (ax2, ax3) = plt.subplots(nrows=2, ncols=1)
ax2.set_title('Average Gas NO Concentration over time (uncorrected ppm)')
ax3.set_title('Average Gas CO Concentration over time (uncorrected ppm)')
plt.xlabel('Time (seconds)')

#for i in range(1):
    #df = bp.particle_list[i].get_time_history(dataframe=True)
    #labelstr = 'Secondary Stream Particle' if i == 0 else 'Entrained Particle ' + str(i)
    #NO_corr = correctNOx(df['X_NO'], df['X_H2O'], df['X_O2'])
    #CO_corr = correctNOx(df['X_CO'], df['X_H2O'], df['X_O2'])
    #ax1.plot(df['age'], df['T'], label=labelstr)
    #ax2.plot(df['age'], NO_corr, label=labelstr)
    #ax3.plot(df['age'], CO_corr, label=labelstr)
    #ax4.plot(df['age'], df['h'], label=labelstr)
ax1.plot(states.t, states.T, label='Entrained Case')
ax2.plot(states.t, states('NO').X*1e6, label='Entrained Case')
ax3.plot(states.t, states('CO').X*1e6, label='Entrained Case')
ax4.plot(states.t, enthalpy, label='Entrained Case')

ax1.plot(states2.t, states2.T, label = 'BR + MFC')
ax2.plot(states2.t, states2('NO').X*1e6, label = 'BR + MFC')
ax3.plot(states2.t, states2('CO').X*1e6, label = 'BR + MFC')
ax4.plot(states2.t, enthalpy2, label = 'BR + MFC')

for ax in [ax1, ax2, ax3, ax4]:
    ax.grid(True)
    ax.label_outer()
    ax.legend()
fig1.savefig('Entrained_Temperature_Enthalpy.png')
fig2.savefig('Entrained_NO_CO_Concentrations.png')
# plt.show()