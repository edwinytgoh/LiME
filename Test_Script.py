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
secondary_part = Particle.fromGas(secondary_gas, particle_mass=mfs+mas)

# Initial PaSBR 
bp = PaSBR(particle_list=[secondary_part], dt=0.001*milliseconds, N_MAX=50)
totalTime = 1*milliseconds
t = np.arange(0, totalTime, bp.dt)

# Constant entrainment rate for now (decide on more realistic with reference)
bp.prepEntrainment(added_gas = vit_reactor.thermo, total_mass_added = (mam+mfm), tau_ent = totalTime, numParticles=30, method='constant')
for i in range(0,t.size):
    bp.entrain(t[i])
    bp.react()
    bp.mix(tau_mix=0.1*milliseconds)

fig1, ax1 = plt.subplots()
fig1.suptitle('Particle Temperature over Time in Second Stage')
fig2, (ax2, ax3) = plt.subplots(nrows=2, ncols=1)
ax2.set_title('NO Concentration over time (15% O2 corr ppm)')
ax3.set_title('CO Concentration over time (15% O2 corr ppm)')

for i in range(1):
    df = bp.particle_list[i].get_timeHistory(dataFrame=True)
    labelstr = 'Secondary Stream Particle' if i == 0 else 'Entrained Particle ' + str(i)
    NO_corr = correctNOx(df['X_NO'], df['X_H2O'], df['X_O2'])
    CO_corr = correctNOx(df['X_CO'], df['X_H2O'], df['X_O2'])
    ax1.plot(df['age'], df['T'], label=labelstr)
    ax2.plot(df['age'], NO_corr, label=labelstr)
    ax3.plot(df['age'], CO_corr, label=labelstr)

for ax in [ax1, ax2, ax3]:
    ax.grid(True)
    ax.label_outer()
    ax.legend()
fig1.savefig('Entrained_Temperature_Plot.png')
fig2.savefig('Entrained_NO_CO_Concentrations_Plot.png')