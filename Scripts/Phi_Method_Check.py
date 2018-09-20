import sys
sys.path.insert(0, "../")
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
vit_part = Particle.fromGas(vit_reactor.thermo, particle_mass=mfm+mam)
vit_gas = ct.Solution('gri30.xml')

# Initial PaSBR 
bp = PaSBR(particle_list=[secondary_part, vit_part], dt=0.001*milliseconds, N_MAX=100)
states = ct.SolutionArray(bp.mean_gas, extra=['t'])
phi_global = []
phi_unburned = []
phi_global_vit = []
phi_unburned_vit = []
enthalpy = []
totalTime = 1*milliseconds
t = np.arange(0, totalTime, bp.dt)

# Constant entrainment rate for now (decide on more realistic with reference)
# bp.prepEntrainment(added_gas = vit_reactor.thermo, total_mass_added = (mam+mfm), tau_ent = totalTime, numParticles=50, method='constant')
for i in range(0,t.size):
#     bp.entrain(t[i])
    bp.react()
    bp.mix(tau_mix=0.01*milliseconds)
    states.append(bp.mean_gas.state, t=t[i])
    phi_global.append(secondary_part.get_global_equivalence_ratio())
    secondary_gas.HPY = [secondary_part.state[0], secondary_part.P, secondary_part.state[1:]]
    phi_unburned.append(secondary_gas.get_equivalence_ratio())
    phi_global_vit.append(vit_part.get_global_equivalence_ratio())
    vit_gas.HPY = [vit_part.state[0], vit_part.P, vit_part.state[1:]]
    phi_unburned_vit.append(vit_gas.get_equivalence_ratio())
    enthalpy.append(bp.mean_gas.enthalpy_mass)



fig1, (ax1, ax4) = plt.subplots(nrows=2, ncols=1)
ax1.set_title('Average Temperature over Time')
ax4.set_title('Specific Enthalpy over Time')
plt.xlabel('Time (ms)')
fig2, (ax2, ax3) = plt.subplots(nrows=2, ncols=1)
ax2.set_title('Average Gas NO Concentration over time (uncorrected ppm)')
ax3.set_title('Average Gas CO Concentration over time (uncorrected ppm)')
plt.xlabel('Time (ms)')
fig3, (ax5) = plt.subplots(nrows=1, ncols=1)
ax5.set_title('Phi of Fuel Particle')
# ax6.set_title('Phi Unburned of Fuel Particle')
plt.xlabel('Time (ms)')
fig4, (ax6) = plt.subplots(nrows=1, ncols=1)
ax6.set_title('Phi of Air Particle')
plt.xlabel('Time (ms)')


ax1.plot(states.t/milliseconds, states.T, label='Entrained Case')
ax2.plot(states.t/milliseconds, states('NO').X*1e6, label='Entrained Case')
ax3.plot(states.t/milliseconds, states('CO').X*1e6, label='Entrained Case')
ax4.plot(states.t/milliseconds, enthalpy, label='Entrained Case')
ax5.plot(states.t/milliseconds, phi_global, label='Phi Global')
ax5.plot(states.t/milliseconds, phi_unburned, label='Phi Unburned')
ax5.set_ylim((0, 5))
ax6.plot(states.t/milliseconds, phi_global_vit, label='Phi Global')
ax6.plot(states.t/milliseconds, phi_unburned_vit, label='Phi Unburned')
ax6.set_ylim((0, 5))


for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
    ax.grid(True)
    ax.label_outer()
    ax.legend()
# fig1.savefig('Entrained_Temperature_Enthalpy.png')
# fig2.savefig('Entrained_NO_CO_Concentrations.png')
plt.show()