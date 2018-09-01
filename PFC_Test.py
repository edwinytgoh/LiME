from BatchPaSR import *
import matplotlib.pyplot as plt
# from BatchPaSR_3_Performance_Assessment import BatchPaSR

milliseconds = 0.001

[mfm, mam, mfs, mas] = solvePhi_airSplit(0.638, 0.4, 100, 1)
[vit_reactor, main_burner_DF] = runMainBurner(0.4, 19*milliseconds)

secondary_gas = ct.Solution('gri30.xml')
secondary_gas.TPX = 300, 25*ct.one_atm, {'CH4':1}
secondary_part = Particle.fromGas(secondary_gas, particle_mass=mfs+mas)

# Initialize PaSBR 
temp_gas = ct.Solution('gri30.xml')
temp_gas.TPX = 300, 25*ct.one_atm, {'AR':1} # I need some sort of filler for the PaSBR to work
temp_part = Particle.fromGas(temp_gas, particle_mass=1e-6)
bp = PaSBR(particle_list=[temp_part], dt=0.001*milliseconds, N_MAX=100)
states = ct.SolutionArray(bp.mean_gas, extra=['t'])
enthalpy = []
totalTime = 2*milliseconds
t = np.arange(0, totalTime, bp.dt)

# Constant entrainment rate for now (decide on more realistic with reference)
mdot_vit = (mfm+mam)/(1*milliseconds)
mdot_sec = (mfm+mam)/(1*milliseconds)
pfc_vit = ParticleFlowController(bp, vit_reactor.thermo, mfm+mam, totalTime/20, method=str(mdot_vit))
pfc_sec = ParticleFlowController(bp, secondary_gas, mfs+mas, totalTime/20, method=str(mdot_sec))
for i in range(0,t.size):
    pfc_vit.entrain(t[i])
    pfc_sec.entrain(t[i])
    bp.react()
    bp.mix(tau_mix=0.01*milliseconds)
    states.append(bp.mean_gas.state, t=t[i])
    enthalpy.append(bp.mean_gas.enthalpy_mass)

# BR + MFC Code
[vit_reactor2, main_burner_DF] = runMainBurner(0.4, 19*milliseconds)

secondary_gas2 = ct.Solution('gri30.xml')
secondary_gas2.TPX = 300, 25*ct.one_atm, {'CH4':1}
sec_reactor = ct.ConstPressureReactor(secondary_gas2, volume = (mfs+mas)/secondary_gas2.density_mass)

interface_gas2 = ct.Solution('gri30.xml')
interface_reactor = ct.ConstPressureReactor(interface_gas2)

mfc_vit = ct.MassFlowController(vit_reactor, interface_reactor, mdot = mdot_vit)
mfc_sec = ct.MassFlowController(sec_reactor, interface_reactor, mdot = mdot_sec)    # Using same flow rate

reactorNet = ct.ReactorNet([vit_reactor, sec_reactor, interface_reactor])
states2 = ct.SolutionArray(interface_gas2, extra=['t'])
enthalpy2 = []

tnow = 0.0
for i in range(0,t.size):
#    if tnow*mdot_sec > (mfs+mas):   # mfc_sec needs to stop when it runs out of mass
#        mfc_sec.set_mass_flow_rate(0.0)
    reactorNet.advance(tnow)
    states2.append(interface_reactor.thermo.state, t=tnow)
    enthalpy2.append(interface_reactor.thermo.enthalpy_mass)
    tnow = t[i]

# Plotting Settings
fig1, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
ax1.set_title('Average Temperature over Time')
ax2.set_title('Specific Enthalpy over Time')
plt.xlabel('Time (seconds)')
fig2, (ax3, ax4) = plt.subplots(nrows=2, ncols=1)
ax3.set_title('Average Gas NO Concentration over time (uncorrected ppm)')
ax4.set_title('Average Gas CO Concentration over time (uncorrected ppm)')
plt.xlabel('Time (seconds)')

ax1.plot(states.t, states.T, label='Entrained Case')
ax2.plot(states.t, enthalpy, label='Entrained Case')
ax3.plot(states.t, states('NO').X*1e6, label='Entrained Case')
ax4.plot(states.t, states('CO').X*1e6, label='Entrained Case')

ax1.plot(states2.t, states2.T, label = 'BR + MFC')
ax2.plot(states2.t, enthalpy2, label = 'BR + MFC')
ax3.plot(states2.t, states2('NO').X*1e6, label = 'BR + MFC')
ax4.plot(states2.t, states2('CO').X*1e6, label = 'BR + MFC')

for ax in [ax1, ax2, ax3, ax4]:
    ax.grid(True)
    ax.label_outer()
    ax.legend()
fig1.savefig('Mult_Entrained_Temperature_Enthalpy.png')
fig2.savefig('Mult_Entrained_NO_CO_Concentrations.png')
plt.show()