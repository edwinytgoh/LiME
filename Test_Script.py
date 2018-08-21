from BatchPaSR import *
import matplotlib.pyplot as plt
# from BatchPaSR_3_Performance_Assessment import BatchPaSR

milliseconds = 0.001;

[mfm, mam, mfs, mas] = solvePhi_airSplit(0.638, 0.4, 100, 1)
[vit_reactor, main_burner_DF] = runMainBurner(0.4, 19*milliseconds)

vit_particle = Particle.fromGas(vit_reactor.thermo, particle_mass = (mfm + mam));
# g1 = ct.Solution('gri30.xml'); 

secondary_gas = ct.Solution('gri30.xml'); 
secondary_gas.TPX = 300, 25*ct.one_atm, {'CH4':1}
secondary_part = Particle.fromGas(secondary_gas, particle_mass=mfs+mas)

bp = PaSBR(particle_list=[vit_particle, secondary_part], dt=0.001*milliseconds)
t = np.arange(0, 1*milliseconds, bp.dt)
for i in range(0,t.size):
    bp.react();
    bp.mix(tau_mix=0.1*milliseconds)

df1 = bp.particle_list[0].get_timeHistory(dataFrame=True)
df2 = bp.particle_list[1].get_timeHistory(dataFrame=True)

plt.plot(df1['age'], df2['T'])
plt.plot(df1['age'], df1['T'])
plt.show()