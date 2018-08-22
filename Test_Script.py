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
bp = PaSBR(particle_list=[secondary_part], dt=0.001*milliseconds, N_MAX=11)
totalTime = 1*milliseconds
t = np.arange(0, totalTime, bp.dt)
# Defining entrainment here
# Using constant entrainment time step
# Constant entrainment rate for now (decide on more realistic with reference)
numParticles = 10
entrainmentMass = (mfm + mam)/numParticles
entrainTime = np.arange(0, totalTime, totalTime/numParticles)
inactiveGas = []
for i in range(0, numParticles):
    tempParticle = Particle.fromGas(vit_reactor.thermo, particle_mass = entrainmentMass)
    tempParticle.react(entrainTime[i])
    inactiveGas.append(tempParticle)
entrainInd = 0

for i in range(0,t.size):
    if entrainInd < numParticles and t[i] >= entrainTime[entrainInd]:
        bp.entrain(inactiveGas, entrainInd)
        entrainInd += 1
    bp.react()
    bp.mix(tau_mix=0.1*milliseconds)

df1 = bp.particle_list[0].get_timeHistory(dataFrame=True)
df2 = bp.particle_list[1].get_timeHistory(dataFrame=True)

ax = plt.subplot(111)
ax.plot(df2['age'], df2['T'], label='First Entrained Particle')
ax.plot(df1['age'], df1['T'], label='Secondary Stream Particle')
ax.legend()
plt.show()