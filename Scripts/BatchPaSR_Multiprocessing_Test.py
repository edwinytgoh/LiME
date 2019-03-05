import sys
sys.path.insert(0, "../")
from CanteraTools import *
print("Hello World")
import BatchPaSR as bp
from Particle import Particle
import os
from tqdm import tqdm
# multiprocessing.set_start_method("spawn")
fuel = ct.Solution('gri30.xml')
fuel.TPX = 300, 25*ct.one_atm, {'CH4':1} # TODO: Come up with a more general way of doing secondary gas so that we can have both fuel and air 
air = ct.Solution('gri30.xml'); 
air.TPX = 650, 25*ct.one_atm, {'O2':0.21, 'N2':0.79}

def run_finite_everything(tau_mix, tau_ent_main, tau_ent_sec, 
                          phi_global = 0.635, phi_main = 0.3719, phi_jet = np.inf, 
                          tau_main = (20-0.158)*1e-3, tau_sec = 5*1e-3, dt = 0.01*1e-3, P = 25*ct.one_atm, mech="gri30.xml"):
    # [mfm, mam, mfs, mas] = solvePhi_airSplit(phi_global, phi_main, 100, airSplit)
    [mfm, mam, mfs, mas] = solveMass_PhiJet(phi_global, phi_main, phi_jet, mdotTotal=100)
    
    mass_main = mfm + mam
    mass_sec = mfs + mas
    mdot_main = mass_main/tau_ent_main
    mdot_sec = mass_sec/tau_ent_sec
    t = np.arange(0, tau_sec, dt)
    print(f"Length of t is {len(t)}; start = 0; stop = {tau_sec:.5f}; step = {dt:.6f}")
    # Calculate main burner:
    [vit_reactor, main_burner_DF] = runMainBurner(phi_main, tau_main, P=P)    

    # Setup PaSBR:
    secondary_gas = mix([fuel, air], [mfs, mas], P = P)
    init_gas = Particle(mech, particle_mass=0.001, P=P)
    init_gas.X = {'AR':1.0}
    
    pasbr = bp.PaSBR.fromGas(init_gas, N_MAX=5000)
    pfc_main = bp.ParticleFlowController(pasbr, vit_reactor.thermo, mass_main, dt, lambda t: mdot_main if t <= tau_ent_main else 0)
    pfc_sec = bp.ParticleFlowController(pasbr, secondary_gas, mass_sec, dt, lambda t: mdot_sec if t <= tau_ent_sec else 0)

    # Setup data management: 
    mean_gas = ct.Solution(mech)
    mean_gas.TPX = secondary_gas.TPX
    system_timeHistory = []
    num_particles_list = []
    # Simulation loop: 
    for t_now in tqdm(t):
        [remaining_main_mass, main_state] = pfc_main.entrain(t_now)
        [remaining_sec_mass, sec_state] = pfc_sec.entrain(t_now)
        total_mass = remaining_main_mass + remaining_sec_mass + pasbr.mass
#         system_state = (remaining_main_mass * main_state + remaining_sec_mass * sec_state + pasbr.mass * pasbr.state)/total_mass
#         mean_gas.HPY = system_state[0], mean_gas.P, system_state[1:]
#         system_timeHistory.append([t_now, total_mass, mean_gas.T, mean_gas.mean_molecular_weight, mean_gas.enthalpy_mass, mean_gas.get_equivalence_ratio()] + mean_gas.Y.tolist() + mean_gas.X.tolist())
        pasbr.react(parallel=True)
        pasbr.mix(tau_mix = tau_mix)
        num_particles_list.append(len(pasbr.particle_list))

    # Output data 
    entrainment_zone_df = pasbr.get_timeHistory()
#     system_timeHistory = np.vstack(system_timeHistory)
#     system_df = pd.DataFrame(columns=entrainment_zone_df.columns, data=system_timeHistory)
#     system_df['num_particles'] = np.array(num_particles_list)

    return entrainment_zone_df, system_df, pasbr

if __name__ == "__main__":
    multiprocessing.set_start_method("spawn")
    run_finite_everything(tau_mix=0.1,tau_ent_main=0.5,tau_ent_sec=0.05)
    gas()