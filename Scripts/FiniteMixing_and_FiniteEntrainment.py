import sys
sys.path.insert(0, "../")
from CanteraTools import *
print("Hello World")
import BatchPaSR as bp

fuel = ct.Solution('gri30.xml')
fuel.TPX = 300, 25*ct.one_atm, {'CH4':1} # TODO: Come up with a more general way of doing secondary gas so that we can have both fuel and air 
air = ct.Solution('gri30.xml'); 
air.TPX = 650, 25*ct.one_atm, {'O2':0.21, 'N2':0.79}

def run_finite_everything(tau_mix, tau_ent_main, tau_ent_sec, phi_global = 0.635, phi_main = 0.3719, airSplit = 1, tau_main = (20-0.158)*1e-3, tau_sec = 5*1e-3, dt = 0.001*1e-3, P = 25*ct.one_atm):
    [mfm, mam, mfs, mas] = solvePhi_airSplit(phi_global, phi_main, 100, airSplit)
    mass_main = mfm + mam
    mass_sec = mfs + mas
    mdot_main = mass_main/tau_ent_main
    mdot_sec = mass_sec/tau_ent_sec
    t = np.arange(0, tau_sec, dt)
    print(f"Length of t is {len(t)}")
    # Calculate main burner:
    [vit_reactor, main_burner_DF] = runMainBurner(phi_main, tau_main, P=P)    

    # Setup PaSBR:
    secondary_gas = mix([fuel, air], [mfs, mas], P = P)
    pasbr = bp.PaSBR([], N_MAX = 500, dt = dt)
    pfc_main = bp.ParticleFlowController(pasbr, vit_reactor.thermo, mass_main, dt, lambda t: mdot_main if t <= tau_ent_main else 0)
    pfc_sec = bp.ParticleFlowController(pasbr, secondary_gas, mass_sec, dt, lambda t: mdot_sec if t <= tau_ent_sec else 0)

    # Setup data management: 
    mean_gas = ct.Solution("gri30.xml")
    mean_gas.TPX = secondary_gas.TPX
    system_timeHistory = []
    num_particles_list = []
    # Simulation loop: 
    for t_now in t:
        [remaining_main_mass, main_state] = pfc_main.entrain(t_now)
        [remaining_sec_mass, sec_state] = pfc_sec.entrain(t_now)
        total_mass = remaining_main_mass + remaining_sec_mass + pasbr.mass
        system_state = (remaining_main_mass * main_state + remaining_sec_mass * sec_state + pasbr.mass * pasbr.state)/total_mass
        mean_gas.HPY = system_state[0], mean_gas.P, system_state[1:]
        system_timeHistory.append([t_now, total_mass, mean_gas.T, mean_gas.mean_molecular_weight, mean_gas.enthalpy_mass] + mean_gas.Y.tolist() + mean_gas.X.tolist())
        pasbr.react()
        pasbr.mix(tau_mix = tau_mix)
        num_particles_list.append(len(pasbr.particle_list))

    # Output data 
    entrainment_zone_df = pasbr.get_timeHistory()
    system_timeHistory = np.vstack(system_timeHistory)
    system_df = pd.DataFrame(columns=entrainment_zone_df.columns, data=system_timeHistory)
    system_df['num_particles'] = np.array(num_particles_list)

    return entrainment_zone_df, system_df, pasbr

def getNOx(sys_df, CO_constraint = 31.82):
    time = sys_df['age']
    NO = sys_df['X_NO']
    CO = sys_df['X_CO']
    H2O = sys_df['X_H2O']
    O2 = sys_df['X_O2']

    NO_corr = correctNOx(NO, H2O, O2)
    CO_corr = correctNOx(CO, H2O, O2)
    # constraint_ind = COLimitInd(CO_corr, 32)
    if len(CO_corr) > 0:
        constraint_ind = np.arange(0,len(CO_corr))[CO_corr >= CO_constraint + 1e-12][-1] + 1
    try:
        tau_sec = time[constraint_ind]
        NO_finalcorr = NO_corr[constraint_ind]
        CO_finalcorr = CO_corr[constraint_ind]
        T_corresponding = sys_df['T'].iloc[constraint_ind]
    except: # if constraint_ind is out of bounds, that means there's no point that meets the constraint, so return "inf" instead
        tau_sec = 10000
        NO_finalcorr = 10000
        CO_finalcorr = 10000
        T_corresponding = 0
    return NO_finalcorr, CO_finalcorr, tau_sec, T_corresponding    

def main():
    milliseconds = 1e-3;
    tau_ent_cf = np.array([0.1, 0.2, 1, 2, 3])*milliseconds
    tau_ent_sec = np.array([0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0])*milliseconds    
    tau_mix = np.array([0.01, 0.05, 0.1, 0.3])*milliseconds

    NO_list = []
    CO_list = []
    T_list = []
    tau_sec_required_list = []
    ent_ratio_list = []
    sys_df_list = []
    pasbr_df_list = []
    pasbr_list = []
    cf_list = []
    sec_list = []
    mix_list = []

    time_taken_list = []
    avg_particles_list = []
    for i in range(0, len(tau_ent_cf)):
        for j in range(0, len(tau_ent_sec)):
            for k in range(0, len(tau_mix)):
                t1 = time.time();
                pasbr_df, sys_df, pasbr = run_finite_everything(tau_mix[k], tau_ent_cf[i], tau_ent_sec[j], tau_sec=5.0*milliseconds, dt=0.001*milliseconds)
                t2 = time.time();
                sys_df_list.append(sys_df)
                pasbr_df_list.append(pasbr_df)
                pasbr_list.append(pasbr)
                NO, CO, tau_sec_required, T_corresponding = getNOx(sys_df)
                NO_list.append(NO)
                CO_list.append(CO)
                T_list.append(T_corresponding)
                tau_sec_required_list.append(tau_sec_required)
                ent_ratio_list.append(tau_ent_cf[i]/tau_ent_sec[j])

                cf_list.append(tau_ent_cf[i])
                sec_list.append(tau_ent_sec[j])
                mix_list.append(tau_mix[k])

                avg_particles_list.append(sys_df['num_particles'].mean())
                time_taken_list.append(t2-t1)
    
    data = np.vstack((cf_list, sec_list, mix_list, ent_ratio_list, NO_list, CO_list, T_list, tau_sec_required_list, avg_particles_list, time_taken_list))
    df = pd.DataFrame(data=np.transpose(data), columns = ['tau_ent_main', 'tau_ent_sec', 'tau_mix', 'ent_ratio', 'NO', 'CO', 'T', 'tau_sec_required', 'avg_num_particles', 'time_taken'])
    df.to_csv("finite_study.csv");
    dataFrame_to_pyarrow(df, "finite_study.pickle")
def test():

    print("Hello, world!")
    tau_mix = 0.01*1e-3
    tau_ent_main = 2.0*1e-3
    tau_ent_sec = 1.0*1e-3
    ent_ratio = tau_ent_main/tau_ent_sec
    # pasbr_df, sys_df, pasbr = run_finite_everything(tau_mix, tau_ent_main, tau_ent_sec)
    # NO, CO, tau_sec_required, T_corresponding = getNOx(sys_df)
    NO = 1
    CO = 1
    tau_sec_required = 1
    T_corresponding = 1
    df = pd.DataFrame(columns=['tau_mix', 'tau_ent_main', 'tau_ent_sec', 'ent_ratio', 'NO', 'CO', 'tau_sec_required', 'T'], data=np.hstack([tau_mix, tau_ent_main, tau_ent_sec, ent_ratio, NO, CO, tau_sec_required, T_corresponding]))
    df.to_csv("limited_everything_test.csv")

if __name__ == "__main__":
    t1 = time.time();    
    main()
    t2 = time.time()
    print(f"Time taken = {(t2-t1)/60:.2f} minutes")


