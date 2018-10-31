import sys
sys.path.insert(0, "../")
import os
from CanteraTools import *

P = 25.0*ct.one_atm
sec_fuel = ct.Solution('gri30.xml')
sec_fuel.TPX = 300, P, {'CH4':1}
sec_air = ct.Solution('gri30.xml')
sec_air.TPX = 650, P, {'O2':0.21, 'N2':0.79}

column_names = ['age', 'mass', 'T', 'MW', 'h', 'phi'] + ["Y_" + sn for sn in sec_fuel.species_names] + ["X_" + sn for sn in sec_fuel.species_names]    

def get_timeHistory(timeHistory_list, dataFrame=True):

    timeHistory_array = np.vstack(timeHistory_list)
    if dataFrame == True:
        df = pd.DataFrame(columns = column_names, data = timeHistory_array)
        df.set_index(['age'])
        return df
    
    return timeHistory_array  

def masssplit(phi_main, phi_global, phi_jet):
    fs = 0.058387057492574147
    mam = 1
    mfm = mam*fs*phi_main
    mfs1 = mam*fs*(phi_global - phi_main)
    mas = mfs1/(fs*(phi_jet-phi_global))
    mfs = mfs1 + mas*fs*phi_global
    # Normalize
    mtot = mam + mfm + mas + mfs
    mam *= 100/mtot
    mfm *= 100/mtot
    mas *= 100/mtot
    mfs *= 100/mtot
    return mam, mfm, mas, mfs

def runCase(tau_main = 20-0.158, tau_sec = 5.0, phi_global = 0.635, phi_main = 0.3719, phi_jet = np.inf, enttype = 'time', ent = (2.0, 0.5)):
    milliseconds = 0.001
    totalTime = tau_sec*milliseconds
    t = np.arange(0, totalTime, 0.001*milliseconds)
    
    if enttype == 'time':
        tau_ent_cf = ent[0]*milliseconds
        tau_ent_sec = ent[1]*milliseconds
    elif enttype == 'mass':
        mdot_vit = ent[0]/milliseconds
        mdot_sec = ent[1]/milliseconds
    else:
        print('Entrainment type must be ''mass'' or ''time''')
        return -1
    
    [vit_reactor, main_burner_DF] = runMainBurner(phi_main, tau_main*milliseconds)    # Mixed temperature around 591 K
    fs = 0.058387057492574147
    secondary_gas = mix([sec_fuel, sec_air], [fs*phi_jet, 1], P = P)
    mam, mfm, mas, mfs = masssplit(phi_main, phi_global, phi_jet)
    main_mass = mfm + mam
    jet_mass = mfs + mas

    # Using mass fractions since we have those
    sec_reactor = ct.ConstPressureReactor(secondary_gas)

    interface_gas = mix([vit_reactor.thermo, secondary_gas], [main_mass, jet_mass])
    initial_mass = 1e-6 # kg
    interface_reactor = ct.ConstPressureReactor(interface_gas, volume = initial_mass/interface_gas.density) # interface reactor is really the secondary stage in an AFS

    mfc_vit = ct.MassFlowController(vit_reactor, interface_reactor)
    if enttype == 'time':
        mdot_vit = main_mass/tau_ent_cf
    elif enttype == 'mass':
        tau_ent_cf = main_mass/mdot_vit
    mfc_vit.set_mass_flow_rate(lambda t: mdot_vit if t <= tau_ent_cf else 0)

    mfc_sec = ct.MassFlowController(sec_reactor, interface_reactor)
    if enttype == 'time':
        mdot_sec = jet_mass/tau_ent_sec
    elif enttype == 'mass':
        tau_ent_sec = jet_mass/mdot_sec
    mfc_sec.set_mass_flow_rate(lambda t: mdot_sec if t <= tau_ent_sec else 0)

    reactorNet = ct.ReactorNet([vit_reactor, sec_reactor, interface_reactor])
    mean_gas = ct.Solution('gri30.xml')
    # Run CRN to completion
    system_timeHistory = []
    reactor_timeHistory = []
    for i in range(len(t)):
        tnow = t[i]
        reactorNet.advance(tnow)
        # Get the number of moles in each reactor since we're working with mole fractions
        mass_vit = max(main_mass - mdot_vit*tnow, 0)
        mass_sec = max(jet_mass - mdot_sec*tnow, 0)
        mass_int = main_mass - mass_vit + jet_mass - mass_sec
        massFrac = ((interface_reactor.thermo.Y * mass_int + vit_reactor.thermo.Y * mass_vit + sec_reactor.thermo.Y * mass_sec) /
                    (mass_int + mass_vit + mass_sec))
        enthalpy = ((interface_reactor.thermo.enthalpy_mass * mass_int + vit_reactor.thermo.enthalpy_mass * mass_vit + sec_reactor.thermo.enthalpy_mass * mass_sec) / 
                    (mass_int + mass_vit + mass_sec))
        
        mean_gas.HPY = enthalpy, P, massFrac
        int_gas = interface_reactor.thermo

        system_timeHistory.append([tnow, mass_int+mass_vit+mass_sec, mean_gas.T, mean_gas.mean_molecular_weight, mean_gas.enthalpy_mass, mean_gas.get_equivalence_ratio()] + mean_gas.Y.tolist() + mean_gas.X.tolist())
        reactor_timeHistory.append([tnow, mass_int, int_gas.T, int_gas.mean_molecular_weight, int_gas.enthalpy_mass, int_gas.get_equivalence_ratio()] + int_gas.Y.tolist() + int_gas.X.tolist())

    # Output data
    entrainment_zone_df = get_timeHistory(reactor_timeHistory)
    system_timeHistory = np.vstack(system_timeHistory)
    system_df = pd.DataFrame(columns=column_names, data=system_timeHistory)

    return entrainment_zone_df, system_df

# Custom Equivalence Ratio test
def get_global_equivalence_ratio(gas, oxidizers = [], ignore = []):
    if not oxidizers:  # Default behavior, find all possible oxidizers
        oxidizers = [s.name for s in gas.species()] # if
                    # all(y not in s.composition for y in ['C', 'H', 'S'])]
        alpha = 0
        mol_O = 0
        for k, s in enumerate(gas.species()):
            if s.name in ignore:
                continue
            else: # elif s.name in oxidizers:
                mol_O += s.composition.get('O', 0) * gas.X[k]
            # else:
                nC = s.composition.get('C', 0)
                nH = s.composition.get('H', 0)
                # nO = s.composition.get('O', 0)
                nS = s.composition.get('S', 0)
                alpha += (2*nC + nH/2 + 2*nS) * gas.X[k]

        if mol_O == 0:
            return float('inf')
        else:
            return alpha / mol_O

def getNOx(sys_df, CO_constraint = 31.82):
    time = sys_df['age']
    NO = sys_df['X_NO']
    CO = sys_df['X_CO']
    H2O = sys_df['X_H2O']
    O2 = sys_df['X_O2']
    NO_corr = correctNOx(NO, H2O, O2)
    CO_corr = correctNOx(CO, H2O, O2)
    if len(CO_corr) > 0:
        mask = np.arange(0,len(CO_corr))[CO_corr >= CO_constraint + 1e-12]
        if len(mask) > 0:
            constraint_ind = mask[-1] + 1
        else:
            constraint_ind = len(CO_corr)
            print("We couldn't find a CO in the entire system time history :(")
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

def outputHandler(enttype, ent_main, ent_sec, out_dir, tau_sec=5.0, phi_jet=np.inf):
    if not out_dir[-1] == "/":
        out_dir += "/"    
    if np.isfinite(phi_jet):
        filename = f"premix-entMain_{ent_main:.3f}-entSec_{ent_sec:.3f}-phiJet_{phi_jet:.1f}"
    else:
        filename = f"premix-entMain_{ent_main:.3f}-entSec_{ent_sec:.3f}"
    if os.path.isfile(out_dir + filename + ".csv"):
        print("Found file " + out_dir + filename + ".csv" + ". Exiting...")
        return

    NO_list = []
    CO_list = []
    T_list = []
    tau_sec_required_list = []
    ent_ratio_list = []
    main_list = []
    sec_list = []
    ent_type_list = []

    reactor_df, sys_df = runCase(tau_sec = tau_sec, phi_jet = phi_jet, enttype = enttype, ent = (ent_main, ent_sec))
    NO, CO, tau_sec_required, T_corresponding = getNOx(sys_df)

    NO_list.append(NO)
    CO_list.append(CO)
    T_list.append(T_corresponding)
    tau_sec_required_list.append(tau_sec_required/milliseconds)
    ent_ratio_list.append(ent_main/ent_sec)
    main_list.append(ent_main)
    sec_list.append(ent_sec)
    ent_type_list.append(enttype)
    
    data = np.vstack((ent_type_list, main_list, sec_list, ent_ratio_list, NO_list, CO_list, T_list, tau_sec_required_list))
    df = pd.DataFrame(data=np.transpose(data), columns = ['ent_type', 'ent_main', 'ent_sec', 'ent_ratio', 'NO', 'CO', 'T', 'tau_sec_required'])
    df.to_csv(out_dir + filename + ".csv")
    dataFrame_to_pyarrow(df, out_dir + filename + ".pickle")
    dataFrame_to_pyarrow(sys_df, out_dir + "sys_df_" + filename + ".pickle")
    dataFrame_to_pyarrow(reactor_df, out_dir + "reactor_df_" + filename + ".pickle")
    return df, sys_df, reactor_df

# Defining run cases here
def main():
    milliseconds = 1e-3
    phi_jet = np.arange(2, 12, 2.0)
    entraintype = 'time'    # use 'time' or 'mass'
    tau_ent_cf = np.array([2.0])*milliseconds
    tau_ent_sec = np.array([0.5])*milliseconds
    mass_ent_cf = np.array([20, 30, 40, 50, 60, 100])*1e3
    mass_ent_sec = np.array([10])*1e3
    if entraintype is 'time':
        ent_cf = tau_ent_cf
        ent_sec = tau_ent_sec
    elif entraintype is 'mass':
        ent_cf = mass_ent_cf
        ent_sec = mass_ent_sec
    else:
        print('Use ''time'' or ''mass''  entraintype only')
        return 1

    ilen = len(phi_jet)
    jlen = len(ent_cf)
    klen = len(ent_sec)

    NOs = np.zeros((ilen, jlen, klen))
    COs = np.zeros((ilen, jlen, klen))
    phi_global = np.zeros((ilen, jlen, klen))
    for i in range(ilen):
        for j in range(jlen):
            for k in range(klen):
                (dummy, NOs[i, j, k], COs[i, j, k], phi_global[i, j, k]) = runCase(phi_jet[i], enttype = entraintype, ent = (ent_cf[j], ent_sec[k]), toPickle = True)
                print('Final Overall 15% O2 Corrected NO concentration: ' + str(NOs[i, j, k]) + ' ppm')

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("enttype", type=str)
    parser.add_argument("ent_main_low", type=float)
    parser.add_argument("ent_main_upp", type=float)
    parser.add_argument("ent_main_step", type=float)
    parser.add_argument("ent_sec_low", type=float)
    parser.add_argument("ent_sec_upp", type=float)
    parser.add_argument("ent_sec_step", type=float)
    parser.add_argument("out_dir", type=str)
    parser.add_argument("tau_sec", type=float)
    parser.add_argument("phi_jet", nargs='?', type=float, default=np.inf)
    args = parser.parse_args()
    enttype = args.enttype
    ent_main_low=args.ent_main_low
    ent_main_upp=args.ent_main_upp
    ent_main_step=args.ent_main_step
    ent_sec_low=args.ent_sec_low
    ent_sec_upp=args.ent_sec_upp
    ent_sec_step=args.ent_sec_step
    out_dir=args.out_dir
    tau_sec=args.tau_sec
    phi_jet=args.phi_jet

    milliseconds = 1e-3
    if enttype == 'time' or enttype == 'mass':
        ent_cf = np.arange(ent_main_low, ent_main_upp+ent_main_step, ent_main_step)
        ent_sec = np.arange(ent_sec_low, ent_sec_upp+ent_sec_step, ent_sec_step)
    else:
        print('Use ''time'' or ''mass''  enttype only')

    ilen = len(ent_cf)
    jlen = len(ent_sec)

    # For full run, ignore tau_ent_cf < tau_ent_sec
    if enttype is main:
        mam, mfm, mas, mfs = masssplit(0.3719, 0.635, phi_jet)
        main_mass = mam+mfm
        sec_mass = mas+mfs

    print(os.getcwd())
    for i in range(ilen):
        for j in range(jlen):
            if (enttype == 'time' and ent_cf[i] >= ent_sec[j]) or (enttype == 'mass' and (main_mass/ent_cf[i]) >= (sec_mass/ent_sec[j])):
                outputHandler(enttype, ent_cf[i], ent_sec[j], out_dir, tau_sec=tau_sec, phi_jet=phi_jet)
 