import sys
sys.path.insert(0, "../")
import os
from CanteraTools import *
import matplotlib.pyplot as plt

P = 1.0*ct.one_atm

sec_fuel = ct.Solution('gri30.xml')
sec_fuel.TPX = 440, P, {'CH4':1}
sec_air = ct.Solution('gri30.xml')
sec_air.TPX = 440, P, {'O2':0.21, 'N2':0.79}

def COLimitInd(COHistory, constraint):
    return min(np.arange(len(COHistory))[COHistory > constraint][-1] + 1, len(COHistory)-1)

def runCase(phi_HE, J_ratio, phi_jet, tau_ent_cf = 0.1*1e-3, tau_ent_sec = 0.1*1e-3, toPickle = False):
    milliseconds = 0.001
    totalTime = 52.0*milliseconds
    t = np.arange(0, totalTime, 0.001*milliseconds)
    filename = 'pickles/exp_cf_{0:0.2f}ms_sec_{1:0.2f}ms.pickle'.format(tau_ent_cf/milliseconds, tau_ent_sec/milliseconds)
    if not os.path.exists('pickles'):
        os.makedirs('pickles')
    
    if os.path.isfile(filename):
        timetraceDF = pyarrow_to_dataFrame(filename)
    else:
        # CO_constraint = 31.82 # corrected ppm
        [vit_reactor, main_burner_DF] = runMainBurner(phi_HE, 1.17/17.3, T_fuel=591, T_ox=591, P=P)    # Mixed temperature around 591 K
        fs = 0.058387057492574147
        secondary_gas = mix([sec_fuel, sec_air], [fs*phi_jet, 1], P = P)

        mam = 1
        mfm = mam*fs*phi_HE
        u_j = (J_ratio*17.3**2*vit_reactor.thermo.density / secondary_gas.density)**0.5
        K = J_ratio*(17.3*np.pi*(0.012/2)**2)/(u_j*0.0635*0.114) * (1+phi_HE*fs)/(1+phi_jet*fs)
        phi_global = (phi_HE+K*phi_jet)/(1+K)
        mas = mam * (phi_global-phi_HE) / (phi_jet-phi_global)
        mfs = mas * phi_jet*fs

        # Normalize
        mtot = mam + mfm + mas + mfs
        mam *= 100/mtot
        mfm *= 100/mtot
        mas *= 100/mtot
        mfs *= 100/mtot

        # [mfm, mam, mfs, mas] = solvePhi_airSplit(phi_HE + delta_phi, phi_HE, 100, airsplit)
        main_mass = mfm + mam
        jet_mass = mfs + mas
        
        # secondary_gas = ct.Solution('gri30.xml')
        # Using mass fractions since we have those
        sec_reactor = ct.ConstPressureReactor(secondary_gas)

        interface_gas = mix([vit_reactor.thermo, secondary_gas], [main_mass, jet_mass])
        initial_secReactor_mass = 1e-6 # kg
        interface_reactor = ct.ConstPressureReactor(interface_gas, volume = initial_secReactor_mass/interface_gas.density) # interface reactor is really the secondary stage in an AFS

        mfc_vit = ct.MassFlowController(vit_reactor, interface_reactor)
        mdot_vit = main_mass/tau_ent_cf
        mfc_vit.set_mass_flow_rate(lambda t: mdot_vit if t <= tau_ent_cf else 0)

        mfc_sec = ct.MassFlowController(sec_reactor, interface_reactor)
        mdot_sec = jet_mass/tau_ent_sec
        mfc_sec.set_mass_flow_rate(lambda t: mdot_sec if t <= tau_ent_sec else 0)

        reactorNet = ct.ReactorNet([vit_reactor, sec_reactor, interface_reactor])
        mean_gas = ct.Solution('gri30.xml')
        columnNames = ['time', 'T', 'NOppmvd', 'COppmvd', 'phi_global']
        dataArray = np.array([None] * len(t) * len(columnNames)).reshape(len(t), len(columnNames))
        for i in range(len(t)):
            tnow = t[i]
            reactorNet.advance(tnow)
            # Get the number of moles in each reactor since we're working with mole fractions
            mass_vit = main_mass - mdot_vit*min(tnow, tau_ent_cf)
            mass_sec = jet_mass - mdot_sec*min(tnow, tau_ent_sec)
            mass_int = main_mass - mass_vit + jet_mass - mass_sec
            massFrac = ((interface_reactor.thermo.Y * mass_int + vit_reactor.thermo.Y * mass_vit + sec_reactor.thermo.Y * mass_sec) /
                        (mass_int + mass_vit + mass_sec))
            enthalpy = ((interface_reactor.thermo.enthalpy_mass * mass_int + vit_reactor.thermo.enthalpy_mass * mass_vit + sec_reactor.thermo.enthalpy_mass * mass_sec) / 
                        (mass_int + mass_vit + mass_sec))
            mean_gas.HPY = enthalpy, 25*ct.one_atm, massFrac

            species_NOCO = mean_gas['NO', 'CO'].X
            species_O2 = mean_gas['O2'].X
            species_H2O = mean_gas['H2O'].X
            NOCOppmvd = correctNOx(species_NOCO, species_H2O, species_O2)
            dataArray[i, :] = np.hstack([tnow, mean_gas.T, NOCOppmvd[0], NOCOppmvd[1], mean_gas.get_equivalence_ratio()])# mean_gas.get_equivalence_ratio()]) 

        timetraceDF = pd.DataFrame(data=dataArray, columns=columnNames)
        timetraceDF.set_index = 'time'
        if toPickle:
            # Pickling data if desired
            dataFrame_to_pyarrow(timetraceDF, filename)

    NO_corr = timetraceDF['NOppmvd']
    CO_corr = timetraceDF['COppmvd']
    # constraint_ind = COLimitInd(CO_corr, CO_constraint)
    tau_sec = t[-1]
    NO_finalcorr = NO_corr.iat[-1]
    CO_finalcorr = CO_corr.iat[-1]
    return tau_sec, NO_finalcorr, CO_finalcorr, phi_global

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

# Defining run cases here
def main():
    milliseconds = 1e-3
    tau_ent_cf = 40*milliseconds
    tau_ent_sec = 25*milliseconds
    
    J_ratio = np.array([1.4, 2.1, 2.8, 3.6])
    phi_jet = np.arange(2, 9, 1)

    NOs = np.zeros((len(phi_jet), len(J_ratio)))
    COs = np.zeros((len(phi_jet), len(J_ratio)))
    for i in range(len(phi_jet)):
        for j in range(len(J_ratio)):
            (dummy, NOs[i, j], COs[i, j], phi_global) = runCase(0.625, J_ratio[j], phi_jet[i], tau_ent_cf, tau_ent_sec, toPickle = False)
            print('Final Overall 15% O2 Corrected NO concentration: ' + str(NOs[i,j]) + ' ppm')
            print('Phi_Global of ' + str(phi_global))    

    # Make some plots
    fig1 = plt.figure()
    ax1 = plt.axes()
    ax1.set_title('NO at Exit over Jet Parameters')
    plt.xlabel('Jet Equivalence Ratio')
    plt.ylabel('Corrected NO Concentration to 15% O2 (ppm)')
    
    for i in range(len(J_ratio)):
        ax1.plot(phi_jet, NOs[:, i], '*', label='J = ' + str(J_ratio[i]))
    ax1.grid(True)
    handles, labels = ax1.get_legend_handles_labels()
    lgd = ax1.legend(handles, labels, loc='upper right', bbox_to_anchor=(1.5,0.3))
    fig1.savefig('Phi_HE_0.625.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    main()