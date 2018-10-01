import sys
sys.path.insert(0, "../")
import os
from CanteraTools import *
import matplotlib.pyplot as plt

P = 25.0*ct.one_atm

sec_fuel = ct.Solution('gri30.xml')
sec_fuel.TPX = 300, P, {'CH4':1}
sec_air = ct.Solution('gri30.xml')
sec_air.TPX = 650, P, {'O2':0.21, 'N2':0.79}

def COLimitInd(COHistory, constraint):
    return min(np.arange(len(COHistory))[COHistory > constraint][-1] + 1, len(COHistory)-1)

def runCase(phi_jet, phi_global = 0.635, tau_ent_cf = 0.3*1e-3, tau_ent_sec = 0.2*1e-3, toPickle = False):
    milliseconds = 0.001
    totalTime = 2.0*milliseconds
    t = np.arange(0, totalTime, 0.001*milliseconds)
    CO_constraint = 31.82 # corrected ppm
    filename = 'pickles/fin_jet_{0:0.2f}_tau_ent_{1:0.2f}_{2:0.2f}.pickle'.format(phi_jet, tau_ent_cf/milliseconds, tau_ent_sec/milliseconds)
    if not os.path.exists('pickles'):
        os.makedirs('pickles')
    
    if os.path.isfile(filename):
        timetraceDF = pyarrow_to_dataFrame(filename)
    else:
        phi_HE = 0.3719
        [vit_reactor, main_burner_DF] = runMainBurner(phi_HE, 19.842*milliseconds)    # Mixed temperature around 591 K
        fs = 0.058387057492574147
        secondary_gas = mix([sec_fuel, sec_air], [fs*phi_jet, 1], P = P)

        mam = 1
        mfm = mam*fs*phi_HE
        mfs1 = mam*fs*(phi_global - phi_HE)
        mas = mfs1/(fs*(phi_jet-phi_global))
        mfs = mfs1 + mas*fs*phi_global
        # Normalize
        mtot = mam + mfm + mas + mfs
        mam *= 100/mtot
        mfm *= 100/mtot
        mas *= 100/mtot
        mfs *= 100/mtot

        main_mass = mfm + mam
        jet_mass = mfs + mas

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
            mean_gas.HPY = enthalpy, P, massFrac

            species_NOCO = mean_gas['NO', 'CO'].X
            species_O2 = mean_gas['O2'].X
            species_H2O = mean_gas['H2O'].X
            NOCOppmvd = correctNOx(species_NOCO, species_H2O, species_O2)
            dataArray[i, :] = np.hstack([tnow, mean_gas.T, NOCOppmvd[0], NOCOppmvd[1], get_global_equivalence_ratio(mean_gas)])# mean_gas.get_equivalence_ratio()]) 

        timetraceDF = pd.DataFrame(data=dataArray, columns=columnNames)
        timetraceDF.set_index = 'time'
        if toPickle:
            # Pickling data if desired
            dataFrame_to_pyarrow(timetraceDF, filename)

    NO_corr = timetraceDF['NOppmvd']
    CO_corr = timetraceDF['COppmvd']
    constraint_ind = COLimitInd(CO_corr, CO_constraint)
    tau_sec = t[constraint_ind]
    NO_finalcorr = NO_corr.iat[constraint_ind]
    CO_finalcorr = CO_corr.iat[constraint_ind]
    phi_global = timetraceDF.phi_global.iat[constraint_ind]
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
    phi_jet = np.arange(2, 12, 2.0)
    tau_ent_cf = np.array([2.0])*milliseconds
    tau_ent_sec = np.array([0.5, 1.0, 2.0])*milliseconds

    NOs = np.zeros((len(phi_jet), len(tau_ent_cf), len(tau_ent_sec)))
    COs = np.zeros((len(phi_jet), len(tau_ent_cf), len(tau_ent_sec)))
    phi_global = np.zeros((len(phi_jet), len(tau_ent_cf), len(tau_ent_sec)))
    for i in range(len(phi_jet)):
        for j in range(len(tau_ent_cf)):
            for k in range(len(tau_ent_sec)):
                (dummy, NOs[i, j, k], COs[i, j, k], phi_global[i, j, k]) = runCase(phi_jet[i], tau_ent_cf = tau_ent_cf[j], tau_ent_sec = tau_ent_sec[k], toPickle = True)
                print('Final Overall 15% O2 Corrected NO concentration: ' + str(NOs[i, j, k]) + ' ppm')
                print('Phi_Global of ' + str(phi_global[i, j, k]))    

    # Make some plots
    fig1 = plt.figure()
    ax1 = plt.axes()
    ax1.set_title('NO at CO Constraint')
    plt.xlabel('Jet Equivalence Ratio')
    plt.ylabel('Corrected NO Concentration to 15% O2 (ppm)')
    
    for j in range(len(tau_ent_cf)):
        for k in range(len(tau_ent_sec)):
            ax1.plot(phi_jet, NOs[:, j, k], '*-', label='cf={0:0.2f} ms; sec={1:0.2f} ms'.format(tau_ent_cf[j]/milliseconds, tau_ent_sec[k]/milliseconds))
    ax1.grid(True)
    handles, labels = ax1.get_legend_handles_labels()
    lgd = ax1.legend(handles, labels, loc='upper right', bbox_to_anchor=(1.5,0.3))
    fig1.savefig('Fin_Phi_Sec_NO.png', bbox_extra_artists=(lgd,), bbox_inches='tight')

    fig2 = plt.figure()
    ax2 = plt.axes()
    ax2.set_title('Phi Global')
    plt.xlabel('Jet Equivalence Ratio')
    plt.ylabel('Phi Global')
    
    for j in range(len(tau_ent_cf)):
        for k in range(len(tau_ent_sec)):
            ax2.plot(phi_jet, phi_global[:, j, k], '*', label='cf={0:0.2f} ms; sec={1:0.2f} ms'.format(tau_ent_cf[j]/milliseconds, tau_ent_sec[k]/milliseconds))
    ax2.grid(True)
    ax2.set_ylim((0, 1))
    handles, labels = ax2.get_legend_handles_labels()
    lgd = ax2.legend(handles, labels, loc='upper right', bbox_to_anchor=(1.5,0.3))
    # fig2.savefig('Phi_HE_0.625_phi_globals.png', bbox_extra_artists=(lgd,), bbox_inches='tight')

    plt.show()

if __name__ == "__main__":
    main()