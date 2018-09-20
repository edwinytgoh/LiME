import sys
sys.path.insert(0, "../")
import os.path
from CanteraTools import *
import matplotlib.pyplot as plt

def COLimitInd(COHistory, constraint):
    return min(np.arange(len(COHistory))[COHistory > constraint][-1] + 1, len(COHistory)-1)

def runCase(tau_ent_cf, tau_ent_sec, toPickle = False):
    milliseconds = 0.001
    totalTime = 5.0*milliseconds
    t = np.arange(0, totalTime, 0.001*milliseconds)
    filename = 'infmix_cf_{0:0.2f}ms_sec_{1:0.2f}ms.pickle'.format(tau_ent_cf/milliseconds, tau_ent_sec/milliseconds)
    if os.path.isfile(filename):
        timetraceDF = pyarrow_to_dataFrame(filename)
    else:
        P = 25*ct.one_atm
        CO_constraint = 31.82 # corrected ppm
        [mfm, mam, mfs, mas] = solvePhi_airSplit(0.635, 0.3719, 100, 1)
        main_mass = mfm + mam
        jet_mass = mfs + mas
        [vit_reactor, main_burner_DF] = runMainBurner(0.3719, 19.842*milliseconds)

        secondary_gas = ct.Solution('gri30.xml')
        secondary_gas.TPX = 300, P, {'CH4':1}
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
            dataArray[i, :] = np.hstack([tnow, mean_gas.T, NOCOppmvd[0], NOCOppmvd[1], get_equivalence_ratio(mean_gas)])# mean_gas.get_equivalence_ratio()]) 

        timetraceDF = pd.DataFrame(data=dataArray, columns=columnNames)
        timetraceDF.set_index = 'time'
        if toPickle:
            # Pickling data if desired
            dataFrame_to_pyarrow(timetraceDF, filename)

    NO_corr = timetraceDF['NOppmvd']
    CO_corr = timetraceDF['COppmvd']
    constraint_ind = COLimitInd(CO_corr, CO_constraint)

    tau_sec = t[constraint_ind]
    NO_finalcorr = NO_corr[constraint_ind]
    CO_finalcorr = CO_corr[constraint_ind]
    return tau_sec, NO_finalcorr, CO_finalcorr

# Custom Equivalence Ratio test
def get_equivalence_ratio(gas, oxidizers = [], ignore = []):
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
    tau_ent_cf = np.array([0.1, 0.2, 0.5, 1, 2, 3])*milliseconds
    tau_ent_sec = np.array([0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0])*milliseconds
    
    NOs = np.zeros((len(tau_ent_cf), len(tau_ent_sec)))
    COs = np.zeros((len(tau_ent_cf), len(tau_ent_sec)))
    tau_sec = np.zeros((len(tau_ent_cf), len(tau_ent_sec)))
    for i in range(len(tau_ent_cf)):
        for j in range(len(tau_ent_sec)):
            (tau_sec[i, j], NOs[i, j], COs[i, j]) = runCase(tau_ent_cf[i], tau_ent_sec[j], toPickle = True)
            print('Reaction Time until constraint: ' + str(tau_sec[i, j]*1e3) + ' ms')
            print('Final Overall 15% O2 Corrected NO concentration: ' + str(NOs[i,j]) + ' ppm')    
    # Write to file
    np.savetxt("InfMix_FinEnt_NO.csv", NOs, delimiter=",")
    np.savetxt("InfMix_FinEnt_CO.csv", COs, delimiter=",")
    np.savetxt("InfMix_FinEnt_tau_sec.csv", tau_sec, delimiter=",")

    # Make some plots
    tau_ent_cf /= milliseconds
    tau_ent_sec /= milliseconds
    # Constant tau_ent_sec plots:
    fig1 = plt.figure()
    ax1 = plt.axes()
    ax1.set_title('Variation of Exit NO at CO Constraint over Constant tau_ent_sec')
    plt.xlabel('tau_ent_cf (ms)')
    plt.ylabel('Corrected NO Concentration to 15% O2 (ppm)')
    
    for i in range(len(tau_ent_sec)):
        ax1.plot(tau_ent_cf, NOs[:, i], '*-', label='tau_ent_sec = ' + str(tau_ent_sec[i]) + ' ms')
        # ax1.plot(tau_ent_cf, COs[:, i], label='tau_ent_sec = ' + str(tau_ent_sec[i]) + ' ms')
    ax1.grid(True)
    handles, labels = ax1.get_legend_handles_labels()
    lgd = ax1.legend(handles, labels, loc='upper right', bbox_to_anchor=(1.5,0.3))
    fig1.savefig('Const_tau_ent_sec.png', bbox_extra_artists=(lgd,), bbox_inches='tight')

    # Constant tau_ent_cf plots:
    fig2 = plt.figure()
    ax2 = plt.axes()
    ax2.set_title('Variation of Exit NO at CO Constraint over Constant tau_ent_cf')
    plt.xlabel('tau_ent_sec (ms)')
    plt.ylabel('Corrected NO Concentration to 15% O2 (ppm)')
    
    for i in range(len(tau_ent_cf)):
        ax2.plot(tau_ent_sec, NOs[i, :], '*-', label='tau_ent_cf = ' + str(tau_ent_cf[i]) + ' ms')
        # ax2.plot(tau_ent_sec, COs[i, :], label='tau_ent_cf = ' + str(tau_ent_cf[i]) + ' ms')
    ax2.grid(True)
    handles, labels = ax2.get_legend_handles_labels()
    lgd = ax2.legend(handles, labels, loc='upper right', bbox_to_anchor=(1.5,0.3))
    fig2.savefig('Const_tau_ent_cf.png', bbox_extra_artists=(lgd,), bbox_inches='tight')

    fig1 = plt.figure()
    ax1 = plt.axes()
    ax1.set_title('Variation of Exit CO over Constant tau_ent_sec')
    plt.xlabel('tau_ent_cf (ms)')
    plt.ylabel('Corrected CO Concentration to 15% O2 (ppm)')
    
    for i in range(len(tau_ent_sec)):
        ax1.plot(tau_ent_cf, COs[:, i], '*-', label='tau_ent_sec = ' + str(tau_ent_sec[i]) + ' ms')
    ax1.grid(True)
    handles, labels = ax1.get_legend_handles_labels()
    lgd = ax1.legend(handles, labels, loc='upper right', bbox_to_anchor=(1.5,0.3))
    fig1.savefig('Const_tau_ent_sec_CO.png', bbox_extra_artists=(lgd,), bbox_inches='tight')

    # Constant tau_ent_cf plots:
    fig2 = plt.figure()
    ax2 = plt.axes()
    ax2.set_title('Variation of Exit CO over Constant tau_ent_cf')
    plt.xlabel('tau_ent_sec (ms)')
    plt.ylabel('Corrected CO Concentration to 15% O2 (ppm)')
    
    for i in range(len(tau_ent_cf)):
        ax2.plot(tau_ent_sec, COs[i, :], '*-', label='tau_ent_cf = ' + str(tau_ent_cf[i]) + ' ms')
    ax2.grid(True)
    handles, labels = ax2.get_legend_handles_labels()
    lgd = ax2.legend(handles, labels, loc='upper right', bbox_to_anchor=(1.5,0.3))
    fig2.savefig('Const_tau_ent_cf_CO.png', bbox_extra_artists=(lgd,), bbox_inches='tight')

    # tau_sec plots
    fig3 = plt.figure()
    ax3 = plt.axes()
    ax3.set_title('Variation of time to meet CO constrained (tau_sec) over Constant tau_ent_sec')
    plt.xlabel('tau_ent_cf (ms)')
    plt.ylabel('tau_sec (ms)')
    for i in range(len(tau_ent_sec)):
        ax3.plot(tau_ent_cf, tau_sec[:, i]/milliseconds, '*-', label='tau_ent_sec = ' + str(tau_ent_sec[i]) + ' ms')
    ax3.grid(True)
    handles, labels = ax3.get_legend_handles_labels()
    lgd = ax3.legend(handles, labels, loc='upper right', bbox_to_anchor=(1.5,0.3))
    fig3.savefig('tau_sec_v_ent_cf.png', bbox_extra_artists=(lgd,), bbox_inches='tight')

    fig4 = plt.figure()
    ax4 = plt.axes()
    ax4.set_title('Variation of time to meet CO constrained (tau_sec) over Constant tau_ent_cf')
    plt.xlabel('tau_ent_sec (ms)')
    plt.ylabel('tau_sec (ms)')
    for i in range(len(tau_ent_cf)):
        ax4.plot(tau_ent_sec, tau_sec[i, :]/milliseconds, '*-', label='tau_ent_cf = ' + str(tau_ent_cf[i]) + ' ms')
    ax4.grid(True)
    handles, labels = ax4.get_legend_handles_labels()
    lgd = ax4.legend(handles, labels, loc='upper right', bbox_to_anchor=(1.5,0.3))
    fig4.savefig('tau_sec_v_ent_sec.png', bbox_extra_artists=(lgd,), bbox_inches='tight')

    # Entrainment ratio plots
    ent_ratio = np.zeros((len(tau_ent_cf), len(tau_ent_sec)))
    for i in range(len(tau_ent_cf)):
        ent_ratio[i, :] = tau_ent_sec/tau_ent_cf[i]
    fig5, (ax5, ax6) = plt.subplots(nrows=2, ncols=1)
    ax5.set_title('NO, CO v. Entrainment Ratio')
    ax6.set_title('tau_sec v. Entrainment Ratio')
    plt.xlabel('Entrainment Ratio')
    
    order = np.argsort(ent_ratio.flatten())
    ax5.plot(ent_ratio.flatten()[order], NOs.flatten()[order], '*-', label='15% O2 Corrected NO')
    ax5.plot(ent_ratio.flatten()[order], COs.flatten()[order], '*-', label='15% O2 Corrected CO')
    ax6.plot(ent_ratio.flatten()[order], (tau_sec.flatten()/milliseconds)[order], '*-', label='tau_sec')
    ax5.grid(True)
    ax6.grid(True)
    handles, labels = ax5.get_legend_handles_labels()
    lgd = ax5.legend(handles, labels, loc='upper right', bbox_to_anchor=(1.5,0.7))
    fig5.savefig('ent_ratio.png', bbox_extra_artists=(lgd,), bbox_inches='tight')

    # plt.show()

if __name__ == "__main__":
    main()