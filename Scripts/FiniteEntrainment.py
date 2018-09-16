import sys
sys.path.insert(0, "../")
from CanteraTools import *
import matplotlib.pyplot as plt

def COLimitInd(COHistory, constraint):

    return np.arange(len(COHistory))[COHistory > constraint][-1] + 1
    # index = len(COHistory) - 1
    # while COHistory[index] < constraint:
    #     index -= 1
    # index += 1  # Go back one so within constraint
    # index = min(index, len(COHistory))    # Avoid index out of bounds
    # return index

fuel = ct.Solution('gri30.xml')
fuel.TPX = 300, P, {'CH4':1} # TODO: Come up with a more general way of doing secondary gas so that we can have both fuel and air 
air = ct.Solution('gri30.xml'); 
air.TPX = 650, P, {'O2':0.21, 'N2':0.79}

def runCase(tau_ent_cf, tau_ent_sec):
    milliseconds = 0.001
    P = 25*ct.one_atm
    CO_constraint = 31.82 # corrected ppm
    [mfm, mam, mfs, mas] = solvePhi_airSplit(0.635, 0.3719, 100, 1)
    main_mass = mfm + mam
    jet_mass = mfs + mas
    [vit_reactor, main_burner_DF] = runMainBurner(0.3719, 19.842*milliseconds)

    secondary_gas = ct.Solution('gri30.xml')
    secondary_gas.TPX = 300, P, {'CH4':1}
    sec_reactor = ct.ConstPressureReactor(secondary_gas)

    totalTime = 3.0*milliseconds
    interface_gas = mix([vit_reactor.thermo, secondary_gas], [main_mass, jet_mass])
    initial_secReactor_mass = 1e-6 # kg
    interface_reactor = ct.ConstPressureReactor(interface_gas, volume = initial_secReactor_mass/interface_gas.density) # interface reactor is really the secondary stage in an AFS


    mfc_vit = ct.MassFlowController(vit_reactor, interface_reactor)
    mdot_vit = (mfm+mam)/tau_ent_cf
    mfc_vit.set_mass_flow_rate(lambda t: mdot_vit if t <= tau_ent_cf else 0)

    mfc_sec = ct.MassFlowController(sec_reactor, interface_reactor)
    mdot_sec = (mfs+mas)/tau_ent_sec
    mfc_sec.set_mass_flow_rate(lambda t: mdot_sec if t <= tau_ent_sec else 0)

    reactorNet = ct.ReactorNet([vit_reactor, sec_reactor, interface_reactor])
    t = np.arange(0, totalTime, 0.001*milliseconds)
    mean_gas = ct.Solution('gri30.xml')
    # states = []
    species_NO = np.array([])
    species_CO = np.array([])
    species_O2 = np.array([])
    species_H2O = np.array([])
    columnNames = ['time', 'T', 'NOppmvd', 'COppmvd']
    dataArray = np.array([None] * len(t) * len(columnNames)).reshape(len(t), len(columnNames))
    # enthalpy = []

    for i in range(len(t)):
        tnow = t[i]
        reactorNet.advance(tnow)

        # Get the number of moles in each reactor since we're working with mole fractions
        mass_vit = (mam + mfm) - mdot_vit*min(tnow, tau_ent_cf)
        mass_sec = (mas + mfs) - mdot_sec*min(tnow, tau_ent_sec)
        mass_int = (mam + mfm) - mass_vit + (mas + mfs) - mass_sec
        massFrac = ((interface_reactor.thermo.Y * mass_int + vit_reactor.thermo.Y * mass_vit + sec_reactor.thermo.Y * mass_sec) /
                    (mass_int + mass_vit + mass_sec))
        enthalpy = ((interface_reactor.thermo.enthalpy_mass * mass_int + vit_reactor.thermo.enthalpy_mass * mass_vit + sec_reactor.thermo.enthalpy_mass * mass_sec) / 
                    (mass_int + mass_vit + mass_sec))
        mean_gas.HPY = enthalpy, 25*ct.one_atm, massFrac

        species_NOCO = mean_gas['NO', 'CO'].X
        species_O2 = mean_gas['O2'].X
        species_H2O = mean_gas['H2O'].X
        NOCOppmvd = correctNOx(species_NOCO, species_H2O, species_O2)
        
        dataArray[i, :] = np.hstack([tnow, mean_gas.T, NOCOppmvd[0], NOCOppmvd[1]]) 

    timetraceDF = pd.DataFrame(data=dataArray, columns=columnNames)
    timetraceDF.set_index = 'time'
    table = pa.Table.from_pandas(timetraceDF)
    filename = 'infmix_cf_{0:.2f}ms_sec_{1:.2f}ms.pickle'.format(tau_ent_cf/milliseconds, tau_ent_sec/milliseconds)
    pq.write_table(table, filename)

    NO_corr = dataArray[:, 2]
    CO_corr = dataArray[:, 3]
    constraint_ind = COLimitInd(CO_corr, 32)

    tau_sec = t[constraint_ind]
    NO_finalcorr = NO_corr[constraint_ind]
    return tau_sec, NO_finalcorr

# Defining run cases here
def main():
    milliseconds = 1e-3
    tau_ent_cf = np.array([0.1, 0.2, 1, 2, 3])*milliseconds
    tau_ent_sec = np.array([0.05, 0.1, 0.2, 0.5])*milliseconds
    
    NOs = np.zeros((len(tau_ent_cf), len(tau_ent_sec)))
    tau_sec = np.zeros((len(tau_ent_cf), len(tau_ent_sec)))
    for i in range(len(tau_ent_cf)):
        for j in range(len(tau_ent_sec)):
            (tau_sec[i, j], NO_finalcorr) = runCase(tau_ent_cf[i], tau_ent_sec[j])
            NOs[i, j] = NO_finalcorr
            print('Reaction Time until constraint: ' + str(tau_sec[i, j]*1e3) + ' ms')
            print('Final Overall 15% O2 Corrected NO concentration: ' + str(NO_finalcorr) + ' ppm')    
    # Write to file
    np.savetxt("InfMix_FinEnt.csv", NOs, delimiter=",")

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
        ax1.plot(tau_ent_cf, NOs[:, i], label='tau_ent_sec = ' + str(tau_ent_sec[i]*1e3) + ' ms')
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
        ax2.plot(tau_ent_sec, NOs[i, :], label='tau_ent_cf = ' + str(tau_ent_cf[i]*1e3) + ' ms')
    ax2.grid(True)
    handles, labels = ax2.get_legend_handles_labels()
    lgd = ax2.legend(handles, labels, loc='upper right', bbox_to_anchor=(1.5,0.3))
    fig2.savefig('Const_tau_ent_cf.png', bbox_extra_artists=(lgd,), bbox_inches='tight')

    # tau_sec plots
    fig3 = plt.figure()
    ax3 = plt.axes()
    ax3.set_title('Variation of time to meet CO constrained (tau_sec) over tau_ent_cf')
    plt.xlabel('tau_ent_cf (ms)')
    plt.ylabel('tau_sec (ms)')
    for i in range(len(tau_ent_sec)):
        ax3.plot(tau_ent_cf, tau_sec[:, i]/milliseconds, label='tau_ent_sec = ' + str(tau_ent_sec[i]) + ' ms')
    ax3.grid(True)
    handles, labels = ax3.get_legend_handles_labels()
    lgd = ax3.legend(handles, labels, loc='upper right', bbox_to_anchor=(1.5,0.3))
    fig3.savefig('tau_sec_v_ent_cf.png', bbox_extra_artists=(lgd,), bbox_inches='tight')

    fig4 = plt.figure()
    ax4 = plt.axes()
    ax4.set_title('Variation of time to meet CO constrained (tau_sec) over tau_ent_sec')
    plt.xlabel('tau_ent_sec (ms)')
    plt.ylabel('tau_sec (ms)')
    for i in range(len(tau_ent_cf)):
        ax4.plot(tau_ent_sec, tau_sec[i, :]/milliseconds, label='tau_ent_cf = ' + str(tau_ent_cf[i]) + ' ms')
    ax4.grid(True)
    handles, labels = ax4.get_legend_handles_labels()
    lgd = ax4.legend(handles, labels, loc='upper right', bbox_to_anchor=(1.5,0.3))
    fig4.savefig('tau_sec_v_ent_sec.png', bbox_extra_artists=(lgd,), bbox_inches='tight')

    # plt.show()

if __name__ == "__main__":
    main()