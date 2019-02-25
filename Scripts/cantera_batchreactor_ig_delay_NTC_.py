#notes
#gri vs klippenstien comparison of ignition times based on OH
# 25 atm, 1 atm, 10 atm
#range of phis (.5,.65,.75,1,5,200)
#edited from, as per cantera 2.40
#https://cantera.org/examples/jupyter/reactors/batch_reactor_ignition_delay_NTC.ipynb.html

import sys
sys.path.append("..")
from CanteraTools import *
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time
import cantera as ct
ms = 0.001 #miliseconds
# Define the ignition delay time (IDT). This function computes the ignitions
# delay from the occurrence of the peak concentration for the specified
# species.
def ignitionDelay(states, species):
    return states[species].idxmax()

def runMainBurner(phi_main, tau_main, T_fuel=300, T_ox=650, P=25*101325, mech="gri30.xml", slope=0.01, curve=0.01): 
    flameGas = premix(phi_main, P=P, mech=mech, T_fuel=T_fuel, T_ox=T_ox) 
    # filename = '{0}_{1}-{2}_{3}-{4}_{5}'.format('phi_main', phi_main, 'P', P, )
    filename = '{0}_{1:.4f}.pickle'.format('phi_main', phi_main);
    if os.path.isfile(filename):
            mainBurnerDF = pd.read_parquet(filename)
            # mainBurnerDF = table.to_pandas()
            flameTime = mainBurnerDF.index.values;
    else:
        flame, flameTime = runFlame(flameGas, slope=slope, curve=curve)
        columnNames = ['x', 'u', 'T', 'n', 'MW'] + ["Y_" + sn for sn in flameGas.species_names] + ["X_" + sn for sn in
                                                                                    flameGas.species_names]
        flameData = np.concatenate(
[np.array([flame.grid]), np.array([flame.u]), np.array([flame.T]), np.array([[0] * len(flame.T)]), np.array([[0] * len(flame.T)]), flame.Y, flame.X], axis=0)
        mainBurnerDF = pd.DataFrame(data=flameData.transpose(), index=flameTime, columns=columnNames)
        mainBurnerDF.index.name = 'Time'
        mainBurnerDF['P'] = flame.P;
        # table = pa.Table.from_pandas(mainBurnerDF);
        # pq.write_table(table, filename);
        mainBurnerDF.to_parquet(filename, compression='gzip')
        
    vitiatedProd, flameCutoffIndex, mainBurnerDF = getStateAtTime(mainBurnerDF, flameTime, tau_main)
    vitReactor = ct.ConstPressureReactor(vitiatedProd)
    return vitReactor, mainBurnerDF, flameGas

def getMixtureTemp(phicase,ctSol):
    [vit_reactor, main_burner_DF, flameGas] = runMainBurner(0.3719, 19.842*ms,mech=ctSol)
    [mfm, mam, mfs, mas] = solvePhi_airSplit(phicase, 0.3719, 100, 1)
    main_mass = mfm + mam
    jet_mass = mfs + mas
    secondaryGas = ct.Solution('gri30.xml')
    secondaryGas.TPX = 300,1,{'CH4':1} 
    mixergas = mix([flameGas,secondaryGas],[main_mass, jet_mass])
    return mixergas.T

def runcase(ctSol, reactorT, reactorP, phicase):

    # Define the reactor temperature and pressure
    gas = ct.Solution(ctSol)
    reactorTemperature = reactorT #in Kelvin
    reactorPressure = reactorP*101325 #in atm
    gas.TP = reactorTemperature, reactorPressure

    # Define the fuel, oxidizer and set the stoichiometry
    #gas.set_equivalence_ratio(phi=phicase, fuel='CH4', oxidier={'02':1.0, 'N2':3.76})
    gas.set_equivalence_ratio(phi=1.0, fuel='CH4',
                              oxidizer={'o2':1.0, 'n2':3.76})
    #creating reactor
    r = ct.ConstPressureReactor(contents = gas)
    reactorNetwork = ct.ReactorNet([r])
    #timeHistory = ct.SolutionArray(gas, extra=['t'])

    stateVariableNames = [r.component_name(item) for item in range(r.n_vars)]
    timeHistory = pd.DataFrame(columns=stateVariableNames)
    #print(stateVariableNames)

    t0 = time.time()

    estimatedIgnitionDelayTime = 1000.0*ms
    t = np.arange(0,estimatedIgnitionDelayTime,0.1*ms)
    

    
    for i in range(len(t)):
        tnow = t[i]
        reactorNetwork.advance(tnow)
        timeHistory.loc[tnow] = reactorNetwork.get_state()
    #print(timeHistory['OH'])

    tau_ig = ignitionDelay(timeHistory,'OH')

    t1 = time.time()
    runtime = t1-t0
    print('Computed Ignition Delay: {:.3e} seconds. Took {:3.2f}s to compute'.format(tau_ig, runtime))
    # plt.rcParams['axes.labelsize'] = 18
    # plt.rcParams['xtick.labelsize'] = 12
    # plt.rcParams['ytick.labelsize'] = 12
    # plt.rcParams['figure.autolayout'] = True

    # plt.style.use('ggplot')
    # plt.style.use('seaborn-pastel')

    # plt.figure()
    # plt.plot(timeHistory.index, timeHistory['OH'],'-o')
    # plt.xlabel('Time (s)')
    # plt.ylabel('$Y_{OH}$')

    # plt.xlim([0,0.005])
    # plt.arrow(0, 0.008, tau_ig, 0, width=0.0001, head_width=0.0005,
    #         head_length=0.001, length_includes_head=True, color='r', shape='full')
    # plt.annotate(r'$Ignition Delay: \tau_{ign}$', xy=(0,0), xytext=(0.01, 0.0082), fontsize=16)
    # plt.show()

    return tau_ig,runtime

def main():
    ct.suppress_thermo_warnings()
    ctSol = ('gri30.cti' ,'Klippenstein.cti') 
    reactorP = (10,25)
    phicase = (.5,.65,.75,1,5,200)
    lengthy = len(reactorP)*len(phicase)
    for x in range(len(ctSol)):
        counter = 0
        dataarray = np.zeros((lengthy,5))
        for y in range(len(reactorP)):
            for z in range(len(phicase)):
                reactorT = getMixtureTemp(phicase[z],ctSol[x])
                print(reactorT)
                print(phicase[z])
                (tau_ig,runtime) = runcase(ctSol[x],reactorT,reactorP[y],phicase[z])
                dataarray[counter,:] = [reactorT,reactorP[y],phicase[z],tau_ig,runtime]
                counter+=1
        pd.DataFrame(dataarray).to_csv('ignitionDelay%s.csv' %(str(x)),header=None,index=None)

if __name__ == "__main__":
    main()