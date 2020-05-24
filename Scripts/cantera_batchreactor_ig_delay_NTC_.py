#gri vs klippenstien comparison of ignition times based on OH
# 25 atm, 1 atm, 10 atm
#range of phis (.5,.65,.75,1,5,200)
#edited from, as per cantera 2.40
#https://cantera.org/examples/jupyter/reactors/batch_reactor_ignition_delay_NTC.ipynb.html

import sys
sys.path.append("../")
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

# def getMixtureTemp(phicase,mech):
#     [vit_reactor, main_burner_DF] = runMainBurner(0.3719, 19.842*ms,mech=mech)
#     [mfm, mam, mfs, mas] = solvePhi_airSplit(0.675, 0.3719, 100, 1)
#     main_mass = mfm + mam
#     jet_mass = mfs + mas
#     secondaryGas = ct.Solution('gri30.xml')
#     secondaryGas.TPX = 300,1,{'CH4':1} 
#     mixergas = mix([vit_reactor.thermo,secondaryGas],[main_mass, jet_mass], mech=mech)
#     return mixergas.T

#takes in Mechanism name, reactor TP, Phi, timestep you want to go through,
#outputs tau_ignition, runtime,Temp at tau_ig, OH at tau_ig, Temp max of Reactor through time
#outputs csv of timehistory contain Temp, OH, CO2, NO2, etc with phi,reactorP,Mechname in title
def runcase(mech, reactorT, reactorP, phicase,timestep):

    # Define the reactor temperature and pressure
    gas = ct.Solution(mech)
    reactorTemperature = reactorT #in Kelvin
    reactorPressure = reactorP*101325 #in atm
    gas.TP = reactorTemperature, reactorPressure

    # Define the fuel, oxidizer and set the stoichiometry
    gas.set_equivalence_ratio(phi=phicase, fuel='CH4',
                              oxidizer={'o2':1.0, 'n2':3.76})
    #creating reactor
    r = ct.ConstPressureReactor(contents = gas)
    reactorNetwork = ct.ReactorNet([r])

    #initializing Dataframe to store data
    stateVariableNames = [r.component_name(item) for item in range(r.n_vars)] + ['T']
    timeHistory = pd.DataFrame(columns=stateVariableNames)

    t0 = time.time()
    runtime = 0
    times = 0;
    #stepping through to our runtimes;
    while runtime < 3:
        times += timestep
        reactorNetwork.advance(times)
        timeHistory.loc[times] = np.hstack([reactorNetwork.get_state(), r.thermo.T])
        t1 = time.time()
        runtime = t1-t0
    #     if timeHistory.loc[times,'T'] > 2000:
    #         tau_ig = timeHistory['OH'].idxmax()
    #         break
    #     else:
    #         tau_ig = 9999
    #         pass
    # #print(timeHistory['OH'])
    tau_ig = timeHistory['OH'].idxmax()
    tatig = timeHistory.loc[tau_ig,'T']
    OHatT = timeHistory.loc[tau_ig,'OH']
    TMax = timeHistory['T'].max()
    NOMax = timeHistory['NO'].max()
    #if reactor did not increase in temperature, it never ignited, ie tau_ig = inf, or 0
    if TMax <= 1620:
        tau_ig = 999
    else:
        pass

        
    ##some plots
    # print('Computed Ignition Delay: {:.3e} seconds. Took {:3.2f}s to compute'.format(tau_ig, runtime))
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
    #filename = input('Please Enter Filename')
    filename = ('timehistory{}{}{}.csv'.format(reactorP,phicase,mech[:mech.rfind(".")]))
    timeHistory.to_csv(filename);

    return tau_ig,runtime,OHatT,TMax,NOMax

#saves csv with reatorT, Pressures, Phi, Tau_ig, runtime, 
def main():
    ct.suppress_thermo_warnings()
    mech = ('gri30.cti','Klippenstein.cti' ) 
    reactorP = (1,10,25) #in atm
    phicase = np.arange(1,50,5)
    lengthy = len(reactorP)*len(phicase)
    reactorT = 1600 #in Kelvin
    for x in range(len(mech)):
        counter = 0
        dataarray = np.zeros((lengthy,8))
        for y in range(len(reactorP)):
            for z in range(len(phicase)):
                #reactorT = getMixtureTemp(phicase[z],mech[x])
                (tau_ig,runtime,OHatT,TMax,NOMax) = runcase(mech[x],reactorT,reactorP[y],phicase[z],.001*milliseconds)
                dataarray[counter,:] = [reactorT,reactorP[y],phicase[z],tau_ig,runtime,OHatT,TMax,NOMax]
                counter+=1
        pd.DataFrame(dataarray).to_csv('ignitiondelay{}.csv'.format(mech[x][:mech[x].rfind(".")]),header=None,index=None)


def main_edwin(): 
    ct.suppress_thermo_warnings()
    T_list = np.linspace(800,1000,4); 
    tau_ign_gri_list = [None]*len(T_list)
    runtime_gri_list = [None]*len(T_list)
    tau_ign_klip_list = [None]*len(T_list)
    runtime_klip_list = [None]*len(T_list)
    gas_gri = ct.Solution("gri30.xml"); 
    gas_klip = ct.Solution("Klippenstein.cti"); 
    dataarrayg = np.zeros((len(T_list),3))
    dataarrayi = np.zeros((len(T_list),3))
    counter = 0
    for i,T in enumerate(T_list):
        tau_ign_gri_list[i], runtime_gri_list[i] = runcase("gri30.xml", T, 25, 0.6)
        print(f"Done running GRI case T = {T:.3f} K. Time taken = {runtime_gri_list[i]:.3f} seconds");
        dataarrayg[counter,:] = [T_list[i],tau_ign_gri_list[i], runtime_gri_list[i]]
        tau_ign_klip_list[i], runtime_klip_list[i] = runcase("Klippenstein.cti", T, 25, 0.6)
        print(f"Done running Klippenstein case T = {T:.3f} K. Time taken = {runtime_klip_list[i]:.3f} seconds");
        dataarrayi[counter,:] = [T_list[i],tau_ign_klip_list[i], runtime_klip_list[i]]
        counter += 1
    pd.DataFrame(dataarrayg).to_csv('ignitionDelayGri.csv')
    pd.DataFrame(dataarrayi).to_csv('ignitionDelayKlip.csv')
    fig = plt.figure();
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    ax1.plot(T_list,tau_ign_gri_list)
    ax2.plot(T_list, runtime_gri_list)
    ax3.plot(T_list,tau_ign_klip_list)
    ax4.plot(T_list, runtime_klip_list)
    plt.show()

# if __name__ == "__main__":
#     main_edwin() #head_width=0.0005,
    #         head_length=0.001, length_includes_head=True, color='r', shape='full')
    # plt.annotate(r'$Ignition Delay: \tau_{ign}$', xy=(0,0), xytext=(0.01, 0.0082), fontsize=16)
    #plt.show()

    #return tau_ig,runtime


if __name__ == "__main__":
    main()