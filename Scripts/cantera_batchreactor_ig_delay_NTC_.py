#notes
#gri vs klippenstien comparison of ignition times based on OH
# 25 atm, 1 atm, 10 atm
#range of phis (.5,.65,.75,1,5,200)
#temp; 300K, 591K
#edited from, as per cantera 2.40
#https://cantera.org/examples/jupyter/reactors/batch_reactor_ignition_delay_NTC.ipynb.html


import numpy as np
import matplotlib.pyplot as plt
import time
import cantera as ct

# Define the ignition delay time (IDT). This function computes the ignition
# delay from the occurrence of the peak concentration for the specified
# species.
def ignitionDelay(states, species):
    i_ign = states(species).Y.argmax()
    return states.t[i_ign]

def runcase(ctSol, reactorT, reactorP, phicase):

    # Define the reactor temperature and pressure
    gas = ct.Solution(ctSol)
    reactorTemperature = reactorT #in Kelvin
    reactorPressure = reactorP*101325 #in atm
    gas.TP = reactorTemperature, reactorPressure
    # Define the fuel, oxidizer and set the stoichiometry
    gas.set_equivalence_ratio(phi = phicase, fuel = 'ch4', oxidier={'02':1.0, 'n2':3.76})
    #creating reactor
    r = ct.Reactor(contents = gas)
    reactorNetwork = ct.ReactorNet([r])
    timeHistory = ct.SolutionArray(gas,extra=['t'])

    t0 = time.time()

    estimatedIgnitionDelayTime = 1.0
    t=0

    counter = 1
    while(t< estimatedIgnitionDelayTime):
        t = reactorNetwork.step()
        if(counter%10==0):
            timeHistory.append(r.thermo.state,t=t)
        counter+=1
    
    tau_ig = ignitionDelay(timeHistory, 'oh')

    t1 = time.time()
    print('Computed Ignition Delay: {:.3e} seconds. Took {:3.2f}s to compute'.format(tau_ig, t1-t0))
    plt.rcParams['axes.labelsize'] = 18
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12
    plt.rcParams['figure.autolayout'] = True

    plt.style.use('ggplot')
    plt.style.use('seaborn-pastel')

    plt.figure()
    plt.plot(timeHistory.index, timeHistory['oh'],'-o')
    plt.xlabel('Time (s)')
    plt.ylabel('$Y_{OH}$')

    plt.xlim([0,0.05])
    plt.arrow(0, 0.008, tau_ig, 0, width=0.0001, head_width=0.0005,
            head_length=0.001, length_includes_head=True, color='r', shape='full')
    plt.annotate(r'$Ignition Delay: \tau_{ign}$', xy=(0,0), xytext=(0.01, 0.0082), fontsize=16)
    plt.show()

    return tau_ig

def main():
    ctSol = 'gri30.xml'
    reactorT = 300
    reactorP = 1
    phicase = 0.5
    (tau_ig) = runcase(ctSol,reactorT,reactorP,phicase)
    print(tau_ig)
