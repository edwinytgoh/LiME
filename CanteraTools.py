import cantera as ct 
import numpy as np 
import pandas as pd
import time 
import pdb 
from argparse import ArgumentParser
import os.path
import pyarrow.parquet as pq 
import pyarrow as pa
import multiprocessing

milliseconds = 0.001 # seconds 
pd.options.mode.chained_assignment = None  # default='warn'
def runFlame(gas, slope=0.01, curve=0.01):
    # Simulation parameters
    width = 0.06  # m
    # Flame object
    f = ct.FreeFlame(gas, width=width)
    f.set_refine_criteria(ratio=2, slope=slope, curve=curve)
    f.transport_model = 'Multi'
    f.solve(loglevel=0, auto=True, refine_grid=True)
    

    # Convert distance into time:
    CH2O = gas.species_index('CH2O');
    X_CH2O = f.X[CH2O]
    maxIndex = np.arange(0, len(X_CH2O))[X_CH2O == max(X_CH2O)][0];
#     startingIndex = np.arange(0, len(X_CH2O))[X_CH2O >= X_CH2O[0] + 5][0]
    startingIndex = maxIndex; 
    #     startingIndex = np.arange(0, len(f.heat_release_rate))[f.heat_release_rate == max(f.heat_release_rate)][0]
    u_avg = np.array(f.u[startingIndex:] + f.u[startingIndex - 1:-1]) * 0.5
    dx = np.array(f.grid[startingIndex:] - f.grid[startingIndex - 1:-1])
    dt = dx / u_avg

    pre_time = [-999] * (startingIndex-1)  # need to add some padding to the time array to account for dist = 0.0
    pre_time.extend([0])
    time = np.hstack(
        (np.array(pre_time), np.cumsum(dt)))  # numerically integrate (read: sum) dt to get array of times for each x location

    return f, time

def getStateAtTime(flame, tList, tRequired, mech='gri30.xml'):
    '''A function that gets the state at a desired point in time of a flame simulation performed using runFlame. 
        Takes in a flame object, its associated time series, and the desired point in time (in seconds).
        Returns a new Cantera gas object with the desired state, and the corresponding index in the flame at which the point was found. 

        Example usage: gas, t_ind = getStateAtTime(flame, time, t_req)'''
    mainBurnerDF = flame
    columnNames = mainBurnerDF.columns.values
    vel_final = mainBurnerDF['u'].iloc[-1]
    moleFracs = mainBurnerDF.columns[mainBurnerDF.columns.str.contains('X_')]
    assert (tRequired > 0)
    newGas = ct.Solution(mech)

    if (tRequired >= max(tList)):
        tau_vit = tRequired - max(tList)
        newGas.TPX = flame['T'].iloc[-1], flame['P'].iloc[-1], flame[moleFracs].iloc[-1]
        tau_vit_start = tList[-1]
    else:
        dt_list = abs(tList - tRequired)
        # Find index at which required time is closest to an entry in tList:     
        minIndex = next(ind for ind, val in enumerate(dt_list) if abs(val - min(dt_list)) < 1e-6)
        flameTime_closestTreq = tList[minIndex] 
        if ((flameTime_closestTreq - tRequired) > 1e-6):
            newGas.TPX = flame['T'].iloc[minIndex-1], flame['P'].iloc[minIndex-1], flame[moleFracs].iloc[minIndex-1]
            assert((newGas['NO'].X - flame['X_NO'].iloc[minIndex-1] <= 1e-6))    
            tau_vit = tRequired - tList[minIndex - 1]
            tau_vit_start = tList[minIndex - 1]
            assert(tau_vit > 0)# if the closest flameTime is larger than tRequired, the second closest should be less than, otherwise /that/ would be the closest...?
        elif ((tRequired - flameTime_closestTreq) > 1e-6): # if closest time is more than 1e-6 s less than tReq
            newGas.TPX = flame['T'].iloc[minIndex], flame['P'].iloc[minIndex], flame[moleFracs].iloc[minIndex] 
            assert((newGas['NO'].X - flame['X_NO'].iloc[minIndex] <= 1e-6))    
            tau_vit = tRequired - tList[minIndex] 
            tau_vit_start = tList[minIndex]
        else:
            newGas.TPX = flame['T'].iloc[minIndex], flame['P'].iloc[minIndex], flame[moleFracs].iloc[minIndex]
            tau_vit = 0
            assert((newGas['NO'].X - flame['X_NO'].iloc[minIndex] <= 1e-6))    
    if tau_vit > 0:
        vitiator = ct.ConstPressureReactor(newGas)
        vitRN = ct.ReactorNet([vitiator])
        dt = 0.0001 * 1e-3
        vit_tList = np.arange(0, tau_vit, dt)
        vitArray = np.array([None] * len(vit_tList) * len(columnNames)).reshape(len(vit_tList), len(columnNames))
        # performance note: we don't need to do this. we can simply just advance to the desired tau_main
        for i in range(0, len(vit_tList)):
            MWi = vitiator.thermo.mean_molecular_weight
            vitArray[i, :] = np.hstack([vel_final * dt, vel_final, vitiator.thermo.T, 0, MWi, vitiator.thermo.Y, vitiator.thermo.X, vitiator.thermo.P])
            vitRN.advance(vit_tList[i])
        vit_tList += tau_vit_start             
        vitDF = pd.DataFrame(data=vitArray, index=vit_tList, columns=columnNames)
        print("Vitiator end time:", vit_tList[-1]/1e-3, "milliseconds") 
        vitiator.syncState()
        '''Call syncState so that newGas has the right state for use in later functions.'''
        minIndex = -1  # last index in time list 
        mainBurnerDF = mainBurnerDF[np.array(mainBurnerDF.index > 0) & np.array(mainBurnerDF.index <= tRequired)]
        print("Initial mainBurnerDF length:", len(mainBurnerDF.index.values))
        mainBurnerDF = pd.concat([mainBurnerDF, vitDF])        
        print("New mainBurnerDF length:", len(mainBurnerDF.index.values))
    else:
        mainBurnerDF = mainBurnerDF[np.array(mainBurnerDF.index > 0) & np.array(mainBurnerDF.index <= tRequired)]
    mainBurnerDF['NOppmvd'] = correctNOx(mainBurnerDF['X_NO'], mainBurnerDF['X_H2O'], mainBurnerDF['X_O2'])
    mainBurnerDF['COppmvd'] = correctNOx(mainBurnerDF['X_CO'], mainBurnerDF['X_H2O'], mainBurnerDF['X_O2'])    
    return newGas, minIndex, mainBurnerDF

def fstoich(fuel={'CH4':1}, ox={'O2':0.21, 'N2':0.79}, mech='gri30.xml'): 
    '''Function that returns the stoichiometric fuel-to-air ratio (by mass) for a given fuel and oxidizer. 
    
    Example usage: f_stoich_ch4_air = fstoich() 
                    f_stoich_c2h4_o2 = fstoich(fuel={'C2H4':1}, ox={'O2':1})'''
    gas = ct.Solution(mech) 
    gas.set_equivalence_ratio(1.0, fuel, ox) 
    return sum(gas[fuel.keys()].Y)/sum(gas[ox.keys()].Y)

def solvePhi_airSplit(phiGlobal, phiMain, mdotTotal=1000, airSplit=1):
    # fs = 
    fs = 0.058387057492574147;
    mfm = airSplit*fs*mdotTotal*(1+fs*phiGlobal)**(-1)*phiMain
    mam = airSplit*mdotTotal*(1+fs*phiGlobal)**(-1)
    mfs = (-1)*(1+fs*phiGlobal)**(-1)*((-1)*fs*mdotTotal*phiGlobal+airSplit*fs*mdotTotal*phiMain)
    mas = (-1)*((-1)+airSplit)*mdotTotal*(1+fs*phiGlobal)**(-1)
    return mfm, mam, mfs, mas

def solveMass_PhiJet(phiGlobal, phiMain, phiJet, mdotTotal=1000):
    fs = 0.058387057492574147
    mam = 1
    mfm = mam*fs*phiMain
    mfs1 = mam*fs*(phiGlobal - phiMain)
    mas = mfs1/(fs*(phiJet-phiGlobal))
    mfs = mfs1 + mas*fs*phiGlobal
    # Normalize
    mtot = mam + mfm + mas + mfs
    mam *= mdotTotal/mtot
    mfm *= mdotTotal/mtot
    mas *= mdotTotal/mtot
    mfs *= mdotTotal/mtot
    return mfm, mam, mfs, mas

def correctNOx(X_i, X_H2O, X_O2):
    dry_i = X_i/(1 - X_H2O)
    dry_O2 = X_O2/(1 - X_H2O)
    corrected = dry_i*(20.9 - 15)/(20.9 - dry_O2*100)
    corrected_ppm = corrected*1e6
    return corrected_ppm     

def mix(streams, mdots, mech="gri30.xml", P=25*101325):
    # Create mixer gas: 
    mixerGas = ct.Solution(mech) 
    mixerGas.TPX = [300, P, 'H2:1']
    
    # Create reactor with CHEMISTRY DISABLED: 
    mixer = ct.ConstPressureReactor(mixerGas) 
    mixer.chemistry_enabled = False # disable chemistry 
    
    # For each stream (and given mass flow rate), connect mass flow controller to mixer: 
    mfcs = [ct.MassFlowController(ct.Reservoir(streams[i]), mixer, mdot=mdots[i]) for i in range(0,len(streams))]
    
    exhaust = ct.Reservoir(mixerGas) 
    exhaust_mfc = ct.MassFlowController(mixer, exhaust, mdot=sum(mdots)) 
    
    rn = ct.ReactorNet([mixer]) 
    rn.advance_to_steady_state() 
    
    return mixer.thermo

def premix(phi=0.4, fuel={'CH4':1}, ox={'N2':0.79, 'O2':0.21}, mech='gri30.xml', P=25*101325, T_fuel=300, T_ox=650, M_total=1):
    
    air = ct.Solution(mech)
    air.TPX = [T_ox, P, ox]
    fuelGas = ct.Solution(mech)
    fuelGas.TPX = T_fuel, P, fuel
    
    # Temporary ThermoPhase object to get mass flow rates:
    temp = ct.Solution('gri30.xml')
    temp.set_equivalence_ratio(phi, fuel, ox)
    mdot_fuel = M_total * sum(temp[fuel.keys()].Y)
    mdot_ox = M_total * sum(temp[ox.keys()].Y)

    # Output mixer gas: 
    return mix([fuelGas, air], [mdot_fuel, mdot_ox], P=P)       

def iem(m, tpArray, rArray, rn, dt, omega):
    # Constant k:
    C_phi = 2
    k = -C_phi * omega * 0.5 * dt
    # Calculate average: 
    m_total_r = 1/sum(m)
    M_species_total = sum([m[i] * tpArray[i].Y for i in range(0, len(tpArray))]) 
    H_total = sum([m[i]*tpArray[i].enthalpy_mass for i in range(0, len(tpArray))])
    Y_avg = M_species_total * m_total_r # Y_species_avg = (M_total_species)/(M_total_system)
    h_avg = H_total * m_total_r # H_avg is the specific mass-weighted average across all reactors of the total enthalpy.
    # Adjust reactor state:     
    for i in range(0, len(tpArray)):
        Y_current = tpArray[i].Y
        Y_new =  Y_current + k * (Y_current - Y_avg)         
        h = tpArray[i].enthalpy_mass
        h_new = h + k * (h - h_avg)
        tpArray[i].HPY = [h_new, tpArray[i].P, Y_new]
        rArray[i].syncState()
    # Reinitialize reactor network solver:         
    rn.reinitialize()    
    return Y_avg, h_avg

def runMainBurner(phi_main, tau_main, T_fuel=300, T_ox=650, P=25*101325, mech="gri30.xml", slope=0.01, curve=0.01): 
    flameGas = premix(phi_main, P=P, mech=mech, T_fuel=T_fuel, T_ox=T_ox) 
    # filename = '{0}_{1}-{2}_{3}-{4}_{5}'.format('phi_main', phi_main, 'P', P, )
    filename = '{0}_{1:.4f}.pickle'.format('phi_main', phi_main);
    if os.path.isfile(filename):
            table = pq.read_table(filename, nthreads=5)
            mainBurnerDF = table.to_pandas()
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
        table = pa.Table.from_pandas(mainBurnerDF);
        pq.write_table(table, filename);
        
    vitiatedProd, flameCutoffIndex, mainBurnerDF = getStateAtTime(mainBurnerDF, flameTime, tau_main)
    vitReactor = ct.ConstPressureReactor(vitiatedProd)
    return vitReactor, mainBurnerDF

def dataFrame_to_pyarrow(df, filename):
    pq.write_table(pa.Table.from_pandas(df), filename)

def pyarrow_to_dataFrame(filename):
    table = pq.read_table(filename, nthreads=4)
    df = table.to_pandas()
    return df

def twoStage_ideal(phi_global,phi_main,tau_global,tau_sec,airSplit=1,phiSec=None,T_fuel=300,T_ox=650,P=25, mech="gri30.xml",trace=False):
    P *= ct.one_atm
    tau_global *= milliseconds 
    tau_sec *= milliseconds 
    tau_main = tau_global - max(0,tau_sec)
    mfm, mam, mfs, mas = solvePhi_airSplit(phi_global, phi_main, airSplit=airSplit)
    if ((mfm < 0) | (mfs < 0) | (mas < 0) | (mam < 0)):
        return np.hstack([0, 0, 100000, 100000, 10000, 10000, 0, 0, 0, np.full(53, 0, dtype=np.float64), np.full(53, 0, dtype=np.float64), 100000, 100000, 100000])
    m = [mfm + mam, mfs + mas] 
    mainMassPerc = 100*(mfm + mam)/(mfm + mam + mfs + mas) 
    secMassPerc = 100*(mfs + mas)/(mfm + mam + mfs + mas)
    condColumns = ['phi_global', 'phi_main', 'tau_global', 'tau_main', 'tau_sec', 'tau_mix', 'dt', 'T_fuel', 'T_ox', 'main_mass_perc', 'sec_mass_perc', 'mfm', 'mam', 'mfs', 'mas', 'airSplit']

    if (secMassPerc <= 0.001):
        tau_main = tau_global;
        tau_sec = 0;
    if (tau_sec <= 0):
        phi_main = phi_global;

    # RUN MAIN BURNER
    tic = time.time() 

    # pdb.set_trace()
    vitReactor, mainBurnerDF = runMainBurner(phi_main, tau_main, T_fuel, T_ox, P=P, mech=mech)
    massFracs = mainBurnerDF.columns[mainBurnerDF.columns.str.contains('Y_')] 
    moleFracs = mainBurnerDF.columns[mainBurnerDF.columns.str.contains('X_')]       
    flameTime = mainBurnerDF.index.values

    toc = time.time() 
    # print("Time taken to get flame:", (toc-tic), "secs")
    if (tau_sec <= 0):
        T = mainBurnerDF['T'].iloc[-1]
        Y_values = mainBurnerDF[massFracs].iloc[-1]

        NOCOppmvd = correctNOx(np.array(mainBurnerDF[['X_NO', 'X_CO', 'X_NO2', 'X_N2O']].iloc[-1], dtype=np.float64), mainBurnerDF['X_H2O'].iloc[-1], mainBurnerDF['X_O2'].iloc[-1])     
        if (trace==True): 
            mainDF = mainBurnerDF[['T', 'NOppmvd', 'COppmvd']];
            mainDF.index.name = 'time';
            return mainDF;
        return np.hstack([flameTime[-1], mainBurnerDF['T'].iloc[-1], NOCOppmvd[0], NOCOppmvd[1], NOCOppmvd[2], NOCOppmvd[3], np.array(mainBurnerDF.iloc[-1],dtype=np.float64)])
    # CREATE FUEL PARTICLE
    fuelStream = ct.Solution(mech) 
    oxStream = ct.Solution(mech) 
    fuelStream.TPX = [T_fuel, P, {'CH4':1.0}]
    oxStream.TPX = [T_ox, P, {'O2':0.21, 'N2':0.79}]     
    secondaryReactor = ct.ConstPressureReactor(mix([fuelStream, oxStream, vitReactor.thermo], [mfs, mas, mfm + mam], mech=mech)) 
    tpArray = np.array([vitReactor.thermo, secondaryReactor.thermo])
    rArray = np.array([vitReactor, secondaryReactor])

    rn = ct.ReactorNet([secondaryReactor])
    if (trace == False):
        # RUN SECONDARY REACTOR
        rn.advance(tau_sec) 

        # GET ENDPOINT DATA
        NOCOppmvd = correctNOx(secondaryReactor.thermo['NO', 'CO', 'NO2', 'N2O'].X, secondaryReactor.thermo['H2O'].X, secondaryReactor.thermo['O2'].X) 
        out = np.hstack([rn.time, secondaryReactor.thermo.T, NOCOppmvd[0], NOCOppmvd[1], NOCOppmvd[2], NOCOppmvd[3], secondaryReactor.thermo.T, 0, secondaryReactor.thermo.mean_molecular_weight, secondaryReactor.thermo.Y, secondaryReactor.thermo.X, secondaryReactor.thermo.P, NOCOppmvd[0], NOCOppmvd[1]])   
        # out = np.hstack([phi_main, tau_sec, NOCOppmvd[0], NOCOppmvd[1], secondaryReactor.thermo.T])
        return out
    else:
        dt = 0.00001 * 1e-3
        sec_tList = np.arange(0, tau_sec, dt)
        columnNames = ['time', 'T', 'NOppmvd', 'COppmvd']
        secArray = np.array([None] * len(sec_tList) * len(columnNames)).reshape(len(sec_tList), len(columnNames))
        # performance note: we don't need to do this. we can simply just advance to the desired tau_main
        mainDF = mainBurnerDF[['T', 'NOppmvd', 'COppmvd']];
        mainDF.index.name = 'time';
        for i in range(0, len(sec_tList)):
            MWi = secondaryReactor.thermo.mean_molecular_weight
            NOCOppmvd = correctNOx(secondaryReactor.thermo['NO', 'CO', 'NO2', 'N2O'].X, secondaryReactor.thermo['H2O'].X, secondaryReactor.thermo['O2'].X) 
            # secArray[i, :] = np.hstack([vel_final * dt, vel_final, secondaryReactor.thermo.T, 0, MWi, secondaryReactor.thermo.Y, secondaryReactor.thermo.X, secondaryReactor.thermo.P])
            rn.advance(sec_tList[i])
            secArray[i,:] = np.hstack([sec_tList[i] + tau_main, secondaryReactor.thermo.T, NOCOppmvd[0], NOCOppmvd[1]])
        sec_tList += tau_main             
        secDF = pd.DataFrame(data=secArray, columns=columnNames)
        secDF = secDF.set_index('time');
        secondaryReactor.syncState() 
        stagedDF = pd.concat([mainDF, secDF]) 
        return stagedDF;

def get_tauHot(timeSeries, out_df, dT = 200): ## REMEMBER TO CHECK UNITS!!!
    """Function that calculates tau_hot, i.e. timescale that represents 
    the duration of high-temperature regions in a staged combustor's secondary stage.

    Parameters
    ----------
    timeSeries : `pandas.DataFrame`
        DataFrame that contains the time history of the current case. Needs to contain columns ['age'] and ['T'] for current time and corresponding temperature, respectively.
    
    out_df : `pandas.DataFrame`
        DataFrame of shape (1,n) that contains parameters of the current case: ['tau_ign_OH'], ['tau_sec_required'], ['T']. Warning: DO NOT pass in a pandas.Series!

    dT : float
        Number that represents drop in temperature allowed before calling the end of the high-temperature region 
    
    Returns
    --------
    tau_hot : float
        Calculated tau_hot in milliseconds 

    tau_hot_start : float 
        Starting time of the high-temperature region in seconds 
    
    tau_hot_end : float 
        Ending time of the high-temperature region in seconds

    """
    tau_hot_start = out_df['tau_ign_OH'].values[0]
    T_max = max(timeSeries['T'])
    T_final = out_df['T'].values
    # print(f"T_max = {T_max:.2f} K;\nIgnition delay based on OH conc: {tau_hot_start/1e-3} ms")
    overshoot = T_max > T_final + 15
    overshoot_value = T_max - T_final
    max_ind = timeSeries['T'].idxmax()
    if overshoot:
        # print(f"Overshoot by: {T_max - out_df['T']:.2f} K")
        # find temperature after peak where T = T_max - dT
        remaining_df = timeSeries.iloc[max_ind:]
        if overshoot_value > dT:
            remaining_df = remaining_df[(remaining_df['T'] <= T_max - dT)]
            tau_hot_end = remaining_df['age'].values[0]
        elif overshoot_value > 0.5*dT:
            remaining_df = remaining_df[(remaining_df['T'] <= T_max - 0.5*dT)] # In this case, T_max - T_final > 0.5*dT, T_max - 0.5*dT > T_final
            tau_hot_end = remaining_df['age'].values[0]
        else:
            tau_hot_end = out_df['tau_sec_required'].values[0]*1e-3
    else: 
        # print(f"No overshoot; using tau_CO = {out_df['tau_sec_required'].iloc[0]:.2} ms")
        tau_hot_end = out_df['tau_sec_required'].values[0]*1e-3
    tau_hot = (tau_hot_end - tau_hot_start)/1e-3
    return tau_hot, tau_hot_start, tau_hot_end

def get_tauNOx(timeSeries, tauHot_start, tauHot_end, P=25*101325):
    """Function that calculates tau_NOx, i.e. timescale that represents 
    the rate of NOx production in a staged combustor's high-temperature region.

    Parameters
    ----------
    timeSeries : `pandas.DataFrame`
        DataFrame that contains the time history of the current case. Needs to contain columns ['age'] and ['T'] for current time and corresponding temperature, respectively.
    
    tauHot_start : float 
        Starting time of the high-temperature region in seconds 
    
    tauHot_end : float 
        Ending time of the high-temperature region in seconds
    
    P : float 
        Pressure of the system in Pascals. Defaults to 25 atm, i.e. 2.5 MPa

    Returns
    --------
    tau_NOx : float
        Calculated tau_NOx in milliseconds 

    """    
    eps = 0.0001*1e-3
    ign_state = timeSeries[(timeSeries['age'] >= tauHot_start - eps) & (timeSeries['age'] <= tauHot_start + eps)] # system state at ignition
    end_state = timeSeries[(timeSeries['age'] >= tauHot_end - eps) & (timeSeries['age'] <= tauHot_end + eps)] # system "end" state

    R_universal = 8.3144598; # J/mol-K or m3-Pa/mol-K
    M = (P/R_universal)/timeSeries['T'] # units of moles/volume
    M = (M[1:].values + M[:-1].values)*0.5

    
    d_XNO = timeSeries['X_NO'].iloc[1:].values - timeSeries['X_NO'].iloc[:-1].values
    d_t = timeSeries['age'].iloc[1:].values - timeSeries['age'].iloc[:-1].values
    dNO_dt = (M*d_XNO)/d_t
    start_ind = ign_state.index.values[0]
    end_ind = end_state.index.values[0] + 1
    tau_NOx = np.mean(M[start_ind:end_ind]/(dNO_dt[start_ind:end_ind]))
    return tau_NOx/1e-3