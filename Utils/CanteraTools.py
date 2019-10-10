import cantera as ct 
import numpy as np 
import pandas as pd
import time 
import pdb 
from argparse import ArgumentParser
import os.path
from functools import partial
from numba import jit
# import pyarrow.parquet as pq 
# import pyarrow as pa
import multiprocessing as mp
fs = 0.058387057492574147288255659304923028685152530670166015625;
milliseconds = 0.001 # seconds 
pd.options.mode.chained_assignment = None  # default='warn'

def air(T:float=300, P=25*ct.one_atm, mech='gri30.xml'):
    air = ct.Solution(mech)
    air.TPX = T, P, {'N2':0.79, 'O2':0.21}
    return air

def CH4(T:float=300,P=25*ct.one_atm,mech='gri30.xml'):
    CH4 = ct.Solution(mech)
    CH4.TPX = T, P, {'CH4':1.0}
    return CH4

def runFlame(gas, slope=0.01, curve=0.01):
    # Simulation parameters
    width = 0.06  # m
    # Flame object
    f = ct.FreeFlame(gas, width=width)
    f.set_refine_criteria(ratio=2, slope=slope, curve=curve, prune=min(slope,curve)/1e3)
    f.transport_model = 'Multi'
    f.soret_enabled = True
    f.max_grid_points = 10000
    f.solve(loglevel=1, auto=True, refine_grid=True)
    # Convert distance into time:
    CH2O = gas.species_index('CH2O');
    X_CH2O = f.X[CH2O]
    maxIndex = np.arange(0, len(X_CH2O))[X_CH2O == max(X_CH2O)][0];
    maxIndex2 = X_CH2O.argmax()
    assert maxIndex == maxIndex2
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
        # print("Vitiator end time:", vit_tList[-1]/1e-3, "milliseconds") 
        vitiator.syncState()
        '''Call syncState so that newGas has the right state for use in later functions.'''
        minIndex = -1  # last index in time list 
        mainBurnerDF = mainBurnerDF[np.array(mainBurnerDF.index > 0) & np.array(mainBurnerDF.index <= tRequired)]
        # print("Initial mainBurnerDF length:", len(mainBurnerDF.index.values))
        mainBurnerDF = pd.concat([mainBurnerDF, vitDF])        
        # print("New mainBurnerDF length:", len(mainBurnerDF.index.values))
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

def calculate_flowRates(phiGlobal, phiMain, phiSec):
    fs = 0.058387057492574147288255659304923028685152530670166015625;
    tol = -1e-14
    mfm = 1
    mam = 1/(phiMain*fs)
    mas = (phiGlobal/phiMain - 1)/((phiSec-phiGlobal)*fs)
    mfs = phiGlobal/phiMain * (1 + (phiGlobal-phiMain)/(phiSec - phiGlobal)) - 1

    # assert mfm >= tol and mfs >= tol and mam >= tol and mas >= tol, f"You have negative mass flow rates! mfm = {mfm}, mam = {mam}, mfs = {mfs}, mas = {mas}"

    m_perc = 100/(mfm + mam + mas + mfs)
    # return {'mfm':mfm*m_perc, 'mam':mam*m_perc, 'mfs':mfs*m_perc, 'mas':mas*m_perc}
    return mfm*m_perc, mam*m_perc, mfs*m_perc, mas*m_perc
def solvePhi_airSplit(phiGlobal, phiMain, mdotTotal=1000, airSplit=1):
    # fs = 
    fs = 0.058387057492574147288255659304923028685152530670166015625;
    mfm = airSplit*fs*mdotTotal*(1+fs*phiGlobal)**(-1)*phiMain
    mam = airSplit*mdotTotal*(1+fs*phiGlobal)**(-1)
    mfs = (-1)*(1+fs*phiGlobal)**(-1)*((-1)*fs*mdotTotal*phiGlobal+airSplit*fs*mdotTotal*phiMain)
    mas = (-1)*((-1)+airSplit)*mdotTotal*(1+fs*phiGlobal)**(-1)
    return mfm, mam, mfs, mas

def solveMass_PhiJet(phiGlobal, phiMain, phiJet, mdotTotal=1000):
    fs = 0.058387057492574147288255659304923028685152530670166015625
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

def runMainBurner(phi_main, tau_main, T_fuel=300, T_ox=650, P=25*101325, mech="gri30.xml", slope=0.01, curve=0.01, filename=None): 
    flameGas = premix(phi_main, P=P, mech=mech, T_fuel=T_fuel, T_ox=T_ox) 
    # filename = '{0}_{1}-{2}_{3}-{4}_{5}'.format('phi_main', phi_main, 'P', P, )
    if filename == None:
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
    vitReactor = ct.ConstPressureReactor(vitiatedProd, name='MainBurner')
    return vitReactor, mainBurnerDF

def dataFrame_to_pyarrow(df, filename):
    # pq.write_table(pa.Table.from_pandas(df), filename)
    df.to_parquet(filename, compression='gzip')

def pyarrow_to_dataFrame(filename):
    # table = pq.read_table(filename, nthreads=4)
    # df = table.to_pandas()
    df = pd.read_parquet(filename)
    return df


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
    # pdb.set_trace()
    tau_hot_start = out_df['tau_ign_OH'].values[0] # seconds
    T_max = max(timeSeries['T']) 
    T_final = out_df['T'].values
    # print(f"T_max = {T_max:.2f} K;\nIgnition delay based on OH conc: {tau_hot_start/1e-3} ms")
    overshoot_value = T_max - T_final
    overshoot = overshoot_value > 15
    max_ind = timeSeries['T'].values.argmax()
    if overshoot:
        # print(f"Overshoot by: {T_max - out_df['T']:.2f} K")
        # find temperature after peak where T = T_max - dT
        remaining_df = timeSeries.iloc[max_ind:]
        if overshoot_value > dT:
            remaining_df = remaining_df[(remaining_df['T'] <= T_max - dT)] #  choose first value where T drops by dT from T_Max
            tau_hot_end = remaining_df['age'].values[0] # in SECONDS
        elif overshoot_value > 0.5*dT:
            remaining_df = remaining_df[(remaining_df['T'] <= T_max - 0.5*dT)] # In this case, T_max - T_final > 0.5*dT, T_max - 0.5*dT > T_final
            tau_hot_end = remaining_df['age'].values[0] # SECONDS
        else:
            tau_hot_end = out_df['tau_sec_required'].values[0]*1e-3 # *1e-3 to convert to SECONDS
    else: 
        # print(f"No overshoot; using tau_CO = {out_df['tau_sec_required'].iloc[0]:.2} ms")
        tau_hot_end = out_df['tau_sec_required'].values[0]*1e-3
    tau_hot = (tau_hot_end - tau_hot_start)/1e-3 # convert to MILLISECONDS
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
    dt = timeSeries['age'].diff()
    dM = M.diff();
    dM_dt = dM/dt
    X_NO = timeSeries['X_NO']
    conc_NO = M*X_NO; 
    d_XNO = X_NO.diff()
    d_XNO_dt = d_XNO/dt
    dNO_dt = (X_NO*dM_dt) + (M*d_XNO_dt)
    
    tau_NOx_column = M/dNO_dt
    tau_NOx_NO_column = conc_NO/dNO_dt
    
    start_ind = int(ign_state.index.values[0])
    end_ind = int(end_state.index.values[0] + 1)
    tau_NOx = np.mean(tau_NOx_column.iloc[start_ind:end_ind+1])
    
    tau_NOx_NO = np.mean(tau_NOx_NO_column.iloc[start_ind:end_ind+1])

    start_time = ign_state['age'].values[0]/1e-3
    end_time = end_state['age'].values[0]/1e-3

    # tau_NOx_2 = np.mean()

    return tau_NOx/1e-3, tau_NOx_NO/1e-3, start_time, end_time

def get_Da(timeSeries, out_df, P=25*101325):
    eps = 0.0001*1e-3
    R_universal = 8.3144598; # J/mol-K or m3-Pa/mol-K
    M = (P/R_universal)/timeSeries['T'] # units of moles/volume
    dt = timeSeries['age'].diff()
    dM = M.diff();
    dM_dt = dM/dt
    X_NO = timeSeries['X_NO']
    conc_NO = M*X_NO; 
    d_XNO = X_NO.diff()
    d_XNO_dt = d_XNO/dt
    dNO_dt = (X_NO*dM_dt) + (M*d_XNO_dt)

    tau_NOx_column = M/dNO_dt
    tau_NOx_NO_column = conc_NO/dNO_dt    
    
    # ign_state = timeSeries[(timeSeries['age'] >= tau_hot_start - eps) & (timeSeries['age'] <= tau_hot_start + eps)] # system state at ignition
    # end_state = timeSeries[(timeSeries['age'] >= tau_hot_end - eps) & (timeSeries['age'] <= tau_hot_end + eps)] # system "end" state
    tau_hot_start = out_df['tau_ign_OH'].values[0]
    start_ind = int(timeSeries[(timeSeries['age'] >= tau_hot_start - eps) & (timeSeries['age'] <= tau_hot_start + eps)].index.values[0])
    # end_ind = int(end_state.index.values[0]) if int(end_state.index.values[0]) <= start_ind else start_ind
    dNO_dt_post_ign = dNO_dt.iloc[start_ind:]
    max_dNO_dt = max(dNO_dt_post_ign)
    dNO_dt_post_max = dNO_dt_post_ign.iloc[dNO_dt_post_ign.values.argmax():] # need to limit ourselves to post_max because don't want to get points before peak NO production
    perc_max = 0.5
    result_list = dNO_dt_post_max[dNO_dt_post_max <= perc_max*max_dNO_dt]
    # pdb.set_trace()
    # while len(result_list) < 1 and perc_max <= 0.5:
    #     perc_max += .1
    #     result_list = dNO_dt_post_max[dNO_dt_post_max <= perc_max*max_dNO_dt]

    end_ind = result_list.index.values[0]
    tau_hot_end = timeSeries['age'].iloc[end_ind]
    assert tau_hot_end > tau_hot_start, "tauHotEnd <= tauHotStart"
    # print(f"end_ind = {end_ind}, end time = {tau_hot_end/1e-3:.3f}, tau_sec_req = {out_df['tau_sec_required'].values[0]:.3f}")
    # print(f"max NO rate = {max_dNO_dt:.2f} = {dNO_dt_post_ign.iloc[dNO_dt_post_ign.values.argmax()]:.2f}")
    # print(f"len(result_list) = {len(result_list)}; end_ind = {end_ind}")

    tau_hot = (tau_hot_end - tau_hot_start)/1e-3 # convert to MILLISECONDS

    tau_NOx_NO = np.mean(tau_NOx_NO_column.iloc[start_ind:end_ind+1])/1e-3 # in milliseconds; note: have to actually do a weighted time-average if our timesteps are not equal (e.g., in the finite-mixing cases)

    Da = tau_hot/tau_NOx_NO
    # print(f"{tau_hot:.2f}, {tau_NOx_NO:.2f}, {Da:.2f}, {tau_hot_start/1e-3:.2f}, {tau_hot_end/1e-3:.2f}")
    return tau_hot, tau_NOx_NO, Da, tau_hot_start/1e-3, tau_hot_end/1e-3

def get_concNO(timeSeries, P):
    eps = 0.0001*1e-3
    R_universal = 8.3144598; # J/mol-K or m3-Pa/mol-K
    M = (P/R_universal)/timeSeries['T'] # units of moles/volume
    dt = timeSeries['age'].diff()
    dM = M.diff();
    dM_dt = dM/dt
    X_NO = timeSeries['X_NO']
    conc_NO = M*X_NO; 
    d_XNO = X_NO.diff()
    d_XNO_dt = d_XNO/dt
    dNO_dt = (X_NO*dM_dt) + (M*d_XNO_dt) 
    timeSeries['M'] = M
    timeSeries['conc_NO'] = conc_NO
    timeSeries['dNO_dt'] = dNO_dt
    timeSeries['dt'] = dt
    return timeSeries   

def get_Da_TempDrop(timeSeries, out_df, P=25*101325):
    timeSeries = get_concNO(timeSeries, P)

    tau_NOx_column = M/timeSeries['dNO_dt']
    tau_NOx_NO_column = timeSeries['conc_NO']/timeSeries['dNO_dt']
    
    # ign_state = timeSeries[(timeSeries['age'] >= tau_hot_start - eps) & (timeSeries['age'] <= tau_hot_start + eps)] # system state at ignition
    # end_state = timeSeries[(timeSeries['age'] >= tau_hot_end - eps) & (timeSeries['age'] <= tau_hot_end + eps)] # system "end" state
    tau_hot_start = out_df['tau_ign_OH'].values[0]
    start_ind = int(timeSeries[(timeSeries['age'] >= tau_hot_start - eps) & (timeSeries['age'] <= tau_hot_start + eps)].index.values[0])
    # end_ind = int(end_state.index.values[0]) if int(end_state.index.values[0]) <= start_ind else start_ind
    dNO_dt_post_ign = dNO_dt.iloc[start_ind:]
    max_dNO_dt = max(dNO_dt_post_ign)
    dNO_dt_post_max = dNO_dt_post_ign.iloc[dNO_dt_post_ign.values.argmax():] # need to limit ourselves to post_max because don't want to get points before peak NO production
    perc_max = 0.5
    result_list = dNO_dt_post_max[dNO_dt_post_max <= perc_max*max_dNO_dt]
    # pdb.set_trace()
    # while len(result_list) < 1 and perc_max <= 0.5:
    #     perc_max += .1
    #     result_list = dNO_dt_post_max[dNO_dt_post_max <= perc_max*max_dNO_dt]

    end_ind = result_list.index.values[0]
    tau_hot_end = timeSeries['age'].iloc[end_ind]
    assert tau_hot_end > tau_hot_start, "tauHotEnd <= tauHotStart"
    # print(f"end_ind = {end_ind}, end time = {tau_hot_end/1e-3:.3f}, tau_sec_req = {out_df['tau_sec_required'].values[0]:.3f}")
    # print(f"max NO rate = {max_dNO_dt:.2f} = {dNO_dt_post_ign.iloc[dNO_dt_post_ign.values.argmax()]:.2f}")
    # print(f"len(result_list) = {len(result_list)}; end_ind = {end_ind}")

    tau_hot = (tau_hot_end - tau_hot_start)/1e-3 # convert to MILLISECONDS

    tau_NOx_NO = np.mean(tau_NOx_NO_column.iloc[start_ind:end_ind+1])/1e-3 # in milliseconds; note: have to actually do a weighted time-average if our timesteps are not equal (e.g., in the finite-mixing cases)

    Da = tau_hot/tau_NOx_NO
    # print(f"{tau_hot:.2f}, {tau_NOx_NO:.2f}, {Da:.2f}, {tau_hot_start/1e-3:.2f}, {tau_hot_end/1e-3:.2f}")
    return tau_hot, tau_NOx_NO, Da, tau_hot_start/1e-3, tau_hot_end/1e-3

def get_avg_dNO(timeSeries, out_df, P = 25*101325): 
    timeSeries = get_concNO(timeSeries, P)    
    if isinstance(out_df, pd.DataFrame) and len(out_df) > 1:
        out_df = out_df.iloc[0]
    if isinstance(out_df, pd.Series):
        out_df = out_df.to_frame()
    if out_df.shape[1] == 1:
        out_df = out_df.T
    timeSeries_postIgn = timeSeries[(timeSeries['age'] >= out_df['tau_ign_OH'].values[0]) & (timeSeries['age'] <= out_df['tau_sec_required'].values[0]*1e-3)]    
    dNO = sum(timeSeries_postIgn['dNO_dt']*timeSeries_postIgn['dt'])
    avg_dNO = np.mean(dNO)
    return dNO, avg_dNO


def get_tempArea(timeSeries, out_df):
    if isinstance(out_df, pd.DataFrame) and len(out_df) > 1:
        out_df = out_df.iloc[0]
    if isinstance(out_df, pd.Series):
        out_df = out_df.to_frame()
    if out_df.shape[1] == 1:
        out_df = out_df.T
    
    timeSeries['dt'] = timeSeries['age'].diff();
    T_max = max(timeSeries['T'])
    timeSeries_postIgn = timeSeries[(timeSeries['age'] >= out_df['tau_ign_OH'].values[0]) & (timeSeries['age'] <= out_df['tau_sec_required'].values[0]*1e-3)]
    T_final = timeSeries['T'].iloc[-1] # assume final temperature = our desired temperature
    
    timeSeries_postIgn['T_frac'] = np.exp(-(T_max-timeSeries_postIgn['T'])/T_final)

    weighted_time = sum(timeSeries_postIgn['T_frac']*timeSeries_postIgn['dt'])
    assert weighted_time >= 0, "time must be positive... Something is wrong..."
    return weighted_time/1e-3
    # if T_max - T_final >= 30: # means overshoot

@jit(nopython=True, fastmath=True, cache=True)
def get_sys_NOCO(main_X, main_MW, sec_mass, sec_X, sec_MW):
    """
    Function to combine main and secondary stage NO and CO to obtain corrected NO and CO. 
    This function exists separate from the modify_timeTrace_mappable function to enable just-in-time compilation and speed up computations.
    
    Parameters:
    -----------
    main_X: nx4 np.ndarray with columns being ['X_NO', 'X_O2', 'X_H2O', 'X_CO'] 
    main_MW: nx1 np.ndarray containing average molecular weight of mixture in main burner 
    sec_mass: nx1 np.ndarray containing time history of secondary stage mass over time 
    sec_X: nx4 np.ndarray with columns being ['X_NO', 'X_O2', 'X_H2O', 'X_CO']
    sec_MW: nx1 np.ndarray containing time history of average molecular weight of secondary mixture
    """
    main_mass = np.max(sec_mass) - sec_mass # np.max(sec_mass) gives total mass of system, total - sec = main
    sec_moles = (sec_mass/sec_MW)
    sec_moles = sec_moles.reshape((len(sec_moles), 1))
    main_moles = (main_mass/main_MW)
    main_moles = main_moles.reshape((len(main_moles), 1))
    X_average = (main_moles * main_X + sec_moles * sec_X)/(sec_moles + main_moles)
    one_over_X_H2O = 1/X_average[:,2] 
    X_O2 = X_average[:,1]
    dry_O2_perc = X_O2 * one_over_X_H2O * 100
    factor = (20.9 - 15)/(20.9 - dry_O2_perc)
    dry_NO = X_average[:,0] * one_over_X_H2O
    dry_CO = X_average[:,3] * one_over_X_H2O
    corrected_NO = dry_NO*factor * 1e6 # convert to ppm
    corrected_CO = dry_CO*factor * 1e6
    return (corrected_NO, corrected_CO, sec_moles[:,0])

def modify_timeTrace_mappable(vit_df, trace_file):
    """
    Function to modify secondary stage time traces to include sys NO and CO. 
    WARNING: This WILL OVERRIDE the provided trace file. 
    
    Parameters:
    -----------
    vit_df: pandas DataFrame obtained from running Particle.get_timeHistory(dataFrame=True)
    trace_file: string containing full path to sec stage file. We're assuming that the file is located in ./SecStage/<trace_file.pickle>
    """
    trace_file = trace_file.split("/")[-1]
    sec_df = pd.read_parquet("SecStage/"+trace_file)
    if 'sys_NO_ppmvd' in sec_df.columns:
        return 0
#     Y_indices_vit = [i for i in range(0,len(vit_df.columns)) if 'Y_' in vit_df.columns[i]]
#     Y_indices_trace = [i for i in range(0, len(ww.columns)) if 'Y_' in ww.columns[i]]
#     Y_cols = sec_df.columns[Y_indices_vit].values
    X_cols = ['X_NO', 'X_O2', 'X_H2O', 'X_CO']
    main_X = vit_df.loc[(vit_df['age'] >= 0) & (vit_df['age'] <= sec_df['age'].iloc[-1]), X_cols].values
    main_MW = vit_df.loc[(vit_df['age'] >= 0) & (vit_df['age'] <= sec_df['age'].iloc[-1]), 'MW'].values
    sec_X = sec_df.loc[:,X_cols].values
    sec_mass = sec_df['mass'].values
    sec_MW = sec_df['MW'].values
    out = get_sys_NOCO(main_X, main_MW, sec_mass, sec_X, sec_MW)
    sec_df['sys_NO_ppmvd'] = out[0]
    sec_df['sys_CO_ppmvd'] = out[1]
    sec_df['n'] = out[2]
    sec_df.to_parquet("SecStage/"+trace_file)
    return trace_file

def modify_timeTrace(out_file, vit_df, num_processes=10):
    out_df = pd.read_parquet(out_file)   
#     for j, trace_file in tqdm(enumerate(out_df['reactor_file'].values)):
    pool = mp.Pool(num_processes)
    a = list(pool.map(partial(modify_timeTrace_mappable, vit_df), out_df['reactor_file'].values))
    return a
