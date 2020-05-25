import os
import os.path

import cantera as ct
import numpy as np
import pandas as pd

from lime.utils import fs

pd.options.mode.chained_assignment = None  # default='warn'


def air(T: float = 300, P=25 * ct.one_atm, mech='gri30.xml'):
    air = ct.Solution(mech)
    air.TPX = T, P, {'N2': 0.79, 'O2': 0.21}
    return air


def CH4(T: float = 300, P=25 * ct.one_atm, mech='gri30.xml'):
    CH4 = ct.Solution(mech)
    CH4.TPX = T, P, {'CH4': 1.0}
    return CH4


def run_flame(gas, slope=0.01, curve=0.01):
    # Simulation parameters
    width = 0.06  # m
    # Flame object
    f = ct.FreeFlame(gas, width=width)
    f.set_refine_criteria(ratio=2, slope=slope, curve=curve, prune=min(slope, curve) / 1e3)
    f.transport_model = 'Multi'
    f.soret_enabled = True
    f.max_grid_points = 10000
    f.solve(loglevel=1, auto=True, refine_grid=True)
    # Convert distance into time:
    CH2O: int = gas.species_index('CH2O')
    X_CH2O = f.X[CH2O]
    maxIndex = np.arange(0, len(X_CH2O))[X_CH2O == max(X_CH2O)][0]
    maxIndex2 = X_CH2O.argmax()
    assert maxIndex == maxIndex2
    #     startingIndex = np.arange(0, len(X_CH2O))[X_CH2O >= X_CH2O[0] + 5][0]
    startingIndex = maxIndex
    #     startingIndex = np.arange(0, len(f.heat_release_rate))[f.heat_release_rate == max(f.heat_release_rate)][0]
    u_avg = np.array(f.u[startingIndex:] + f.u[startingIndex - 1:-1]) * 0.5
    dx = np.array(f.grid[startingIndex:] - f.grid[startingIndex - 1:-1])
    dt = dx / u_avg

    pre_time = [-999] * (startingIndex - 1)  # need to add some padding to the time array to account for dist = 0.0
    pre_time.extend([0])
    time = np.hstack(
        (np.array(pre_time),
         np.cumsum(dt)))  # numerically integrate (read: sum) dt to get array of times for each x location

    return f, time


def get_state_at_time(flame, tList, tRequired, mech='gri30.xml'):
    '''A function that gets the state at a desired point in time of a flame simulation performed using run_flame.
        Takes in a flame object, its associated time series, and the desired point in time (in seconds).
        Returns a new Cantera gas object with the desired state, and the corresponding index in the flame at which the point was found. 

        Example usage: gas, t_ind = get_state_at_time(flame, time, t_req)'''
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
            newGas.TPX = flame['T'].iloc[minIndex - 1], flame['P'].iloc[minIndex - 1], flame[moleFracs].iloc[
                minIndex - 1]
            assert ((newGas['NO'].X - flame['X_NO'].iloc[minIndex - 1] <= 1e-6))
            tau_vit = tRequired - tList[minIndex - 1]
            tau_vit_start = tList[minIndex - 1]
            assert (
                    tau_vit > 0)  # if the closest flameTime is larger than tRequired, the second closest should be less than, otherwise /that/ would be the closest...?
        elif ((tRequired - flameTime_closestTreq) > 1e-6):  # if closest time is more than 1e-6 s less than tReq
            newGas.TPX = flame['T'].iloc[minIndex], flame['P'].iloc[minIndex], flame[moleFracs].iloc[minIndex]
            assert ((newGas['NO'].X - flame['X_NO'].iloc[minIndex] <= 1e-6))
            tau_vit = tRequired - tList[minIndex]
            tau_vit_start = tList[minIndex]
        else:
            newGas.TPX = flame['T'].iloc[minIndex], flame['P'].iloc[minIndex], flame[moleFracs].iloc[minIndex]
            tau_vit = 0
            assert ((newGas['NO'].X - flame['X_NO'].iloc[minIndex] <= 1e-6))
    if tau_vit > 0:
        vitiator = ct.ConstPressureReactor(newGas)
        vitRN = ct.ReactorNet([vitiator])
        dt = 0.0001 * 1e-3
        vit_tList = np.arange(0, tau_vit, dt)
        vitArray = np.array([None] * len(vit_tList) * len(columnNames)).reshape(len(vit_tList), len(columnNames))
        # performance note: we don't need to do this. we can simply just advance to the desired tau_main
        for i in range(0, len(vit_tList)):
            MWi = vitiator.thermo.mean_molecular_weight
            vitArray[i, :] = np.hstack(
                [vel_final * dt, vel_final, vitiator.thermo.T, 0, MWi, vitiator.thermo.Y, vitiator.thermo.X,
                 vitiator.thermo.P])
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
    mainBurnerDF['NOppmvd'] = correct_nox(mainBurnerDF['X_NO'], mainBurnerDF['X_H2O'], mainBurnerDF['X_O2'])
    mainBurnerDF['COppmvd'] = correct_nox(mainBurnerDF['X_CO'], mainBurnerDF['X_H2O'], mainBurnerDF['X_O2'])
    return newGas, minIndex, mainBurnerDF


def fstoich(fuel={'CH4': 1}, ox={'O2': 0.21, 'N2': 0.79}, mech='gri30.xml'):
    '''Function that returns the stoichiometric fuel-to-air ratio (by mass) for a given fuel and oxidizer. 
    
    Example usage: f_stoich_ch4_air = fstoich() 
                    f_stoich_c2h4_o2 = fstoich(fuel={'C2H4':1}, ox={'O2':1})'''

    gas = ct.Solution(mech)
    gas.set_equivalence_ratio(1.0, fuel, ox)
    return np.sum(gas[fuel.keys()].Y) / np.sum(gas[ox.keys()].Y)


def calculate_flowrates(phi_global, phi_main, phi_sec):
    tol = -1e-14
    mfm = 1
    mam = 1 / (phi_main * fs)
    mas = (phi_global / phi_main - 1) / ((phi_sec - phi_global) * fs)
    mfs = phi_global / phi_main * (1 + (phi_global - phi_main) / (phi_sec - phi_global)) - 1

    # assert mfm >= tol and mfs >= tol and mam >= tol and mas >= tol, f"You have negative mass flow rates! mfm = {mfm}, mam = {mam}, mfs = {mfs}, mas = {mas}"

    m_perc = 100 / (mfm + mam + mas + mfs)
    # return {'mfm':mfm*m_perc, 'mam':mam*m_perc, 'mfs':mfs*m_perc, 'mas':mas*m_perc}
    return mfm * m_perc, mam * m_perc, mfs * m_perc, mas * m_perc


def solve_mass_airsplit(phi_global, phi_main, mdotTotal=1000, airSplit=1):
    # fs = 
    mfm = airSplit * fs * mdotTotal * (1 + fs * phi_global) ** (-1) * phi_main
    mam = airSplit * mdotTotal * (1 + fs * phi_global) ** (-1)
    mfs = (-1) * (1 + fs * phi_global) ** (-1) * (
            (-1) * fs * mdotTotal * phi_global + airSplit * fs * mdotTotal * phi_main)
    mas = (-1) * ((-1) + airSplit) * mdotTotal * (1 + fs * phi_global) ** (-1)
    return mfm, mam, mfs, mas


def solve_mass_phi_jet(phi_global, phi_main, phiJet, mdotTotal=1000):
    mam = 1
    mfm = mam * fs * phi_main
    mfs1 = mam * fs * (phi_global - phi_main)
    mas = mfs1 / (fs * (phiJet - phi_global))
    mfs = mfs1 + mas * fs * phi_global
    # Normalize
    mtot = mam + mfm + mas + mfs
    mam *= mdotTotal / mtot
    mfm *= mdotTotal / mtot
    mas *= mdotTotal / mtot
    mfs *= mdotTotal / mtot
    return mfm, mam, mfs, mas


def correct_nox(X_i, X_H2O, X_O2):
    dry_i = X_i / (1 - X_H2O)
    dry_O2 = X_O2 / (1 - X_H2O)
    corrected = dry_i * (20.9 - 15) / (20.9 - dry_O2 * 100)
    corrected_ppm = corrected * 1e6
    return corrected_ppm


def mix(streams, mdots, mech="gri30.xml", P=25 * 101325):
    # Create mixer gas: 
    mixerGas = ct.Solution(mech)
    mixerGas.TPX = [300, P, 'H2:1']

    # Create reactor with CHEMISTRY DISABLED: 
    mixer = ct.ConstPressureReactor(mixerGas)
    mixer.chemistry_enabled = False  # disable chemistry

    # For each stream (and given mass flow rate), connect mass flow controller to mixer: 
    mfcs = [ct.MassFlowController(ct.Reservoir(streams[i]), mixer, mdot=mdots[i]) for i in range(0, len(streams))]

    exhaust = ct.Reservoir(mixerGas)
    exhaust_mfc = ct.MassFlowController(mixer, exhaust, mdot=sum(mdots))

    rn = ct.ReactorNet([mixer])
    rn.advance_to_steady_state()

    return mixer.thermo


def premix(phi=0.4, fuel={'CH4': 1}, ox={'N2': 0.79, 'O2': 0.21}, mech='gri30.xml', P=25 * 101325, T_fuel=300, T_ox=650,
           M_total=1):
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


def run_main_burner(phi_main, tau_main, T_fuel=300, T_ox=650, P=25 * 101325, mech="gri30.xml", slope=0.01, curve=0.01,
                    filename=None):
    flameGas = premix(phi_main, P=P, mech=mech, T_fuel=T_fuel, T_ox=T_ox)
    # filename = '{0}_{1}-{2}_{3}-{4}_{5}'.format('phi_main', phi_main, 'P', P, )
    if filename == None:
        filename = '{0}_{1:.4f}.pickle'.format('phi_main', phi_main)
    if os.path.isfile(filename):
        mainBurnerDF = pd.read_parquet(filename)
        # mainBurnerDF = table.to_pandas()
        flameTime = mainBurnerDF.index.values
    else:
        flame, flameTime = run_flame(flameGas, slope=slope, curve=curve)
        columnNames = ['x', 'u', 'T', 'n', 'MW'] + ["Y_" + sn for sn in flameGas.species_names] + ["X_" + sn for sn in
                                                                                                   flameGas.species_names]
        flameData = np.concatenate(
            [np.array([flame.grid]), np.array([flame.u]), np.array([flame.T]), np.array([[0] * len(flame.T)]),
             np.array([[0] * len(flame.T)]), flame.Y, flame.X], axis=0)
        mainBurnerDF = pd.DataFrame(data=flameData.transpose(), index=flameTime, columns=columnNames)
        mainBurnerDF.index.name = 'Time'
        mainBurnerDF['P'] = flame.P
        # table = pa.Table.from_pandas(mainBurnerDF)
        # pq.write_table(table, filename)
        mainBurnerDF.to_parquet(filename, compression='gzip')

    vitiatedProd, flameCutoffIndex, mainBurnerDF = get_state_at_time(mainBurnerDF, flameTime, tau_main)
    vitReactor = ct.ConstPressureReactor(vitiatedProd, name='MainBurner')
    return vitReactor, mainBurnerDF


def get_global_equivalence_ratio(species_names, gas_template, oxidizers=[], ignore=[]):
    if not oxidizers:  # Default behavior, find all possible oxidizers
        oxidizers = [sn for sn in species_names]  # if
        # all(y not in s.composition for y in ['C', 'H', 'S'])]
        alpha = 0
        mol_O = 0
        for k, s in enumerate(gas_template.species):
            if s.name in ignore:
                continue
            else:  # elif s.name in oxidizers:
                mol_O += s.composition.get('O', 0) * gas_template.X[k]
                # else:
                nC = s.composition.get('C', 0)
                nH = s.composition.get('H', 0)
                # nO = s.composition.get('O', 0)
                nS = s.composition.get('S', 0)

                alpha += (2 * nC + nH / 2 + 2 * nS) * gas_template.X[k]

        if mol_O == 0:
            return float('inf')
        else:
            return alpha / mol_O