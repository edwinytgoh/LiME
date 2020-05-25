import multiprocessing as mp
import os
import shutil
from functools import partial

import numpy as np
import pandas as pd
from numba import jit


def get_tauHot(timeSeries, out_df, dT=200):  ## REMEMBER TO CHECK UNITS!!!
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
    tau_hot_start = out_df['tau_ign_OH'].values[0]  # seconds
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
            remaining_df = remaining_df[
                (remaining_df['T'] <= T_max - dT)]  # choose first value where T drops by dT from T_Max
            tau_hot_end = remaining_df['age'].values[0]  # in SECONDS
        elif overshoot_value > 0.5 * dT:
            remaining_df = remaining_df[(remaining_df[
                                             'T'] <= T_max - 0.5 * dT)]  # In this case, T_max - T_final > 0.5*dT, T_max - 0.5*dT > T_final
            tau_hot_end = remaining_df['age'].values[0]  # SECONDS
        else:
            tau_hot_end = out_df['tau_sec_required'].values[0] * 1e-3  # *1e-3 to convert to SECONDS
    else:
        # print(f"No overshoot; using tau_CO = {out_df['tau_sec_required'].iloc[0]:.2} ms")
        tau_hot_end = out_df['tau_sec_required'].values[0] * 1e-3
    tau_hot = (tau_hot_end - tau_hot_start) / 1e-3  # convert to MILLISECONDS
    return tau_hot, tau_hot_start, tau_hot_end


def get_tauNOx(timeSeries, tauHot_start, tauHot_end, P=25 * 101325):
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
    eps = 0.0001 * 1e-3
    ign_state = timeSeries[(timeSeries['age'] >= tauHot_start - eps) & (
            timeSeries['age'] <= tauHot_start + eps)]  # system state at ignition
    end_state = timeSeries[
        (timeSeries['age'] >= tauHot_end - eps) & (timeSeries['age'] <= tauHot_end + eps)]  # system "end" state

    R_universal = 8.3144598  # J/mol-K or m3-Pa/mol-K
    M = (P / R_universal) / timeSeries['T']  # units of moles/volume
    dt = timeSeries['age'].diff()
    dM = M.diff();
    dM_dt = dM / dt
    X_NO = timeSeries['X_NO']
    conc_NO = M * X_NO;
    d_XNO = X_NO.diff()
    d_XNO_dt = d_XNO / dt
    dNO_dt = (X_NO * dM_dt) + (M * d_XNO_dt)

    tau_NOx_column = M / dNO_dt
    tau_NOx_NO_column = conc_NO / dNO_dt

    start_ind = int(ign_state.index.values[0])
    end_ind = int(end_state.index.values[0] + 1)
    tau_NOx = np.mean(tau_NOx_column.iloc[start_ind:end_ind + 1])

    tau_NOx_NO = np.mean(tau_NOx_NO_column.iloc[start_ind:end_ind + 1])

    start_time = ign_state['age'].values[0] / 1e-3
    end_time = end_state['age'].values[0] / 1e-3

    # tau_NOx_2 = np.mean()

    return tau_NOx / 1e-3, tau_NOx_NO / 1e-3, start_time, end_time


def get_Da(timeSeries, out_df, P=25 * 101325):
    eps = 0.0001 * 1e-3
    R_universal = 8.3144598 # J/mol-K or m3-Pa/mol-K
    M = (P / R_universal) / timeSeries['T']  # units of moles/volume
    dt = timeSeries['age'].diff()
    dM = M.diff()
    dM_dt = dM / dt
    X_NO = timeSeries['X_NO']
    conc_NO = M * X_NO;
    d_XNO = X_NO.diff()
    d_XNO_dt = d_XNO / dt
    dNO_dt = (X_NO * dM_dt) + (M * d_XNO_dt)

    tau_NOx_column = M / dNO_dt
    tau_NOx_NO_column = conc_NO / dNO_dt

    # ign_state = timeSeries[(timeSeries['age'] >= tau_hot_start - eps) & (timeSeries['age'] <= tau_hot_start + eps)] # system state at ignition
    # end_state = timeSeries[(timeSeries['age'] >= tau_hot_end - eps) & (timeSeries['age'] <= tau_hot_end + eps)] # system "end" state
    tau_hot_start = out_df['tau_ign_OH'].values[0]
    start_ind = int(timeSeries[(timeSeries['age'] >= tau_hot_start - eps) & (
            timeSeries['age'] <= tau_hot_start + eps)].index.values[0])
    # end_ind = int(end_state.index.values[0]) if int(end_state.index.values[0]) <= start_ind else start_ind
    dNO_dt_post_ign = dNO_dt.iloc[start_ind:]
    max_dNO_dt = max(dNO_dt_post_ign)
    dNO_dt_post_max = dNO_dt_post_ign.iloc[
                      dNO_dt_post_ign.values.argmax():]  # need to limit ourselves to post_max because don't want to get points before peak NO production
    perc_max = 0.5
    result_list = dNO_dt_post_max[dNO_dt_post_max <= perc_max * max_dNO_dt]
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

    tau_hot = (tau_hot_end - tau_hot_start) / 1e-3  # convert to MILLISECONDS

    tau_NOx_NO = np.mean(tau_NOx_NO_column.iloc[
                         start_ind:end_ind + 1]) / 1e-3  # in milliseconds; note: have to actually do a weighted time-average if our timesteps are not equal (e.g., in the finite-mixing cases)

    Da = tau_hot / tau_NOx_NO
    # print(f"{tau_hot:.2f}, {tau_NOx_NO:.2f}, {Da:.2f}, {tau_hot_start/1e-3:.2f}, {tau_hot_end/1e-3:.2f}")
    return tau_hot, tau_NOx_NO, Da, tau_hot_start / 1e-3, tau_hot_end / 1e-3


def get_concNO(timeSeries, P):
    eps = 0.0001 * 1e-3
    R_universal = 8.3144598;  # J/mol-K or m3-Pa/mol-K
    M = (P / R_universal) / timeSeries['T']  # units of moles/volume
    dt = timeSeries['age'].diff()
    dM = M.diff();
    dM_dt = dM / dt
    X_NO = timeSeries['X_NO']
    conc_NO = M * X_NO;
    d_XNO = X_NO.diff()
    d_XNO_dt = d_XNO / dt
    dNO_dt = (X_NO * dM_dt) + (M * d_XNO_dt)
    timeSeries['M'] = M
    timeSeries['conc_NO'] = conc_NO
    timeSeries['dNO_dt'] = dNO_dt
    timeSeries['dt'] = dt
    return timeSeries


def get_Da_TempDrop(timeSeries, out_df, P=25 * 101325):
    timeSeries = get_concNO(timeSeries, P)

    tau_NOx_column = M / timeSeries['dNO_dt']
    tau_NOx_NO_column = timeSeries['conc_NO'] / timeSeries['dNO_dt']

    # ign_state = timeSeries[(timeSeries['age'] >= tau_hot_start - eps) & (timeSeries['age'] <= tau_hot_start + eps)] # system state at ignition
    # end_state = timeSeries[(timeSeries['age'] >= tau_hot_end - eps) & (timeSeries['age'] <= tau_hot_end + eps)] # system "end" state
    tau_hot_start = out_df['tau_ign_OH'].values[0]
    start_ind = int(timeSeries[(timeSeries['age'] >= tau_hot_start - eps) & (
            timeSeries['age'] <= tau_hot_start + eps)].index.values[0])
    # end_ind = int(end_state.index.values[0]) if int(end_state.index.values[0]) <= start_ind else start_ind
    dNO_dt_post_ign = dNO_dt.iloc[start_ind:]
    max_dNO_dt = max(dNO_dt_post_ign)
    dNO_dt_post_max = dNO_dt_post_ign.iloc[
                      dNO_dt_post_ign.values.argmax():]  # need to limit ourselves to post_max because don't want to get points before peak NO production
    perc_max = 0.5
    result_list = dNO_dt_post_max[dNO_dt_post_max <= perc_max * max_dNO_dt]
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

    tau_hot = (tau_hot_end - tau_hot_start) / 1e-3  # convert to MILLISECONDS

    tau_NOx_NO = np.mean(tau_NOx_NO_column.iloc[
                         start_ind:end_ind + 1]) / 1e-3  # in milliseconds; note: have to actually do a weighted time-average if our timesteps are not equal (e.g., in the finite-mixing cases)

    Da = tau_hot / tau_NOx_NO
    # print(f"{tau_hot:.2f}, {tau_NOx_NO:.2f}, {Da:.2f}, {tau_hot_start/1e-3:.2f}, {tau_hot_end/1e-3:.2f}")
    return tau_hot, tau_NOx_NO, Da, tau_hot_start / 1e-3, tau_hot_end / 1e-3


def get_avg_dNO(timeSeries, out_df, P=25 * 101325):
    timeSeries = get_concNO(timeSeries, P)
    if isinstance(out_df, pd.DataFrame) and len(out_df) > 1:
        out_df = out_df.iloc[0]
    if isinstance(out_df, pd.Series):
        out_df = out_df.to_frame()
    if out_df.shape[1] == 1:
        out_df = out_df.T
    timeSeries_postIgn = timeSeries[(timeSeries['age'] >= out_df['tau_ign_OH'].values[0]) & (
            timeSeries['age'] <= out_df['tau_sec_required'].values[0] * 1e-3)]
    dNO = sum(timeSeries_postIgn['dNO_dt'] * timeSeries_postIgn['dt'])
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
    timeSeries_postIgn = timeSeries[(timeSeries['age'] >= out_df['tau_ign_OH'].values[0]) & (
            timeSeries['age'] <= out_df['tau_sec_required'].values[0] * 1e-3)]
    T_final = timeSeries['T'].iloc[-1]  # assume final temperature = our desired temperature

    timeSeries_postIgn['T_frac'] = np.exp(-(T_max - timeSeries_postIgn['T']) / T_final)

    weighted_time = sum(timeSeries_postIgn['T_frac'] * timeSeries_postIgn['dt'])
    assert weighted_time >= 0, "time must be positive... Something is wrong..."
    return weighted_time / 1e-3
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
    main_mass = np.max(sec_mass) - sec_mass  # np.max(sec_mass) gives total mass of system, total - sec = main
    sec_moles = (sec_mass / sec_MW)
    sec_moles = sec_moles.reshape((len(sec_moles), 1))

    main_moles = (main_mass / main_MW)
    main_moles = main_moles.reshape((len(main_moles), 1))
    X_average = (main_moles * main_X + sec_moles * sec_X) / (sec_moles + main_moles)
    one_over_X_H2O = 1 / (1 - X_average[:, 2])
    X_O2 = X_average[:, 1]
    dry_O2_perc = X_O2 * one_over_X_H2O * 100
    factor = (20.9 - 15) / (20.9 - dry_O2_perc)
    dry_NO = X_average[:, 0] * one_over_X_H2O
    dry_CO = X_average[:, 3] * one_over_X_H2O
    corrected_NO = dry_NO * factor * 1e6  # convert to ppm
    corrected_CO = dry_CO * factor * 1e6
    return (corrected_NO, corrected_CO, sec_moles[:, 0])


def modify_timeTrace_mappable(vit_df, trace_file, sec_stage_path="SecStage/"):
    """
    Function to modify secondary stage time traces to include sys NO and CO.
    WARNING: This WILL OVERRIDE the provided trace file.

    Parameters:
    -----------
    vit_df: pandas DataFrame obtained from running Particle.get_time_history(dataframe=True)
    trace_file: string containing full path to sec stage file. We're assuming that the file is located in ./SecStage/<trace_file.pickle>
    """
    trace_file = trace_file.split(os.pathsep)[-1]
    low_temp_path = sec_stage_path + "LowTemp/"
    os.makedirs(low_temp_path, exist_ok=True)
    try:
        sec_df = pd.read_parquet(sec_stage_path + trace_file)
    except Exception as e:
        return f"{e}: {trace_file}"
    if sec_df['T'].iloc[-1] < 1900:  # TODO: don't hardcode this
        shutil.move(sec_stage_path + trace_file, low_temp_path + trace_file)
        return f"Temp too low ({sec_df['T'].iloc[-1]:.2f} K): {trace_file}"
        # if 'sys_NO_ppmvd' in sec_df.columns:
    #     return f"Done: {trace_file}"
    #     Y_indices_vit = [i for i in range(0,len(vit_df.columns)) if 'Y_' in vit_df.columns[i]]
    #     Y_indices_trace = [i for i in range(0, len(ww.columns)) if 'Y_' in ww.columns[i]]
    #     Y_cols = sec_df.columns[Y_indices_vit].values
    X_cols = ['X_NO', 'X_O2', 'X_H2O', 'X_CO']

    # Extract main burner state until end of secondary reactor
    main_X = vit_df.loc[(vit_df['age'] >= 0) & (vit_df['age'] <= sec_df['age'].iloc[-1]), X_cols].values
    main_MW = vit_df.loc[(vit_df['age'] >= 0) & (vit_df['age'] <= sec_df['age'].iloc[-1]), 'MW'].values

    sec_X = sec_df.loc[:, X_cols].values
    sec_mass = sec_df['mass'].values
    sec_MW = sec_df['MW'].values

    if len(main_MW) < len(sec_mass):
        dt = vit_df['age'].diff().mean()
        missing_timesteps = np.arange(vit_df['age'].iloc[-1] + dt, sec_df['age'].iloc[-1], dt)
        new_MW = np.interp(missing_timesteps, vit_df['age'].iloc[-100:], vit_df['MW'].iloc[-100:])
        main_NW = np.append(main_MW, new_MW)
        new_X = np.vstack(
            [np.interp(missing_timesteps, vit_df['age'].iloc[-100:], vit_df['X_NO'].iloc[-100:]),
             np.interp(missing_timesteps, vit_df['age'].iloc[-100:], vit_df['X_O2'].iloc[-100:]),
             np.interp(missing_timesteps, vit_df['age'].iloc[-100:], vit_df['X_H2O'].iloc[-100:]),
             np.interp(missing_timesteps, vit_df['age'].iloc[-100:], vit_df['X_CO'].iloc[-100:]),
             ]).T
        main_X = np.append(main_X, new_X)

    out = get_sys_NOCO(main_X, main_MW, sec_mass, sec_X, sec_MW)
    sec_df['sys_NO_ppmvd'] = out[0]
    sec_df['sys_CO_ppmvd'] = out[1]
    sec_df['n'] = out[2]
    sec_df.to_parquet(sec_stage_path + trace_file)
    return trace_file


def modify_timeTrace(trace_files, vit_df, num_processes=10):
    # out_df = pd.read_parquet(out_file)
    #     for j, trace_file in tqdm(enumerate(out_df['reactor_file'].values)):
    pool = mp.Pool(num_processes)
    a = list(pool.map(partial(modify_timeTrace_mappable, vit_df), trace_files))
    return a


@jit(nopython=True)
def get_cons_idx(COHistory, CO_constraint) -> int:
    # ! if all(COHistory < CO_constraint), get error: index -1 is out of bounds for axis 0 with size 0
    # try:
    CO_greater_idx = np.arange(len(COHistory))[
        COHistory > CO_constraint]  # all indices where CO is greater than constraint
    if len(CO_greater_idx) >= 1:
        return max(min(CO_greater_idx[-1] + 1, len(COHistory) - 1), 0)
    else:
        return -1
    # except:
    #     return -1


@jit(nopython=True)
def get_weighted_temperature(time, temp, factor=1):
    weighted_temp = np.sum(np.exp(-factor / temp)[1:] * (time[1:] - time[0:-1]))
    return weighted_temp


@jit(nopython=True)
def central_difference(y, x):
    sorted_idx = np.argsort(x)
    x_sorted = x[sorted_idx]
    y_sorted = y[sorted_idx]
    dydx = np.zeros(y.shape)
    for i in range(1, len(y)):
        dydx[sorted_idx[i]] = (y_sorted[i + 1] - y_sorted[i - 1]) / (x_sorted[i + 1] - x_sorted[i - 1])
    #     for i in range(1,len(y)):
    #         dydx[i] = (y[i + 1] - y[i - 1])/(x[i + 1] - x[i - 1])
    dydx[sorted_idx[0]] = np.nan
    dydx[sorted_idx[-1]] = np.nan
    return dydx.astype(np.float64)


def post_process_reactor_trace(trace_file, sec_stage_path="SecStage/", CO_constraint=32, dataframe=False):
    trace_file = trace_file.split(os.pathsep)[-1]
    try:
        sec_df = pd.read_parquet(sec_stage_path + trace_file)
    except Exception as e:
        return f"{e}: {trace_file}"
    cons_idx = get_cons_idx(sec_df['CO_ppmvd'].values, CO_constraint)
    cons_row = sec_df.iloc[cons_idx] if cons_idx >= 0 else None
    columns = [
        'tau_ent_main (s)',
        'tau_ent_sec (s)',
        'tau_ent_ratio',
        'phi_global',
        'phi_main',
        'phi_jet',
        'phi_jet_norm',
        'sys_NO_ppmvd',
        'sys_CO_ppmvd',
        'NO_ppmvd_constraint',
        'CO_ppmvd_constraint',
        'T_constraint',
        'tau_sec_required (s)',
        'phi_constraint',
        'tau_ign_maxGradT (s)',
        'NO_ppmvd_ign',
        'CO_ppmvd_ign',
        'T_ign',
        'phi_ign',
        'X_CH4_ign',
        'X_CO2_ign',
        'X_O2_ign',
        'X_OH_ign',
        'phi_init',
        'T_init',
        'T_max',
        'T_median',
        'T_mean',
        'T_weighted',  # integral of exp(-38,379/T)*dt
        'T_weighted_eq',
        'T_weighted_1950',
        'tau_hot',  # time where T > 1950 K
    ]
    tau_ent_main = float(
        trace_file[trace_file.find('tauEntMain') + len('tauEntMain'):trace_file.find('tauEntSec') - 1]) * 1e-3
    tau_ent_sec = float(
        trace_file[trace_file.find('tauEntSec') + len('tauEntSec'):trace_file.find('phiJetNorm') - 1]) * 1e-3
    tau_ent_ratio = tau_ent_main / tau_ent_sec
    phi_global = float(trace_file[trace_file.find('phiGlobal') + len('phiGlobal'):trace_file.find('phiMain') - 1])
    phi_main = float(trace_file[trace_file.find('phiMain') + len('phiMain'):trace_file.find('tauEntMain') - 1])
    phi_jet_norm = float(trace_file[trace_file.find('phiJetNorm') + len('phiJetNorm'):trace_file.find('dt') - 1])
    phi_jet = phi_jet_norm / (1 - phi_jet_norm) if phi_jet_norm < 1 else np.inf
    if isinstance(cons_row, pd.Series):
        sys_NO_ppmvd = cons_row['sys_NO_ppmvd']
        sys_CO_ppmvd = cons_row['sys_CO_ppmvd']
        NO_ppmvd_constraint = cons_row['NO_ppmvd']
        CO_ppmvd_constraint = cons_row['CO_ppmvd']
        T_constraint = cons_row['T']
        tau_sec_required = cons_row['age']
        phi_constraint = cons_row['phi']
    elif cons_row == None:
        sys_NO_ppmvd = -1000
        sys_CO_ppmvd = -1000
        NO_ppmvd_constraint = -1000
        CO_ppmvd_constraint = -1000
        T_constraint = -1000
        tau_sec_required = -1000
        phi_constraint = -1000
    # define ignition to be max temp gradient
    ign_idx = np.nanargmax(central_difference(sec_df['T'].values, sec_df['age'].values))
    ign_row = sec_df.iloc[ign_idx]
    NO_ppmvd_ign = ign_row['NO_ppmvd']
    CO_ppmvd_ign = ign_row['CO_ppmvd']
    T_ign = ign_row['T']
    phi_ign = ign_row['phi']
    X_CH4_ign = ign_row['X_CH4']
    X_CO2_ign = ign_row['X_CO2']
    X_O2_ign = ign_row['X_O2']
    X_OH_ign = ign_row['X_OH']
    tau_ign_maxGradT = ign_row['age']

    phi_init = sec_df['phi'].iloc[0]
    T_init = sec_df['T'].iloc[0]
    T_max = sec_df['T'].max()
    T_median = sec_df['T'].median()
    T_mean = sec_df['T'].mean()
    T_weighted = get_weighted_temperature(sec_df['age'].values, sec_df['T'].values)
    T_weighted_eq = get_weighted_temperature(sec_df['age'].values, sec_df['T'].values, factor=sec_df['T'].iloc[-1])
    T_weighted_1950 = get_weighted_temperature(sec_df['age'].values, sec_df['T'].values, factor=1950)
    tau_hot = len(sec_df[sec_df['T'] >= 1950]) * sec_df['age'].diff().mean()

    outputs = [
        tau_ent_main,
        tau_ent_sec,
        tau_ent_ratio,
        phi_global,
        phi_main,
        phi_jet,
        phi_jet_norm,
        sys_NO_ppmvd,
        sys_CO_ppmvd,
        NO_ppmvd_constraint,
        CO_ppmvd_constraint,
        T_constraint,
        tau_sec_required,
        phi_constraint,
        tau_ign_maxGradT,
        NO_ppmvd_ign,
        CO_ppmvd_ign,
        T_ign,
        phi_ign,
        X_CH4_ign,
        X_CO2_ign,
        X_O2_ign,
        X_OH_ign,
        phi_init,
        T_init,
        T_max,
        T_median,
        T_mean,
        T_weighted,
        T_weighted_eq,
        T_weighted_1950,
        tau_hot,
    ]
    if dataframe:
        return pd.DataFrame(data=np.array(outputs)[np.newaxis], columns=columns)
    else:
        return outputs
