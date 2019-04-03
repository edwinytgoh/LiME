import Utils.CanteraTools as ctools 
import cantera as ct 
import numpy as np 
import os
import pandas as pd
fs_CH4 = 0.058387057492574147288255659304923028685152530670166015625 

equil_file = os.path.join(os.path.dirname(__file__), 'equil_results.pickle')

def equil_CO(phi, T_air = 650, T_fuel = 300, P = 25*ct.one_atm, mech="gri30.xml"):
    if os.path.isfile(equil_file):
        df = pd.read_parquet(equil_file)
        return np.interp(phi, df['phi'].values, df['CO_ppmvd'].values)
    else:
        return equil(phi, T_air, T_fuel, P, mech)[1]

def equil(phi, T_air = 650, T_fuel = 300, P = 25*ct.one_atm, mech="gri30.xml"): 
    if os.path.isfile(equil_file):
        df = pd.read_parquet(equil_file)
        T = np.interp(phi, df['phi'].values, df['T_eq'].values)
        CO = np.interp(phi, df['phi'].values, df['CO_ppmvd'].values)
        NO = np.interp(phi, df['phi'].values, df['NO_ppmvd'].values)
        return np.hstack([T, CO, NO])
    
    gas = ct.Solution(mech);  
    m_air = 1
    m_fuel = phi*fs_CH4
    mixture =  ctools.mix([ctools.air(T_air, P, mech), ctools.CH4(T_fuel, P, mech)], [m_air, m_fuel], mech, P)
    mixture.equilibrate('HP');  
    CO_ppmvd = ctools.correctNOx(mixture['CO'].X, mixture['H2O'].X, mixture['O2'].X) 
    NO_ppmvd = ctools.correctNOx(mixture['NO'].X, mixture['H2O'].X, mixture['O2'].X) 
    return np.hstack([mixture.T, CO_ppmvd, NO_ppmvd]) 


if __name__ == "__main__":
    import multiprocessing as mp 
    import pandas as pd
    print(f"Running on {mp.cpu_count()} cores.")
    pool = mp.Pool(mp.cpu_count())
    phi_list = np.hstack([np.arange(0.4,1.1,0.001), np.arange(1.1001, 5, 0.1)])
    phi_list = phi_list.reshape([len(phi_list), 1])
    results_array = np.vstack(list(pool.map(equil, phi_list)))
    pool.close()
    results = pd.DataFrame(data=np.hstack([phi_list,results_array]), columns=['phi', 'T_eq', 'CO_ppmvd', 'NO_ppmvd'])
    results.to_parquet(equil_file) # to_parquet is faster than to_csv