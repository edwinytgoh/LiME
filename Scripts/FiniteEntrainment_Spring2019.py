""" 
Python script to run a finite-entrainment staged combustor. 
Designed to be run on an HPC cluster. 
"""
#* IMPORTS
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import multiprocessing as mp
import numpy as np 
import cantera as ct
import pandas as pd
from Utils.EquilTools import equil_CO
from Utils.finiteEntrainment import finite_entrainment_mappable, finite_entrainment_mappable, get_cons_idx
import Utils.CanteraTools as ctools
#* CONSTANTS
P = 25*ct.one_atm
T_fuel = 300 
T_air = 650
mech = 'gri30.xml'
fs_CH4 = 0.058387057492574147288255659304923028685152530670166015625
phi_global = 0.635
phi_main = 0.3719
CO_constraint = 1.25*equil_CO(phi_main)
phi_jet = np.inf
if __name__ == "__main__":
    #* STUDY PARAMS
    # tau_ent_sec_list = np.linspace(0.1, 4.0, 500)
    tau_ent_main_list = [float(sys.argv[1])]
    tau_ent_sec_list = np.linspace(tau_ent_main_list[0]/10, 0.1, num=250, endpoint=False)

    #* GENERATE TEST MATRIX
    argsList = []
    for i, tau_ent_sec in enumerate(tau_ent_sec_list):
        for j, tau_ent_main in enumerate(tau_ent_main_list):
            tau_sec = max(tau_ent_sec,tau_ent_main) + 2.0 # milliseconds
            tau_global = 15 + tau_sec
            argsList.append((phi_global, phi_main, tau_sec, tau_ent_main, tau_ent_sec, phi_jet, None, None, tau_global, T_fuel, T_air, P))

    #* RUN CASE
    # print(f"Running on {mp.cpu_count()} cores.")
    pool = mp.Pool(mp.cpu_count())
    out_list = list(pool.map(finite_entrainment_mappable, argsList))
    pool.close()
    out_df = pd.concat(out_list)
    out_df.to_csv(f"OutFiles/{sys.argv[2]}.csv")
    out_df.to_parquet(f"OutFiles/{sys.argv[2]}.pickle")