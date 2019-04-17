from .CanteraTools import *
import numpy as np 
import sys 
import os
import pdb
import cantera as ct
from Particle import Particle
from .EquilTools import equil
milliseconds = 1e-3

def finite_entrainment(phi_global, phi_main, tau_sec, tau_ent_main, tau_ent_sec,
                        phi_jet=np.inf, main_flowFunc=None, sec_flowFunc=None,
                        tau_global=10, T_fuel=300, T_ox=650, P = 25*ct.one_atm, dt:float=0.001*1e-3,mech="gri30.xml", CO_constraint=None, CO_constraint_active=False, write_df=True, out_dir=os.getcwd(), flame_file=None):
    """
    Example flame_file: flame_file=os.path.join(flame_library_location, f"phi_main_{phi_main:.4f}_GRI30_25atm.pickle")
    """
    tau_global *= milliseconds 
    tau_sec *= milliseconds 
    tau_ent_main *= milliseconds 
    tau_ent_sec *= milliseconds
    tau_main = tau_global - max(0,tau_sec) 
    phi_jet_norm = 1 if phi_jet == np.inf else phi_jet/(1+phi_jet)
    mfm, mam, mfs, mas = calculate_flowRates(phi_global, phi_main, phi_jet)
    T_eq, CO_eq, NO_eq = equil(phi_global)
    if CO_constraint == None: 
        CO_constraint = 1.25*CO_eq

    [main_reactor, main_burner_DF] = runMainBurner(phi_main, tau_main, filename=flame_file) #* Run main burner
    main_reservoir = ct.Reservoir(contents=main_reactor.thermo)
    main_particle = Particle.fromReactor(main_reactor, mech=mech, particle_mass = mfm + mam)
    filename = f"phiGlobal{phi_global:.4f}_phiMain{phi_main:.5f}_tauEntMain{tau_ent_main/milliseconds:.6f}_tauEntSec{tau_ent_sec/milliseconds:.6f}_phiJetNorm{phi_jet_norm:.6f}_dt{dt:.5e}_P{P:.4f}_mech-{mech}"
#* Create Secondary Reservoir:  
    if np.isinf(phi_jet):
        secondary_gas = CH4(T_fuel, P, mech)
    else:
        secondary_gas = mix([CH4(T_fuel, P, mech), air(T_ox, P, mech)], [fs*phi_jet, 1], P = P)
    
    sec_particle = Particle.fromGas(secondary_gas, particle_mass = mfs + mas)
    sec_reservoir = ct.Reservoir(contents=sec_particle)

#* Create secondary STAGE:  
    sec_stage = Particle(mech, P=P, state_vec=sec_particle.state, particle_mass=1e-8)
    sec_stage.net.add_reactor(main_reactor) #TODO: check to see if this line actually changes main_reservoir contents
    sec_stage.net.add_reactor(sec_particle.reactor)
    sec_stage_file = os.path.join(out_dir, 'SecStage', f"{filename}.pickle")
    sec_stage_rate_file = os.path.join(out_dir, 'SecStage_Rate', f"{filename}.pickle")
    if not os.path.isdir(os.path.join(out_dir, 'SecStage')): 
        os.mkdir(os.path.join(out_dir, 'SecStage'))
    if not os.path.isdir(os.path.join(out_dir, 'SecStage_Rate')):
        os.mkdir(os.path.join(out_dir, 'SecStage_Rate'))
#* Create main and secondary mass flow controllers and connect them to secondary STAGE 
    mfc_main = ct.MassFlowController(main_reservoir, sec_stage.reactor)
    if (main_flowFunc == None):
        mdot_main = (mfm + mam)/tau_ent_main
        main_flowFunc = lambda t: mdot_main if t <= tau_ent_main else 0    
    mfc_main.set_mass_flow_rate(main_flowFunc)

    mfc_sec = ct.MassFlowController(sec_reservoir, sec_stage.reactor)
    if (sec_flowFunc == None):        
        mdot_sec = (mfs + mas)/tau_ent_sec
        sec_flowFunc = lambda t: mdot_sec if t <= tau_ent_sec else 0    
    mfc_sec.set_mass_flow_rate(sec_flowFunc)

#* Run simulation: 
    files_exist = False
    if (os.path.isfile(sec_stage_file)) and (os.path.isfile(sec_stage_rate_file)): # don't run simulation if files already exist
        files_exist = True
        try:
            sec_stage_DF = pd.read_parquet(sec_stage_file)
            rate_DF = pd.read_parquet(sec_stage_rate_file)
            assert sec_stage_DF.shape[0] > 0, 'Length of sec_stage_DF = 0, deleting this file'
            assert (T_eq - sec_stage_DF['T'].iloc[-1]) <= 100, "file exists but hasn't ignited" #TODO: Removing the file might be a bit drastic, but I'm a little lazy to create a sec_stage Particle and resume running at the moment >w<
        except Exception: # Delete files if corrupt 
            os.remove(sec_stage_file)
            os.remove(sec_stage_rate_file)
            files_exist = False
    
    if not files_exist: # i.e. if files_exist == False
        main_mass_injected = 0; 
        sec_mass_injected = 0; 
        for i,t in enumerate(np.arange(0, tau_sec, dt)):
            sec_stage.react(dt)
            main_reservoir.syncState()
            main_mass_injected += mfc_main.mdot(t)*dt; 
            sec_mass_injected += mfc_sec.mdot(t)*dt            
            main_particle.mass = max(main_particle.mass - mfc_main.mdot(t)*dt, 0)
            main_particle.age += dt
            main_particle.timeHistory_list.append(main_particle.outState)
            sec_particle.mass = max(sec_particle.mass - mfc_sec.mdot(t)*dt, 0)
            sec_particle.age += dt
            sec_particle.timeHistory_list.append(sec_particle.outState)

            
        #! VERY ARBITRARY. NOT GOOD. 
        #TODO: See if there's a better way to do this:
        # while ((T_eq - sec_stage.T > 100) and (sec_stage.age <= 20*1e-3)): # Run extra steps if sec_stage.T is still less than T_eq by more than 50 K
        #     sec_stage.react(dt)
        #     main_mass_injected += mfc_main.mdot(t)*dt;
        #     sec_mass_injected += mfc_sec.mdot(t)*dt
        sec_stage_DF = sec_stage.get_timeHistory(dataFrame=True, deleteFirstElem=False)
        rate_DF = sec_stage.get_rateHistory(dataFrame=True, deleteFirstElem=False)
        jet_DF = sec_particle.get_timeHistory(dataFrame=True, deleteFirstElem=False)
        vit_df = main_particle.get_timeHistory(dataFrame=True, deleteFirstElem=False)
#* Post-process:

    #* Add system NO and CO: 
    #Y_cols = [f"X_{s}" for s in main_particle.species_names]
        Y_cols = ['X_NO', 'X_O2', 'X_H2O', 'X_CO']
        sec_Y = sec_stage_DF.loc[:,Y_cols].values
        sec_mass = sec_stage_DF['mass'].values
        sec_moles = (sec_mass/sec_stage_DF['MW'].values)[:, np.newaxis]
        
        # jet_mass = (np.zeros(len(sec_stage_DF['age'])) + (mfs + mas)) - np.array([sec_flowFunc(t)*dt for t in sec_stage_DF['age'].values])
        jet_Y = jet_DF.loc[:,Y_cols].values
        jet_mass = jet_DF['mass'].values
        jet_moles = (jet_mass/jet_DF['MW'].values)[:,np.newaxis]
    
    
        main_Y = vit_df.loc[(vit_df['age'] >= 0) & (vit_df['age'] <= sec_stage_DF['age'].iloc[-1]), Y_cols].values
        main_MW = vit_df.loc[(vit_df['age'] >= 0) & (vit_df['age'] <= sec_stage_DF['age'].iloc[-1]), 'MW'].values
        # main_mass = sec_stage_DF['mass'].max() - sec_mass
        main_mass = vit_df['mass'].values
        main_moles = (main_mass/main_MW)[:, np.newaxis]
    
        X_average = (main_moles * main_Y + sec_moles * sec_Y + jet_moles * jet_Y)/(sec_moles + main_moles + jet_moles)
        sec_stage_DF['sys_NO_ppmvd'] = correctNOx(X_average[:,0], X_average[:,2], X_average[:,1])
        sec_stage_DF['sys_CO_ppmvd'] = correctNOx(X_average[:,3], X_average[:,2], X_average[:,1])
        sec_stage_DF['n'] = sec_moles[:,0]
        sec_stage_DF['sys_mass'] = main_mass + jet_mass + sec_mass

    if not files_exist: # i.e. if files already existed, they already have 'CO_ppmvd' and 'NO_ppmvd'
        sec_stage_DF['CO_ppmvd'] = correctNOx(sec_stage_DF['X_CO'].values, sec_stage_DF['X_H2O'].values, sec_stage_DF['X_O2'].values)
        sec_stage_DF['NO_ppmvd'] = correctNOx(sec_stage_DF['X_NO'].values, sec_stage_DF['X_H2O'].values, sec_stage_DF['X_O2'].values)
        sec_stage_DF['conc_NO'] = sec_stage_DF['density_mole'].values * sec_stage_DF['X_NO'].values
        if write_df == True:
            sec_stage_DF.to_parquet(sec_stage_file)
            rate_DF.to_parquet(sec_stage_rate_file)        
    ign_idx = sec_stage_DF['X_OH'].values.argmax()
    ign_row = sec_stage_DF.iloc[ign_idx]
    if CO_constraint_active: 
        cons_idx = get_cons_idx(sec_stage_DF['CO_ppmvd'], CO_constraint)
        cons_row = sec_stage_DF.iloc[cons_idx] if cons_idx >= 0 else pd.Series(zero_dict)    
    else: # if CO_constraint_active = False, then just take first item after end time
        try:
            cons_row = sec_stage_DF.loc[sec_stage_DF['age'] >= tau_sec].iloc[0]
        except Exception: 
            cons_row = sec_stage_DF.iloc[-1]

    #* Compile output data into DataFrame
    zero_dict = dict(zip(sec_stage_DF.columns, np.zeros(len(sec_stage_DF.columns))-1000))
    columns = ['tau_ent_main (s)', 'tau_ent_sec (s)', 'tau_ent_ratio', 'phi_global', 'phi_main', 'phi_jet', 'phi_jet_norm', 'NO_ppmvd_constraint', 'CO_ppmvd_constraint', 'T_constraint', 'tau_sec_required (s)', 'phi_constraint', 'tau_ign_OH (s)', 'NO_ppmvd_ign', 'CO_ppmvd_ign', 'T_ign', 'phi_ign', 'X_CH4_ign', 'X_CO2_ign', 'X_O2_ign', 'X_OH_ign', 'phi_init', 'T_init', 'T_max']
    data = np.hstack([tau_ent_main, tau_ent_sec, tau_ent_main/tau_ent_sec, phi_global, phi_main, phi_jet, phi_jet_norm, cons_row[['NO_ppmvd', 'CO_ppmvd', 'T', 'age', 'phi']], ign_row[['age', 'NO_ppmvd', 'CO_ppmvd', 'T', 'phi', 'X_CH4', 'X_CO2', 'X_O2', 'X_OH']], sec_stage_DF['phi'].iloc[0], sec_stage_DF['T'].iloc[0], sec_stage_DF['T'].max()])
    out_DF = pd.DataFrame(columns=columns,data=data.reshape((1,len(data))))
    out_DF['reactor_file'] = sec_stage_file
    out_DF['rate_file'] = sec_stage_rate_file

    return out_DF

def finite_entrainment_mappable(arg_tuple, **kwargs): 
    return finite_entrainment(*arg_tuple, **kwargs)    

def get_cons_idx(COHistory, CO_constraint) -> int:
#! if all(COHistory < CO_constraint), get error: index -1 is out of bounds for axis 0 with size 0
    # idx_list = (COHistory > CO_constraint) if not all(COHistory < CO_constraint) else (len(COHistory) - 1)
    try:
        return max(min(np.arange(len(COHistory))[COHistory > CO_constraint][-1] + 1, len(COHistory)-1), 0)
    except:
        return -1
