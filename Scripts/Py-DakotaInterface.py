import sys
sys.path.insert(0, "../")
import os.path
import math

from Utils.CanteraTools import *
from Utils.finiteEntrainment import *

import dakota.interfacing as di

#parsing Dakota param File & Creating Variables
#miliseconds = 0.001
tau_ent_main = 10#*miliseconds
tau_ent_sec = 1#*miliseconds
phi_overall = 0.635
flame_lib_loc = "C:\\Users\\Hazelle\\dakota_gui_workspace\\NOxOpt3(Phim,s+Trs)\\FlameLibrary"
params, results = di.read_parameters_file()

#* Input parameters (design variables) from Dakota Input File
phi_sec = params["phi_sec"]
phi_main = params["phi_main"]
tau_res_sec = params["T_res_sec"]
phi_main = round(phi_main,3) # Round phi_main to 3 d.p. because flame library only has 3 dp

#* Calculate mass split and entrainment ratios for non-linear constraints
(mfm,mam,mfs,mas) = solveMass_PhiJet(phi_overall,phi_main, phi_sec)
mmain = mfm+mam
msec = mfs+mas
Msratio = msec/(msec+mmain) # ratio of secondary mass to total mass must be 20% or less
# MdotsRatio =(msec/tau_ent_sec)/((msec/tau_ent_sec)+(mmain/tau_ent_main)) 
constraint_violated = Msratio >= 0.2 # constraint violated if secondary mass > 20% of total mass #TODO: Change this to depend on Dakota inputs!!


#running case
if not constraint_violated:
    (out_dF) = finite_entrainment(phi_overall, phi_main, tau_res_sec, tau_ent_main, tau_ent_sec, phi_jet=phi_sec,tau_global=20,write_df='False',flame_file=os.path.join(flame_lib_loc, f"phi_main_{phi_main:.4f}_GRI30_25atm.pickle")) #input everything in miliseconds
    temp = out_dF["T_constraint"]; 
    NOx = out_dF['sys_NO_ppmvd'];
    COs = out_dF['sys_CO_ppmvd'];
else:
    temp = 0
    NOx = 1e6
    COs = 1e6

#writing to output of function
for i, r in enumerate(results.responses()):
    if r.asv.function:
        if i == 0:
            r.function = NOx
        if i == 1:
            r.function = temp
        if i == 2:
            r.function = Msratio
        if i == 3:
            r.function = COs
#need no gradients, hessians
results.write()
