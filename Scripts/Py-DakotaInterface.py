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

params, results = di.read_parameters_file()
phi_sec = params["phi_sec"]
phi_main = params["phi_main"]
tau_res_sec = params["T_res_sec"]
phi_main = round(phi_main,3)

(mfm,mam,mfs,mas) = solveMass_PhiJet(phi_overall,phi_main, phi_sec)
mmain = mfm+mam
msec = mfs+mas
Msratio = msec/(msec+mmain)
MdotsRatio =(msec/tau_ent_sec)/((msec/tau_ent_sec)+(mmain/tau_ent_main)) 

NOx = 0; Temp = 0; #COs = 0;

#running case
(out_dF) = finite_entrainment(phi_overall, phi_main, tau_res_sec, tau_ent_main, tau_ent_sec, phi_jet=phi_sec,tau_global=20,write_df='False') #input everything in miliseconds
Temp = out_dF["T_Max"]; NOs = out_dF['NO_ppmvd_constraint']; #COs = out_dF['CO_ppmvd_constraint'], for the weighted optimization

#writing to output of function
for i, r in enumerate(results.responses()):
    if r.asv.function:
        if i == 0:
            r.function = NOs
        if i == 1:
            r.function = Temp
        if i == 2:
            r.function = Msratio
        if i == 3:
            r.function = MdotsRatio
#need no gradients, hessians
results.write()
