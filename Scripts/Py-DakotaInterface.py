import sys
sys.path.insert(0, "../")
import os.path
import math

from CanteraTools import *
from FiniteEntrainment import *

import dakota.interfacing as di

#parsing Dakota param File & Creating Variables
miliseconds = 0.001
tau_ent_main = 10*miliseconds
tau_ent_sec = 1*miliseconds
phi_overall = 0.635

params, results = di.read_parameters_file()
phi_sec = params["phi_sec"]
phi_main = params["phi_main"]
tau_res_sec = param["T_res_sec"]*miliseconds


(mfm,mam,mfs,mas) = solveMass_PhiJet(phi_overall,phi_main, phi_sec)
mmain = mfm+mam
msec = mfs+mas
Msratio = msec/(msec+mmain)
MdotsRatio =(msec/tau_ent_sec)/((msec/tau_ent_sec)+(mmain/tau_ent_main)) 

COs = 0;
NOx = 0;
T = 0;

#running case
(tau_sec, NOs, COs) = runCase(tau_ent_main, tau_ent_sec, toPickle = False)

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
