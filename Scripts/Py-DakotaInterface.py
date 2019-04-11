import sys
sys.path.insert(0, "../")
import os.path
import math

from CanteraTools import *
from FiniteEntrainment import *

import dakota.interfacing as di

#parsing Dakota param File & Creating Variables
miliseconds = 0.001
params, results = di.read_parameters_file()
#continuous_vars = [ params['tau_ent_sec'], params['tau_ent_cf'] ]
tau_ent_sec = params["tau_ent_sec"]*miliseconds
tau_ent_main = params["tau_ent_cf"]*miliseconds

tau_sec = 0;
NOs = 0;
COs = 0;

#running case
(tau_sec, NOs, COs) = runCase(tau_ent_main, tau_ent_sec, toPickle = False)
print('Reaction Time until constraint: ' + str(tau_sec*1e3) + ' ms')
print('Final Overall 15% O2 Corrected NO concentration: ' + str(NOs) + ' ppm')

#writing to output of function
for i, r in enumerate(results.responses()):
    if r.asv.function:
        r.function = NOs
#need no gradients, hessians
results.write()