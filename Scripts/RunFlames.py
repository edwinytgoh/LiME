import sys
import os
import gc as garbage_collector
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import multiprocessing as mp 
import numpy as np 
from Utils.CanteraTools import runFlame, runMainBurner
# mb_func = lambda phi_main: runMainBurner(phi_main, tau_main=25*1e-3, filename=f"")
current_dir = os.getcwd()
def mb_func(phi_main): 
    try:
        
        vitReactor, mainBurnerDF = runMainBurner(phi_main, tau_main=25*1e-3, filename=os.path.join(current_dir, 'Flames', f"phi_main_{phi_main:.4f}_GRI30_25atm.pickle"))
        del vitReactor
        del mainBurnerDF        
        print(f"phi_main: {phi_main:.3f} done\n")        
    except:
        print(f"phi_main: {phi_main:.3f} failed\n")
    garbage_collector.collect()

if __name__ == "__main__":
    phi_main_list = np.arange(float(sys.argv[1]), float(sys.argv[2]), 0.001)
    N = mp.cpu_count() if mp.cpu_count() > 15 else np.round(mp.cpu_count()*0.75).astype(int)
    print(f"Creating a pool of size {N}...")
    pool = mp.Pool(N)
    out_list = list(pool.map(mb_func, phi_main_list))
    print(f"Done running. Closing pool...")
    pool.close()
    print(f"Done!")
