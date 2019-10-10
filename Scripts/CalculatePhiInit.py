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
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.collections import LineCollection

def multiline(xs, ys, c, ax=None, **kwargs):
    """Plot lines with different colorings

    Parameters
    ----------
    xs : iterable container of x coordinates
    ys : iterable container of y coordinates
    c : iterable container of numbers mapped to colormap
    ax (optional): Axes to plot on.
    kwargs (optional): passed to LineCollection

    Notes:
        len(xs) == len(ys) == len(c) is the number of line segments
        len(xs[i]) == len(ys[i]) is the number of points for each line (indexed by i)

    Returns
    -------
    lc : LineCollection instance.
    """

    # find axes
    ax = plt.gca() if ax is None else ax

    # create LineCollection
    segments = [np.column_stack([x, y]) for x, y in zip(xs, ys)]
    lc = LineCollection(segments, **kwargs)

    # set coloring of line segments
    #    Note: I get an error if I pass c as a list here... not sure why.
    lc.set_array(np.asarray(c))

    # add lines to axes and rescale 
    #    Note: adding a collection doesn't autoscalee xlim/ylim
    ax.add_collection(lc)
    ax.autoscale()
    return lc

#* CONSTANTS
P = 25*ct.one_atm
T_fuel = 300 
T_air = 650
mech = 'gri30.xml'
fs_CH4 = 0.058387057492574147288255659304923028685152530670166015625
phi_global = 0.635
phi_main = 0.3719
phi_sec_norm = np.linspace(0.3,0.9,200)
phi_sec = np.divide(phi_sec_norm, (1 - phi_sec_norm))
results = []
for j, ps in enumerate(phi_sec):
    mfm, mam, mfs, mas = ctools.calculate_flowRates(phi_global, phi_main, ps) 
    if any(map(lambda x: x < 0, [mfm, mam, mfs, mas])):
        continue
    tau_ent_ratio = np.linspace(0.1, 10, 500)
    tau_ent_main = np.zeros(tau_ent_ratio.shape) + 1.5*0.001 
    tau_ent_sec = np.divide(tau_ent_main,tau_ent_ratio)
    phi_init = np.zeros(tau_ent_ratio.shape)
    for i, (tem, tes) in enumerate(zip(tau_ent_main, tau_ent_sec)):
        mdot_fuel_init = mfm/tem + mfs/tes
        mdot_air_init = mam/tem + mas/tes
        phi_init[i] = (mdot_fuel_init/mdot_air_init)/fs_CH4
        results.append([phi_sec_norm[j], ps, tem, tes, tau_ent_ratio[i], phi_init[i]])

results_df = pd.DataFrame(columns=['phi_sec_norm', 'phi_sec', 'tau_ent_main', 'tau_ent_sec', 'tau_ent_ratio', 'phi_init'], data=results)

cmap = matplotlib.cm.get_cmap('viridis')
cVar='phi_sec_norm'
y='phi_init'
x='tau_ent_sec'

f, ax = plt.subplots()
plt.cla()
vMax = results_df[cVar].max() if cVar != 'NO_ppmvd_constraint' else 50
norm = matplotlib.colors.Normalize(vmin=results_df[cVar].min(), vmax=vMax)    
c_vals = []
x_data = []
y_data = []
for cVal, data in results_df.groupby(cVar):
#         ax.plot(data[x].values, data[y].values, c=cmap(norm(cVal)), linewidth=2.0)
#         data.sort_values(by=x,inplace=True)
    c_vals.append(cVal)
    x_data.append(data[x].values)
    y_data.append(data[y].values)
#     scat = ax.scatter(results_df[x].values,results_df[y].values, c=results_df[cVar].values, vmin=results_df[cVar].min(), vmax=vMax, cmap=cmap)
lc = multiline(np.array(x_data), np.array(y_data), np.array(c_vals), cmap=cmap, lw=2)
plt.xlabel(x)
plt.ylabel(y)
cb = plt.colorbar(lc)
cb.set_label(cVar)
# plt.grid()
plt.show()
# tau_ent_ratio = 1.515
