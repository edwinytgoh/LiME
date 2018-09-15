import sys
sys.path.insert(0, "../")
from CanteraTools import *
# import matplotlib.pyplot as plt

def COLimitInd(COHistory, constraint):

    index = len(COHistory) - 1
    while COHistory[index] < constraint:
        index -= 1
    index += 1  # Go back one so within constraint
    index = min(index, len(COHistory))    # Avoid index out of bounds
    return index

def runCase(tau_ent_cf, tau_ent_sec):
    milliseconds = 0.001

    [mfm, mam, mfs, mas] = solvePhi_airSplit(0.635, 0.3719, 100, 1)
    [vit_reactor, main_burner_DF] = runMainBurner(0.3719, 19.842*milliseconds)

    secondary_gas = ct.Solution('gri30.xml')
    secondary_gas.TPX = 300, 25*ct.one_atm, {'CH4':1}
    sec_reactor = ct.ConstPressureReactor(secondary_gas)

    totalTime = 3.0*milliseconds
    interface_gas = ct.Solution('gri30.xml')
    interface_gas.TPX = 300, 25*ct.one_atm, {'AR':1}
    interface_reactor = ct.ConstPressureReactor(interface_gas, volume = 1/interface_gas.density_mass*1e-6)

    mfc_vit = ct.MassFlowController(vit_reactor, interface_reactor)
    mdot_vit = (mfm+mam)/tau_ent_cf
    mfc_vit.set_mass_flow_rate(lambda t: mdot_vit if t <= tau_ent_cf else 0)

    mfc_sec = ct.MassFlowController(sec_reactor, interface_reactor)
    mdot_sec = (mfs+mas)/tau_ent_sec
    mfc_sec.set_mass_flow_rate(lambda t: mdot_sec if t <= tau_ent_sec else 0)

    reactorNet = ct.ReactorNet([vit_reactor, sec_reactor, interface_reactor])
    # states = []
    species_NO = np.array([])
    species_CO = np.array([])
    species_O2 = np.array([])
    species_H2O = np.array([])
    # enthalpy = []

    t = np.arange(0, totalTime, 0.001*milliseconds)
    for tnow in t:
        reactorNet.advance(tnow)

        # Get the number of moles in each reactor since we're working with mole fractions
        mass_vit = (mam + mfm) - mdot_vit*min(tnow, tau_ent_cf)
        mole_vit = mass_vit/vit_reactor.thermo.mean_molecular_weight
        mass_sec = (mas + mfs) - mdot_sec*min(tnow, tau_ent_sec)
        mole_sec = mass_sec/sec_reactor.thermo.mean_molecular_weight
        mass_int = (mam + mfm) - mass_vit + (mas + mfs) - mass_sec
        mole_int = mass_int/interface_reactor.thermo.mean_molecular_weight

        temp_arr = interface_reactor.thermo.X * mole_int + vit_reactor.thermo.X * mole_vit + sec_reactor.thermo.X * mole_sec
        temp_arr /= (mole_int + mole_vit + mole_sec)
        species_NO = np.append(species_NO, temp_arr[interface_reactor.thermo.species_index('NO')])
        species_CO = np.append(species_CO, temp_arr[interface_reactor.thermo.species_index('CO')])
        species_O2 = np.append(species_O2, temp_arr[interface_reactor.thermo.species_index('O2')])
        species_H2O = np.append(species_H2O, temp_arr[interface_reactor.thermo.species_index('H2O')])
        # enthalpy.append(sec_reactor.thermo.enthalpy_mass)

    NO_corr = correctNOx(species_NO, species_H2O, species_O2)
    CO_corr = correctNOx(species_CO, species_H2O, species_O2)
    constraint_ind = COLimitInd(CO_corr, 32)

    tau_sec = t[constraint_ind]
    NO_finalcorr = NO_corr[constraint_ind]
    return tau_sec, NO_finalcorr

(tau_sec, NO_finalcorr) = runCase(2*1e-3, 0.01*1e-3)
print('Reaction Time until constraint: ' + str(tau_sec*1e3) + ' ms')
print('Final Overall 15% O2 Corrected NO concentration: ' + str(NO_finalcorr) + ' ppm')