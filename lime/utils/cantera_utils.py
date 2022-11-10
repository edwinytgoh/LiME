import cantera as ct


def premix(phi, fuel, ox, mech='gri30.xml', P=25*101325, T_fuel=300, T_ox=650, M_total=1):

    ox_gas = ct.Solution(mech)
    ox_gas.TPX = [T_ox, P, ox]
    fuel_gas = ct.Solution(mech)
    fuel_gas.TPX = T_fuel, P, fuel
    
    # Temporary ThermoPhase object to get mass flow rates:
    temp = ct.Solution('gri30.xml')
    temp.set_equivalence_ratio(phi, fuel, ox)
    mdot_fuel = M_total * sum(temp[fuel.keys()].Y)
    mdot_ox = M_total * sum(temp[ox.keys()].Y)

    # Output mixer gas: 
    return mix([fuel_gas, ox_gas], [mdot_fuel, mdot_ox], P=P) 


def calculate_f_stoich(fuel, oxidizer, basis="mole", mech="gri30.yaml"):
    temp_gas = ct.Solution(mech)
    f_stoich_air_fuel = temp_gas.stoich_air_fuel_ratio(fuel, oxidizer, basis)
    return 1 / f_stoich_air_fuel


def mix(streams, mdots, mech="gri30.xml", P=25*101325):
    # Create mixer gas: 
    mixerGas = ct.Solution(mech) 
    mixerGas.TPX = [300, P, 'H2:1']
    
    # Create reactor with CHEMISTRY DISABLED: 
    mixer = ct.ConstPressureReactor(mixerGas) 
    mixer.chemistry_enabled = False # disable chemistry 
    
    # For each stream (and given mass flow rate), connect mass flow controller to mixer: 
    mfcs = [ct.MassFlowController(ct.Reservoir(streams[i]), mixer, mdot=mdots[i]) for i in range(0,len(streams))]
    
    exhaust = ct.Reservoir(mixerGas) 
    exhaust_mfc = ct.MassFlowController(mixer, exhaust, mdot=sum(mdots)) 
    
    rn = ct.ReactorNet([mixer]) 
    rn.advance_to_steady_state() 
    
    return mixer.thermo