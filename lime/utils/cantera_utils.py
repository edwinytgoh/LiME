import cantera as ct
import numpy as np


def calculate_f_stoich(fuel, oxidizer, basis="mole", mech="gri30.yaml"):
    """
    Calculate the stoichiometric fuel/air ratio for a given fuel and oxidizer.

    Args:
        fuel (Dict[str, float]): Fuel composition
        oxidizer (Dict[str, float]): Oxidizer composition
        basis (str, optional): Stoichiometric basis. Defaults to "mole".
        mech (str, optional): Chemical mechanism to use. Defaults to "gri30.yaml".

    Returns:
        f_stoich: Stoichiometric fuel/air ratio
    """
    temp_gas = ct.Solution(mech)
    f_stoich_air_fuel = temp_gas.stoich_air_fuel_ratio(fuel, oxidizer, basis)
    return 1 / f_stoich_air_fuel


def premix(
    phi, fuel, ox, mech="gri30.yaml", P=25 * 101325, T_fuel=300, T_ox=650, M_total=1
):
    """
    Wrapper for the `mix()` function to premix a fuel and oxidizer at a given equivalence ratio.
    Args:
        phi (float): Equivalence ratio
        fuel (Dict[str, float]): Fuel composition
        ox (Dict[str, float]): Oxidizer composition
        mech (str, optional): Chemical mechanism to use. Defaults to "gri30.yaml".
        P (float, optional): Pressure in Pascals. Defaults to 25*101325.
        T_fuel (int, optional): Fuel preheat temperature in Kelvin. Defaults to 300.
        T_ox (int, optional): Oxidizer preheat temperature in Kelvin. Defaults to 650.
        M_total (int, optional): Total mass of mixture. Defaults to 1.

    Returns:
        mixed_gas: Cantera solution object for the mixed gas
    """
    ox_gas = ct.Solution(mech)
    ox_gas.TPX = [T_ox, P, ox]
    fuel_gas = ct.Solution(mech)
    fuel_gas.TPX = T_fuel, P, fuel

    # Temporary ThermoPhase object to get mass flow rates:
    temp = ct.Solution(mech)
    temp.set_equivalence_ratio(phi, fuel, ox)
    mdot_fuel = M_total * sum(temp[fuel.keys()].Y)
    mdot_ox = M_total * sum(temp[ox.keys()].Y)

    # Output mixer gas:
    return mix([fuel_gas, ox_gas], [mdot_fuel, mdot_ox], P=P)


def mix(streams, mdots, mech="gri30.yaml", P=25 * 101325):
    """
    Mix streams of gas (ct.Solution) at constant pressure to produce a single gas object.

    Args:
        streams (List[ct.Solution]): List of Cantera solution objects to be mixed
        mdots (List[float]): List of mass flow rates for each stream
        mech (str, optional): Chemical mechanism to use. Defaults to "gri30.yaml".
        P (float, optional): Pressure in Pascals. Defaults to 25*101325.

    Returns:
        mixed_gas: Cantera solution object for the mixed gas
    """
    
    # Create mixer gas:
    mixer_gas = ct.Solution(mech)
    mixer_gas.TPX = [300, P, "H2:1"]

    # Create reactor with CHEMISTRY DISABLED:
    mixer = ct.ConstPressureReactor(mixer_gas)
    mixer.chemistry_enabled = False  # disable chemistry

    # For each stream (and given mass flow rate), connect mass flow controller to mixer:
    _ = [
        ct.MassFlowController(ct.Reservoir(streams[i]), mixer, mdot=mdots[i])
        for i in range(0, len(streams))
    ]

    exhaust = ct.Reservoir(mixer_gas)
    exhaust_mfc = ct.MassFlowController(mixer, exhaust, mdot=sum(mdots))

    rn = ct.ReactorNet([mixer])
    rn.advance_to_steady_state()

    return mixer.thermo


def correct_nox(X_i, X_H2O, X_O2):
    dry_i = X_i / (1 - X_H2O)
    dry_O2 = X_O2 / (1 - X_H2O)
    corrected = dry_i * (20.9 - 15) / (20.9 - dry_O2 * 100)
    corrected_ppm = corrected * 1e6
    return corrected_ppm
