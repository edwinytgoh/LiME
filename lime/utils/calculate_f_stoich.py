import cantera as ct


def calculate_f_stoich(fuel, oxidizer, basis="mole", mech="gri30.yaml"):
    temp_gas = ct.Solution(mech)
    f_stoich_air_fuel = temp_gas.stoich_air_fuel_ratio(fuel, oxidizer, basis)
    return 1 / f_stoich_air_fuel
