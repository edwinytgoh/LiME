import os
import sys

import cantera as ct

file_path = os.path.abspath(__file__)
sys.path.insert(1, os.path.dirname(os.path.dirname(file_path)))

from lime.utils import parse_args, solve_mass_flow_symbolic
from lime.utils.cantera_utils import correct_nox, mix, premix

if __name__ == "__main__":
    args = parse_args()
    main_burner_gas = premix(
        args.phi_main,
        args.main_fuel,
        args.oxidizer,
        args.mech,
        args.P_atm * ct.one_atm,
        args.T_fuel_main,
        args.T_oxidizer,
    )
    secondary_gas = ct.Solution(args.mech)
    secondary_gas.TPX = (
        args.T_fuel_secondary,
        args.P_atm * ct.one_atm,
        args.secondary_fuel,
    )

    mfm, mom, mfs, mos = solve_mass_flow_symbolic(args)
    print([mfm, mom, mfs, mos])
    mixed_gas = mix([main_burner_gas, secondary_gas], [mfm + mom, mfs + mos])
    mixed_gas.equilibrate("HP")
    mixed_gas()  # print mixed gas
    corrected_CO = correct_nox(mixed_gas["CO"].X, mixed_gas["H2O"].X, mixed_gas["O2"].X)
    corrected_NO = correct_nox(mixed_gas["NO"].X, mixed_gas["H2O"].X, mixed_gas["O2"].X)
    print(f"Corrected CO (dry, 15% O2): {corrected_CO} ppm")
    print(f"Corrected NO (dry, 15% O2): {corrected_NO} ppm")
