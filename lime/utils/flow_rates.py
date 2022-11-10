import numpy as np
from sympy import Eq, Symbol, symbols
from sympy.solvers import solve

# Potential tests:
# phi_secondary=inf, ox_split=None == phi_secondary=None, ox_split=1


def solve_mass_flow_symbolic(args):
    phi_global, phi_main, phi_sec = symbols(
        "phi_global phi_main phi_sec", real=True, nonnegative=True
    )
    mfm, mom = symbols("mdot_fuel_main mdot_ox_main", real=True, nonnegative=True)
    mfs, mos = symbols("mdot_fuel_sec mdot_ox_sec", real=True, nonnegative=True)
    m_total = symbols("mdot_total", real=True, nonnegative=True)
    fs = symbols("f_stoich", real=True, nonnegative=True)
    ox_split = Symbol("ox_split", real=True, nonnegative=True)

    # These equations always hold regardless of which flow rates are known
    pg = Eq(((mfm + mfs) / (mom + mos)) / fs, phi_global)
    pm = Eq(phi_main, (mfm / mom) / fs)
    total_mass = Eq(m_total, mfm + mfs + mom + mos)

    # Depending on what variables are provided, solve for mfm, mom, mfs and mos.
    if args.ox_split is not None:
        eq4 = Eq(ox_split, mom / (mom + mos))
        assert args.phi_secondary is None, "Cannot specify both ox_split and phi_sec"
    elif args.phi_secondary is not None:
        # check if args.phi_secondary is infinity
        if args.phi_secondary == np.inf or args.phi_secondary == float("inf"):
            args.ox_split = 1
            ox_split = Symbol("ox_split", real=True, nonnegative=True)
            eq4 = Eq(ox_split, mom / (mom + mos))
        else:
            eq4 = Eq(phi_sec, (mfs / mos) / fs)
    else:
        raise NotImplementedError(
            f"We only know how to calculate mass flow rates if either phi_sec or ox_split is provided."
        )

    # Solve for the unknowns
    out_dict = solve([pg, pm, total_mass, eq4], [mfm, mom, mfs, mos], dict=True)[0]

    # Calculate float values for the outputs
    input_dict = {
        phi_global: args.phi_global,
        phi_main: args.phi_main,
        m_total: args.total_mass,
        fs: args.f_stoich,
    }
    if args.ox_split is not None:
        input_dict[ox_split] = args.ox_split
    elif args.phi_secondary is not None:
        input_dict[phi_sec] = args.phi_secondary
    output_values = [expr.evalf(subs=input_dict) for expr in out_dict.values()]
    return output_values


def calculate_mass_flow_rates(args):
    pass
