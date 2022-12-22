import numpy as np
from sympy import Eq, Symbol, symbols
from sympy.solvers import solve

# Potential tests:
# Φ_secondary=inf, ox_split=None == Φ_secondary=None, ox_split=1


def solve_mass_flow_symbolic(args):
    """
    Calculate the mass flow rates of fuel and oxidizer for the main and secondary burners using the
    global equivalence ratio, main burner equivalence ratio, secondary burner equivalence ratio,
    total mass flow rate, and the stoichiometric fuel-to-oxidizer ratio. The mass flow rates are
    calculated either by specifying the fraction of oxidizer that goes to the secondary burner or
    by specifying the secondary burner equivalence ratio.
    
    Parameters:
    - args: argparse.Namespace object with the following attributes:
        - phi_global: float, global equivalence ratio
        - phi_main: float, main burner equivalence ratio
        - phi_secondary: float, secondary burner equivalence ratio
        - total_mass: float, total mass flow rate
        - f_stoich: float, stoichiometric fuel-to-oxidizer ratio
        - ox_split: float, fraction of oxidizer that goes to the secondary burner
    
    Returns:
    - mfm: float, mass flow rate of fuel for the main burner
    - mom: float, mass flow rate of oxidizer for the main burner
    - mfs: float, mass flow rate of fuel for the secondary burner
    - mos: float, mass flow rate of oxidizer for the secondary burner
    """    
    Φ_global, Φ_main, Φ_sec = symbols(
        "Φ_global Φ_main Φ_sec", real=True, nonnegative=True
    )
    mfm, mom = symbols("mdot_fuel_main mdot_ox_main", real=True, nonnegative=True)
    mfs, mos = symbols("mdot_fuel_sec mdot_ox_sec", real=True, nonnegative=True)
    m_total = symbols("mdot_total", real=True, nonnegative=True)
    fs = symbols("f_stoich", real=True, nonnegative=True)
    ox_split = Symbol("ox_split", real=True, nonnegative=True)

    # These equations always hold regardless of which flow rates are known
    pg = Eq(((mfm + mfs) / (mom + mos)) / fs, Φ_global)
    pm = Eq(Φ_main, (mfm / mom) / fs)
    total_mass = Eq(m_total, mfm + mfs + mom + mos)

    # Depending on what variables are provided, solve for mfm, mom, mfs and mos.
    if args.ox_split is not None:
        eq4 = Eq(ox_split, mom / (mom + mos))
        assert args.phi_secondary is None, "Cannot specify both ox_split and Φ_sec"
    elif args.phi_secondary is not None:
        # check if args.phi_secondary is infinity
        if args.phi_secondary == np.inf or args.phi_secondary == float("inf"):
            args.ox_split = 1
            ox_split = Symbol("ox_split", real=True, nonnegative=True)
            eq4 = Eq(ox_split, mom / (mom + mos))
        else:
            eq4 = Eq(Φ_sec, (mfs / mos) / fs)
    else:
        raise NotImplementedError(
            f"We only know how to calculate mass flow rates if either Φ_sec or ox_split is provided."
        )

    # Solve for the unknowns
    out_dict = solve([pg, pm, total_mass, eq4], [mfm, mom, mfs, mos], dict=True)[0]

    # Calculate float values for the outputs
    input_dict = {
        Φ_global: args.phi_global,
        Φ_main: args.phi_main,
        m_total: args.total_mass,
        fs: args.f_stoich,
    }
    if args.ox_split is not None:
        input_dict[ox_split] = args.ox_split
    elif args.phi_secondary is not None:
        input_dict[Φ_sec] = args.phi_secondary
    mfm, mom, mfs, mos = [expr.evalf(subs=input_dict) for expr in out_dict.values()]
    return mfm, mom, mfs, mos