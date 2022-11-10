import os
import sys

import cantera as ct

file_path = os.path.abspath(__file__)
sys.path.insert(1, os.path.dirname(os.path.dirname(file_path)))

from lime.utils import parse_args, solve_mass_flow_symbolic

if __name__ == "__main__":
    args = parse_args()
    print(solve_mass_flow_symbolic(args))
