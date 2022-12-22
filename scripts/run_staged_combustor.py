# Script to run a staged combustor model

# Example usage (from LiME repo root directory):
# python scripts/run_staged_combustor.py --config config/default_config.yaml

import os
import sys
import time

file_path = os.path.abspath(__file__)
sys.path.insert(1, os.path.dirname(os.path.dirname(file_path)))

from lime.utils import parse_args, solve_mass_flow_symbolic
from lime.utils.axial_staged_combustor import StagedCombustor

if __name__ == "__main__":
    t1 = time.time()
    args = parse_args()
    print(solve_mass_flow_symbolic(args))
    staged_combustor = StagedCombustor(args)
    staged_combustor.run_main_burner()
    staged_combustor.run_secondary_stage()
    print(f"Finished in {time.time() - t1:.2f} seconds")