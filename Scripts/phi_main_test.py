import sys
sys.path.insert(0, "../")
import os
from CanteraTools import *

if __name__ == "__main__":
    [vit_reactor, main_burner_DF] = runMainBurner(0.36, (20-0.158)*1e-3)    # Mixed temperature around 591 K
    [vit_reactor, main_burner_DF] = runMainBurner(0.3719, (20-0.158)*1e-3)    # Mixed temperature around 591 K
    [vit_reactor, main_burner_DF] = runMainBurner(0.42, (20-0.158)*1e-3)    # Mixed temperature around 591 K
    print("Done!")