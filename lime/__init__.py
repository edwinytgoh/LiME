import cantera as ct
import os
DATA_DIR = os.path.abspath('../data')
MECH_DIR = os.path.join(DATA_DIR, 'mechanisms')
ct.add_directory(MECH_DIR) # for custom mechanisms
