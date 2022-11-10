import os
import sys

file_path = os.path.abspath(__file__)
sys.path.insert(1, os.path.dirname(os.path.dirname(file_path)))
