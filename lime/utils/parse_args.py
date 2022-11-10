import argparse
import os
import math
from warnings import warn

import yaml

from .cantera_utils import calculate_f_stoich

d = os.path.dirname
default_config_path = os.path.join(
    d(d(d(os.path.abspath(__file__)))), "config", "default_config.yaml"
)


def parse_args(config=default_config_path, **kwargs):
    parser = argparse.ArgumentParser(description="LiME")
    
    # parse config file first, then add arguments from config file
    parser.add_argument("--config", default=config)
    args, unknown = parser.parse_known_args()
    config = yaml_config_hook(config)

    # add arguments from `config` dictionary into parser, handling boolean args too
    bool_args = []
    float_args = [
        "phi_global",
        "tau_global_ms",
        "phi_main",
        "phi_secondary",
        "total_mass",
        "ox_split",
        "tau_main_ms",
        "tau_secondary_ms",
    ]

    for k, v in config.items():
        if k == "config":  # already added config earlier, so skip
            continue
        v = kwargs.get(k, v)
        if k in bool_args:
            parser.add_argument(f"--{k}", default=v, type=bool_or_str)
        elif k in float_args:
            parser.add_argument(f"--{k}", default=v, type=none_or_float)
        else:
            parser.add_argument(f"--{k}", default=v, type=type(v))
    for k, v in kwargs.items():
        if k not in config:
            parser.add_argument(f"--{k}", default=v, type=type(v))

    # parse added arguments
    args, _ = parser.parse_known_args()
    for k, v in vars(args).items():
        # Convert boolean args that happen to be str to bool
        if k in bool_args and isinstance(v, str):
            if v.lower() in ["yes", "no", "true", "false", "none"]:
                exec(f'args.{k} = v.lower() in ["yes", "true"]')

    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path, exist_ok=True)

    args.f_stoich = calculate_f_stoich(args.main_fuel, args.oxidizer, mech=args.mech)

    return args


def yaml_config_hook(config_file):
    """
    Custom YAML config loader, which can include other yaml files (I like using config files
    insteaad of using argparser)
    """

    # load yaml files in the nested 'defaults' section, which include defaults for experiments
    with open(config_file) as f:
        cfg = yaml.safe_load(f)
        for d in cfg.get("defaults", []):
            config_dir, cf = d.popitem()
            ext = ".yaml" if len(os.path.splitext(cf)) == 1 else ""
            cf = os.path.join(os.path.dirname(config_file), config_dir, cf + ext)
            if not os.path.exists(cf):
                cf = os.path.basename(cf)
                repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
                cf = os.path.join(repo_root, "config", config_dir, cf + ext)
            with open(cf) as f:
                l = yaml.safe_load(f)
                cfg = dict(l, **cfg)

    if "defaults" in cfg.keys():
        del cfg["defaults"]

    return cfg


def none_or_other(value):
    if isinstance(value, str) and value.lower() in ["none", "null"]:
        return None
    elif isinstance(value, str) and value.lower() in ["inf", "infty", "infinity"]:
        return math.inf
    return value


none_or_float = none_or_other  # alias


def bool_or_str(value):
    # https://stackoverflow.com/a/48295546/9045125
    value = none_or_other(value)
    if value is None:
        return False
    elif value.lower() in ["yes", "y", "true", "t", "1"]:
        return True
    elif value.lower() in ["no", "n", "null", "false", "f", "0"]:
        return False
    return value
