# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

# flake8: noqa

# Import the config from "/home/user/.osyris/config if it exists.
# If it doesn't, try to create one by copying the default from the source.
# If that fails, just load the default.
import os
import sys
from pint import UnitRegistry

units = UnitRegistry(system="cgs")

config_dir = os.path.join(os.path.expanduser("~"), ".osyris")
sys.path.append(config_dir)
try:
    import config_osyris as config
except ImportError:
    from shutil import copyfile
    this_dir = os.path.dirname(os.path.abspath(__file__))
    if not os.path.exists(config_dir):
        os.mkdir(config_dir)
    copyfile(os.path.join(this_dir, "config.py"),
             os.path.join(config_dir, "config_osyris.py"))
    try:
        import config_osyris as config
    except ImportError:
        from . import config

config.additional_units(units)

from .io import load
from .plot import histogram, plane
from .core import Array, Dataset, Plot
