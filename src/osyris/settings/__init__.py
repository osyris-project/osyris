# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

# Import the config from "/home/user/.osyris/config.py if it exists.
# If it doesn't, try to create one by copying the default from the source.
# If that fails, just load the default.
import os
import sys
from shutil import copyfile

user_config_dir = os.path.join(os.path.expanduser("~"), ".osyris")
this_dir = os.path.dirname(os.path.abspath(__file__))
if not os.path.exists(user_config_dir):
    os.mkdir(user_config_dir)
user_config_file = os.path.join(user_config_dir, "config.py")

if not os.path.exists(user_config_file):
    copyfile(os.path.join(this_dir, "config_default.py"), user_config_file)

sys.path.append(user_config_dir)

# Create a dummy class
class Config:
    pass

# Try to import list of object from config.
# If import fails from the user config file, import from the default.

# import config_osyris as config

objects = ['parameters', 'variable_units', 'additional_units', 'additional_variables']

for obj in objects:
    try:
        from config 



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

from .io import Loader
from .plot import histogram, plane
from .core import Array, Dataset, Plot
