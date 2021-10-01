# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

# Import the config from "/home/user/.osyris/conf.py if it exists.
# If it doesn't, try to create one by copying the default from the source.
# If that fails, just load the default.
import os
import sys
from shutil import copyfile
from . import defaults as default_config

user_config_dir = os.path.join(os.path.expanduser("~"), ".osyris")
this_dir = os.path.dirname(os.path.abspath(__file__))
if not os.path.exists(user_config_dir):
    os.mkdir(user_config_dir)
user_config_file = os.path.join(user_config_dir, "config_osyris.py")
if not os.path.exists(user_config_file):
    copyfile(os.path.join(this_dir, "defaults.py"), user_config_file)
sys.path.append(user_config_dir)

import config_osyris as user_config  # noqa: E402


# Create a dummy class
class Config:
    pass


if hasattr(user_config, 'ureg'):
    units = getattr(user_config, 'ureg')
else:
    units = getattr(default_config, 'ureg')

config = Config()

# Import list of object from user config if present, if not, load from defaults.
objects = ['parameters', 'get_unit', 'additional_units', 'additional_variables']
for obj in objects:
    if hasattr(user_config, obj):
        setattr(config, obj, getattr(user_config, obj))
    else:
        setattr(config, obj, getattr(default_config, obj))

config.additional_units()
