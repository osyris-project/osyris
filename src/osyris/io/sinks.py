# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np
import os
from . import utils
from ..core import Dataset
from .units import get_unit


def read_sinks(nout, path):

    sink_file = utils.generate_fname(nout, path, ftype="sink", cpuid=0, ext=".csv")
    print(sink_file)

    sink_data = np.loadtxt(sink_file, delimiter=',', skiprows=2)
    # print(sink_data)

    with open(sink_file, 'r') as f:
        keys = f.readline()
        unit_combinations = f.readline()

    keys = keys.lstrip(' #').rstrip('\n').split(',')
    unit_combinations = unit_combinations.lstrip(' #').rstrip('\n').split(',')

    print(keys)
    print(unit_combinations)

    # self.data = Dataset()

    # if scale is None:
