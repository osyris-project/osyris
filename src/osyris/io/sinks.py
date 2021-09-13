# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np
import os
from ..core import Array, Dataset
from .. import config
from .. import units
from . import utils


def read_sinks(nout, path, code_units):

    sink_file = utils.generate_fname(nout, path, ftype="sink", cpuid=0, ext=".csv")
    print(sink_file)

    sink_data = np.loadtxt(sink_file, delimiter=',', skiprows=2)
    # print(sink_data)

    with open(sink_file, 'r') as f:
        key_list = f.readline()
        unit_combinations = f.readline()

    key_list = key_list.lstrip(' #').rstrip('\n').split(',')
    unit_combinations = unit_combinations.lstrip(' #').rstrip('\n').split(',')

    print(key_list)
    print(unit_combinations)

    # Parse units
    unit_list = []
    for u in unit_combinations:
        print(u)
        m = code_units['ud'] * code_units['ul']**3 * units.g
        l = code_units['ul'] * units.cm
        t = code_units['ut'] * units.s
        if u == '1':
            unit_list.append(1.0 * units.dimensionless)
        else:
            unit_list.append(eval(u.replace(' ', '*')))

    print(unit_list)

    sinks = Dataset()

    for i, (key, unit) in enumerate(zip(key_list, unit_list)):
        sinks[key] = Array(values=sink_data[:, i] * unit.magnitude, unit=unit.units)

    print(sinks)

    return sinks

    # self.data = Dataset()

    # if scale is None:
