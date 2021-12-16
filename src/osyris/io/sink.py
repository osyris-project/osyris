# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np
import os
from ..core import Array, Datagroup
from .reader import ReaderKind
from .. import units
from . import utils


class SinkReader:
    def __init__(self):
        self.kind = ReaderKind.SINK
        self.initialized = False

    def initialize(self, meta, select, ramses_ism):
        sink = Datagroup()
        if select is False:
            return sink
        sink_file = utils.generate_fname(meta["nout"],
                                         meta["path"],
                                         ftype="sink",
                                         cpuid=0,
                                         ext=".csv")
        if not os.path.exists(sink_file):
            return

        try:
            if ramses_ism:
                print("Reading "+sink_file)
                sink_data = np.loadtxt(sink_file, delimiter=',', skiprows=0)  # do not skip rows
                print(sink_data)
            else:
                sink_data = np.loadtxt(sink_file, delimiter=',', skiprows=2)
        except StopIteration:
            # This is an empty sink file
            return sink

        if ramses_ism:
            key_list = ["number","mass","dmf","x","y","z","vx","vy","vz","period","lx","ly","lz","acc_rate","acc_lum","age","int_lum","Teff"]
            Q = units.Quantity
            units.define("solar_luminosity = 3.83e26 * watt = lsol")  #  defining lsol units
            unit_list = [1.*units.dimensionless, 1.*units.msun, 1.*units.msun, 1.*units.au, 1.*units.au, 1.*units.au, 1.*units.cmps, 1.*units.cmps, 1.*units.cmps, 1.*units.year,
            1.*units.dimensionless, 1.*units.dimensionless, 1.*units.dimensionless, 1.*units.msun/units.year, 1.*Q("lsol"), 1.*units.year, 1.*Q("lsol"), 1.*units.K]
        else:
            with open(sink_file, 'r') as f:
                key_list = f.readline()
                unit_combinations = f.readline()

            key_list = key_list.lstrip(' #').rstrip('\n').split(',')
            unit_combinations = unit_combinations.lstrip(' #').rstrip('\n').split(',')

            # Parse units
            unit_list = []
            for u in unit_combinations:
                m = meta['unit_d'] * meta['unit_l']**3 * units.g  # noqa: F841
                l = meta['unit_l'] * units.cm  # noqa: F841, E741
                t = meta['unit_t'] * units.s  # noqa: F841
                if u == '1':
                    unit_list.append(1.0 * units.dimensionless)
                else:
                    unit_list.append(eval(u.replace(' ', '*')))

        sink = Datagroup()
        for i, (key, unit) in enumerate(zip(key_list, unit_list)):
            print(i, key, unit)
            sink[key] = Array(values=sink_data[:, i] * unit.magnitude, unit=unit.units)
            if not ramses_ism and unit_combinations[i] == 'l':
                sink[key] = sink[key].to(meta["scale"])
        utils.make_vector_arrays(sink, ndim=meta["ndim"])
        return sink
