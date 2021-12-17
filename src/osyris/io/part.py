# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)
import numpy as np
import os
from .reader import Reader, ReaderKind
from .. import config
from ..core import Array
from . import utils


class PartReader(Reader):
    def __init__(self):
        super().__init__(kind=ReaderKind.PART)

    def initialize(self, meta, select):
        # Read the number of variables from the hydro_file_descriptor.txt
        # and select the ones to be read if specified by user
        fname = os.path.join(meta["infile"], "part_file_descriptor.txt")
        try:
            descriptor = np.loadtxt(fname, dtype=str, delimiter=",")
        except IOError:
            return

        scaling = utils.get_spatial_scaling(meta["unit_d"], meta["unit_l"],
                                            meta["unit_t"], meta["scale"])

        part_units = {
            'position_x': scaling,
            'position_y': scaling,
            'position_z': scaling
        }

        for i in range(len(descriptor)):
            key = descriptor[i, 1].strip()
            read = True
            if isinstance(select, bool):
                read = select
            elif key in select:
                if isinstance(select[key], bool):
                    read = select[key]
            self.variables[key] = {
                "read":
                read,
                "type":
                descriptor[i, 2].strip(),
                "buffer":
                None,
                "pieces": {},
                "unit":
                part_units[key] if key in part_units else config.get_unit(
                    key, meta["unit_d"], meta["unit_l"], meta["unit_t"])
            }
        self.initialized = True

    def read_header(self, info):
        self.offsets["i"] += 2
        self.offsets["n"] += 2
        [nparticles] = utils.read_binary_data(fmt="i",
                                              content=self.bytes,
                                              offsets=self.offsets)
        for i in range(5):
            self.offsets["b"] += utils.skip_binary_line(content=self.bytes,
                                                        offsets=self.offsets)
        for item in self.variables.values():
            if item["read"]:
                npieces = len(item["pieces"])
                item["pieces"][npieces] = Array(values=np.array(
                    utils.read_binary_data(fmt="{}{}".format(nparticles, item["type"]),
                                           content=self.bytes,
                                           offsets=self.offsets)) *
                                                item["unit"].magnitude,
                                                unit=1.0 * item["unit"].units)
            else:
                self.offsets[item["type"]] += nparticles
                self.offsets["n"] += 1
        info["nparticles"] += nparticles

    def allocate_buffers(self, ngridmax, twotondim):
        return

    def read_variables(self, ncache, ind, ilevel, cpuid, info):
        return

    def make_conditions(self, select, ncache):
        return {}

    def read_footer(self, *args, **kwargs):
        return

    def step_over(self, ncache, twotondim, ndim):
        return
