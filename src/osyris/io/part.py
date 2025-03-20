# SPDX-License-Identifier: BSD-3-Clause
import os

import numpy as np

from ..core import Array
from . import utils
from .reader import Reader


class PartReader(Reader):
    def __init__(self):
        super().__init__(kind="part")

    def initialize(self, meta, units, select):
        self.initialized = False
        if select is False:
            return

        # Read the number of variables from the hydro_file_descriptor.txt
        # and select the ones to be read if specified by user
        fname = os.path.join(meta["infile"], "part_file_descriptor.txt")
        try:
            desc_from_file = np.loadtxt(fname, dtype=str, delimiter=",")
        except IOError:
            return

        descriptor = {
            desc_from_file[i, 1].strip(): desc_from_file[i, 2].strip()
            for i in range(len(desc_from_file))
        }

        self.descriptor_to_variables(
            descriptor=descriptor, meta=meta, units=units, select=select
        )
        self.initialized = True

    def read_header(self, info):
        if not self.initialized:
            return
        self.offsets["i"] += 2
        self.offsets["n"] += 2
        [nparticles] = utils.read_binary_data(
            fmt="i", content=self.bytes, offsets=self.offsets
        )
        for i in range(5):
            self.offsets["b"] += utils.skip_binary_line(
                content=self.bytes, offsets=self.offsets
            )
        for item in self.variables.values():
            if item["read"]:
                npieces = len(item["pieces"])
                item["pieces"][npieces] = Array(
                    values=np.array(
                        utils.read_binary_data(
                            fmt="{}{}".format(nparticles, item["type"]),
                            content=self.bytes,
                            offsets=self.offsets,
                        )
                    )
                    * item["unit"].magnitude,
                    unit=item["unit"].units,
                )
            else:
                self.offsets[item["type"]] += nparticles
                self.offsets["n"] += 1
        info["nparticles"] += nparticles

    def allocate_buffers(self, ncache, twotondim):
        return

    def read_variables(self, ncache, ind, ilevel, cpuid, info):
        return

    def make_conditions(self, select):
        return {}

    def read_footer(self, *args, **kwargs):
        return

    def step_over(self, ncache, twotondim, ndim):
        return
