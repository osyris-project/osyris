# SPDX-License-Identifier: BSD-3-Clause

import numpy as np

from ..core import Array
from . import utils


class Reader:
    def __init__(self, kind: str = None):
        self.variables = {}
        self.offsets = {}
        self.meta = {}
        self.bytes = None
        self.initialized = False
        self.kind = kind

    def descriptor_to_variables(self, descriptor, meta, units, select):
        read = {key: False for key in descriptor}
        if isinstance(select, dict):
            for key in read:
                read[key] = select.get(key, True)
        elif isinstance(select, bool):
            for key in read:
                read[key] = select
        else:
            for key in select:
                read[key] = True

        for key in descriptor:
            self.variables[key] = {
                "read": read[key],
                "type": descriptor[key],
                "buffer": None,
                "pieces": {},
                "unit": units[key],
            }

    def allocate_buffers(self, ncache, twotondim):
        for item in self.variables.values():
            if item["read"]:
                item["buffer"] = Array(
                    values=np.empty([ncache * twotondim], dtype=np.dtype(item["type"])),
                    unit=item["unit"].units,
                )

    def read_header(self, *args, **kwargs):
        return

    def read_level_header(self, *args, **kwargs):
        return

    def read_domain_header(self, *args, **kwargs):
        return

    def read_cacheline_header(self, *args, **kwargs):
        return

    def read_variables(self, ncache, ind, ilevel, cpuid, info):
        for item in self.variables.values():
            if item["read"]:
                item["buffer"]._array[ind * ncache : (ind + 1) * ncache] = (
                    np.array(
                        utils.read_binary_data(
                            fmt="{}{}".format(ncache, item["type"]),
                            content=self.bytes,
                            offsets=self.offsets,
                        )
                    )
                    * item["unit"].magnitude
                )
            else:
                self.offsets[item["type"]] += ncache
                self.offsets["n"] += 1

    def make_conditions(self, select):
        conditions = {}
        if isinstance(select, dict):
            for key, func in select.items():
                if key in self.variables:
                    conditions[key] = func(self.variables[key]["buffer"])
        return conditions

    def read_footer(self, *args, **kwargs):
        return

    def step_over(self, ncache, twotondim, ndim):
        self.offsets["d"] += ncache * twotondim * len(self.variables)
        self.offsets["n"] += twotondim * len(self.variables)
