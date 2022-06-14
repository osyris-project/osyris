# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)

from enum import Enum
import numpy as np
from . import utils
from ..core import Array


class ReaderKind(Enum):
    AMR = 0
    SINK = 1
    PART = 2


class Reader():

    def __init__(self, kind=None):
        self.variables = {}
        self.offsets = {}
        self.meta = {}
        self.bytes = None
        self.initialized = False
        self.kind = kind

    def descriptor_to_variables(self, descriptor, meta, units, select):
        drop_others = False
        if isinstance(select, dict):
            for key, value in select.items():
                if value is True:
                    drop_others = True

        for key in descriptor:
            read = True
            if isinstance(select, bool):
                read = select
            elif key in select:
                if isinstance(select[key], bool):
                    read = select[key]
            elif drop_others:
                read = False
            self.variables[key] = {
                "read": read,
                "type": descriptor[key],
                "buffer": None,
                "pieces": {},
                "unit": units[key]
            }

    def allocate_buffers(self, ncache, twotondim):
        for item in self.variables.values():
            if item["read"]:
                item["buffer"] = Array(values=np.empty([ncache * twotondim],
                                                       dtype=np.dtype(item["type"])),
                                       unit=item["unit"].units)

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
                item["buffer"]._array[ind * ncache:(ind + 1) * ncache] = np.array(
                    utils.read_binary_data(
                        fmt="{}{}".format(ncache, item["type"]),
                        content=self.bytes,
                        offsets=self.offsets)) * item["unit"].magnitude
            else:
                self.offsets[item["type"]] += ncache
                self.offsets["n"] += 1

    def make_conditions(self, select):
        conditions = {}
        if not isinstance(select, bool):
            for key, func in select.items():
                if not isinstance(func, bool):
                    if key in self.variables:
                        conditions[key] = func(self.variables[key]["buffer"])
        return conditions

    def read_footer(self, *args, **kwargs):
        return

    def step_over(self, ncache, twotondim, ndim):
        self.offsets['d'] += ncache * twotondim * len(self.variables)
        self.offsets['n'] += twotondim * len(self.variables)
