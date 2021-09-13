import numpy as np
from . import utils
from ..core import Array


class Reader():
    def __init__(self, code_units=None):
        self.variables = {}
        self.offsets = {}
        self.meta = {}
        self.bytes = None
        self.code_units = code_units
        # self.select = None
        return

    def allocate_buffers(self, ngridmax, twotondim):
        for item in self.variables.values():
            if item["read"]:
                item["buffer"] = Array(values=np.zeros([ngridmax, twotondim],
                                                       dtype=np.dtype(item["type"])),
                                       unit=1.0 * item["unit"].units)

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
                item["buffer"]._array[:ncache, ind] = np.array(
                    utils.read_binary_data(
                        fmt="{}{}".format(ncache, item["type"]),
                        content=self.bytes,
                        offsets=self.offsets)) * item["unit"].magnitude
            else:
                self.offsets[item["type"]] += ncache
                self.offsets["n"] += 1

    def make_conditions(self, select, ncache):
        conditions = {}
        for key, func in select.items():
            if not isinstance(func, bool):
                if key in self.variables:
                    conditions[key] = func(self.variables[key]["buffer"][:ncache, :])
        return conditions

    def read_footer(self, *args, **kwargs):
        return

    def step_over(self, ncache, twotondim, ndim):
        self.offsets['d'] += ncache * twotondim * len(self.variables)
        self.offsets['n'] += twotondim * len(self.variables)