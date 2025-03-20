# SPDX-License-Identifier: BSD-3-Clause
import os

from . import utils
from .reader import Reader


class GravReader(Reader):
    def __init__(self):
        super().__init__(kind="mesh")

    def initialize(self, meta, units, select):
        self.initialized = False
        if select is False:
            return

        fname = utils.generate_fname(meta["nout"], meta["path"], ftype="grav", cpuid=1)
        # Check if self-gravity files exist
        if not os.path.exists(fname):
            return
        # Add gravity fields
        descriptor = {"grav_potential": "d"}
        for n in range(meta["ndim"]):
            descriptor["grav_acceleration_" + "xyz"[n]] = "d"

        self.descriptor_to_variables(
            descriptor=descriptor, meta=meta, units=units, select=select
        )
        self.initialized = True

    def read_header(self, info):
        self.offsets["i"] += 4
        self.offsets["n"] += 4

    def read_domain_header(self):
        self.offsets["n"] += 2
        self.offsets["i"] += 2
