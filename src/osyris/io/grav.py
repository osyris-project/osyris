# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)
import os
from .reader import Reader, ReaderKind
from . import utils


class GravReader(Reader):

    def __init__(self):
        super().__init__(kind=ReaderKind.AMR)

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

        self.descriptor_to_variables(descriptor=descriptor,
                                     meta=meta,
                                     units=units,
                                     select=select)
        self.initialized = True

    def read_header(self, info):
        self.offsets["i"] += 4
        self.offsets["n"] += 4

    def read_domain_header(self):
        self.offsets['n'] += 2
        self.offsets['i'] += 2
