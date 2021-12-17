# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)
import os
from .reader import Reader, ReaderKind
from .. import config
from . import utils


class GravReader(Reader):
    def __init__(self):
        super().__init__(kind=ReaderKind.AMR)

    def initialize(self, meta, select):
        fname = utils.generate_fname(meta["nout"], meta["path"], ftype="grav", cpuid=1)
        # Check if self-gravity files exist
        if not os.path.exists(fname):
            return
        # Add gravity fields
        descriptor = {"potential": "d"}
        for n in range(meta["ndim"]):
            descriptor["acceleration_" + "xyz"[n]] = "d"
        # Now add to the list of variables to be read
        for key in descriptor:
            read = True
            if isinstance(select, bool):
                read = select
            elif key in select:
                if isinstance(select[key], bool):
                    read = select[key]
            self.variables[key] = {
                "read": read,
                "type": descriptor[key],
                "buffer": None,
                "pieces": {},
                "unit": config.get_unit(key, meta["unit_d"], meta["unit_l"],
                                        meta["unit_t"])
            }
        self.initialized = True

    def read_header(self, info):
        self.offsets["i"] += 4
        self.offsets["n"] += 4

    def read_domain_header(self):
        self.offsets['n'] += 2
        self.offsets['i'] += 2
