# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)
import numpy as np
import os
from .reader import Reader, ReaderKind
from .. import config


class RtReader(Reader):
    def __init__(self):
        super().__init__(kind=ReaderKind.AMR)

    def initialize(self, meta, select):
        # Read the number of variables from the rt_file_descriptor.txt
        # and select the ones to be read if specified by user
        fname = os.path.join(meta["infile"], "rt_file_descriptor.txt")
        try:
            descriptor = np.loadtxt(fname, dtype=str, delimiter=",")
        except IOError:
            return

        for i in range(len(descriptor)):
            key = descriptor[i, 1].strip()
            read = True
            if isinstance(select, bool):
                read = select
            elif key in select:
                if isinstance(select[key], bool):
                    read = select[key]
            self.variables[key] = {
                "read": read,
                "type": descriptor[i, 2].strip(),
                "buffer": None,
                "pieces": {},
                "unit": config.get_unit(key, meta["unit_d"], meta["unit_l"],
                                        meta["unit_t"])
            }
        self.initialized = True

    def read_header(self, info):
        self.offsets["i"] += 5
        self.offsets["n"] += 6
        self.offsets["d"] += 1

    def read_domain_header(self):
        self.offsets['n'] += 2
        self.offsets['i'] += 2
