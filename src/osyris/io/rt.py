# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2024 Osyris contributors (https://github.com/osyris-project/osyris)
import numpy as np
import os
from .reader import Reader, ReaderKind


class RtReader(Reader):

    def __init__(self):
        super().__init__(kind=ReaderKind.AMR)

    def initialize(self, meta, units, select):
        self.initialized = False
        if select is False:
            return

        # Read the number of variables from the rt_file_descriptor.txt
        # and select the ones to be read if specified by user
        fname = os.path.join(meta["infile"], "rt_file_descriptor.txt")
        try:
            desc_from_file = np.loadtxt(fname, dtype=str, delimiter=",")
        except IOError:
            return

        descriptor = {
            desc_from_file[i, 1].strip(): desc_from_file[i, 2].strip()
            for i in range(len(desc_from_file))
        }

        self.descriptor_to_variables(descriptor=descriptor,
                                     meta=meta,
                                     units=units,
                                     select=select)
        self.initialized = True

    def read_header(self, info):
        self.offsets["i"] += 5
        self.offsets["n"] += 6
        self.offsets["d"] += 1

    def read_domain_header(self):
        self.offsets['n'] += 2
        self.offsets['i'] += 2
