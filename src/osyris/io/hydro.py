# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/nvaytet/osyris)
import numpy as np
import os
from .reader import Reader, ReaderKind
from .. import config
from . import utils


class HydroReader(Reader):
    def __init__(self):
        super().__init__(kind=ReaderKind.AMR)

    def initialize(self, meta, select):
        self.initialized = False
        if select is False:
            return

        # Read the number of variables from the hydro_file_descriptor.txt
        # and select the ones to be read if specified by user
        fname = os.path.join(meta["infile"], "hydro_file_descriptor.txt")
        try:
            desc_from_file = np.loadtxt(fname, dtype=str, delimiter=",")
        except IOError:
            return

        descriptor = {
            desc_from_file[i, 1].strip(): desc_from_file[i, 2].strip()
            for i in range(len(desc_from_file))
        }

        self.descriptor_to_variables(descriptor=descriptor, meta=meta, select=select)

        # drop_others = False
        # if isinstance(select, dict):
        #     for key, value in select.items():
        #         if value is True:
        #             drop_others = True

        # for i in range(len(descriptor)):
        #     key = descriptor[i, 1].strip()
        #     read = True
        #     if isinstance(select, bool):
        #         read = select
        #     elif key in select:
        #         if isinstance(select[key], bool):
        #             read = select[key]
        #     elif drop_others:
        #         read = False
        #     self.variables[key] = {
        #         "read":
        #         read,
        #         "type":
        #         descriptor[i, 2].strip(),
        #         "buffer":
        #         None,
        #         "pieces": {},
        #         "unit":
        #         config.get_unit(key, meta["unit_d"], meta["unit_l"], meta["unit_t"],
        #                         meta["scale"])
        #     }
        self.initialized = True

    def read_header(self, info):
        # hydro gamma
        self.offsets["i"] += 5
        self.offsets["n"] += 5
        [info["gamma"]] = utils.read_binary_data(fmt="d",
                                                 content=self.bytes,
                                                 offsets=self.offsets)

    def read_domain_header(self):
        self.offsets['n'] += 2
        self.offsets['i'] += 2
