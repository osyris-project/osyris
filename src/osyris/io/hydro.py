import numpy as np
from .loader import Loader
from .units import get_unit
from . import utils


class HydroLoader(Loader):
    def __init__(self, infile, select, code_units):

        super().__init__()

        # Read the number of variables from the hydro_file_descriptor.txt
        # and select the ones to be read if specified by user
        self.initialized = True
        fname = infile + "/hydro_file_descriptor.txt"
        try:
            descriptor = np.loadtxt(fname, dtype=str, delimiter=",")
        except IOError:
            self.initialized = False

        if self.initialized:
            for i in range(len(descriptor)):
                key = descriptor[i, 1].strip()
                read = True
                if key in select:
                    if select[key] is False:
                        read = False
                if "hydro" in select:
                    if select["hydro"] is False:
                        read = False
                self.variables[key] = {
                    "read":
                    read,
                    "type":
                    descriptor[i, 2].strip(),
                    "buffer":
                    None,
                    "pieces": {},
                    "unit":
                    get_unit(key, code_units["ud"], code_units["ul"],
                             code_units["ut"])
                }

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
