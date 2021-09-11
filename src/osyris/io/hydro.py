import numpy as np
from .reader import Reader
from .units import get_unit
from . import utils


class HydroReader(Reader):
    def __init__(self, infile, code_units):
        super().__init__(code_units=code_units)
        # self.infile = infile
        self.fname = infile + "/hydro_file_descriptor.txt"
        # self.code_units = code_units

    def initialize(self, select):
        # Read the number of variables from the hydro_file_descriptor.txt
        # and select the ones to be read if specified by user
        # self.initialized = True
        # fname = self.infile + "/hydro_file_descriptor.txt"
        try:
            descriptor = np.loadtxt(self.fname, dtype=str, delimiter=",")
        except IOError:
            return False

        # if self.initialized:
        for i in range(len(descriptor)):
            key = descriptor[i, 1].strip()
            read = True
            if "hydro" in select:
                if select["hydro"] is False:
                    read = False
            if key in select:
                if isinstance(select[key], bool):
                    read = select[key]
            self.variables[key] = {
                "read":
                read,
                "type":
                descriptor[i, 2].strip(),
                "buffer":
                None,
                "pieces": {},
                "unit":
                get_unit(key, self.code_units["ud"], self.code_units["ul"],
                         self.code_units["ut"])
            }
        return True

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
