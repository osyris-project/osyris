import numpy as np
from .loader import Loader
from .units import get_unit


class RtLoader(Loader):
    def __init__(self, infile, select, units):

        super().__init__()

        # Read the number of variables from the rt_file_descriptor.txt
        # and select the ones to be read if specified by user
        self.initialized = True
        fname = infile + "/rt_file_descriptor.txt"
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
                if "rt" in select:
                    if select["rt"] is False:
                        read = False
                self.variables[key] = {
                    "read": read,
                    "type": descriptor[i, 2].strip(),
                    "buffer": None,
                    "pieces": {},
                    "unit": get_unit(key, units["ud"], units["ul"], units["ut"])
                }

    def read_header(self, info):
        self.offsets["i"] += 5
        self.offsets["n"] += 6
        self.offsets["d"] += 1

    def read_domain_header(self):
        self.offsets['n'] += 2
        self.offsets['i'] += 2
