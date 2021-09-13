import numpy as np
from .reader import Reader
from .. import config


class RtReader(Reader):
    def __init__(self, infile, code_units):
        super().__init__(code_units=code_units)
        self.fname = infile + "/rt_file_descriptor.txt"

    def initialize(self, select):
        # Read the number of variables from the rt_file_descriptor.txt
        # and select the ones to be read if specified by user
        # fname = infile + "/rt_file_descriptor.txt"
        try:
            descriptor = np.loadtxt(self.fname, dtype=str, delimiter=",")
        except IOError:
            return False

        for i in range(len(descriptor)):
            key = descriptor[i, 1].strip()
            read = True
            if "rt" in select:
                if select["rt"] is False:
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
                config.get_unit(key, self.code_units["ud"], self.code_units["ul"],
                                self.code_units["ut"])
            }
        return True

    def read_header(self, info):
        self.offsets["i"] += 5
        self.offsets["n"] += 6
        self.offsets["d"] += 1

    def read_domain_header(self):
        self.offsets['n'] += 2
        self.offsets['i'] += 2
