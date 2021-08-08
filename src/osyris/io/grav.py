import os
from .loader import Loader
from .units import get_unit
from . import utils


class GravLoader(Loader):
    def __init__(self, nout, path, select, units, ndim):

        super().__init__()

        # Check if self-gravity files exist
        fname = utils.generate_fname(nout, path, ftype="grav", cpuid=1)
        self.initialized = os.path.exists(fname)
        # Add gravity fields
        if self.initialized:
            descriptor = {"grav_potential": "d"}
            for n in range(ndim):
                descriptor["grav_acceleration_" + "xyz"[n]] = "d"
            # Now add to the list of variables to be read
            for key in descriptor:
                read = True
                if key in select:
                    if select[key] is False:
                        read = False
                if "gravity" in select:
                    if select["gravity"] is False:
                        read = False
                if "grav" in select:
                    if select["grav"] is False:
                        read = False
                self.variables[key] = {
                    "read": read,
                    "type": descriptor[key],
                    "buffer": None,
                    "pieces": {},
                    "unit": get_unit(key, units["ud"], units["ul"], units["ut"])
                }

    def read_header(self, info):
        self.offsets["i"] += 4
        self.offsets["n"] += 4

    def read_domain_header(self):
        self.offsets['n'] += 2
        self.offsets['i'] += 2
