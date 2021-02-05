from .loader import Loader
from .units import get_unit

class RtLoader(Loader):

    def __init__(self, infile, select, units):

        super().__init__()

        # Read the number of variables from the hydro_file_descriptor.txt
        # and select the ones to be read if specified by user
        self.initialized = True
        fname = infile+"/hydro_file_descriptor.txt"
        try:
            descriptor = np.loadtxt(fname, dtype=str, delimiter=",")
        except IOError:
            self.initialized = False

        if self.initialized:
            # loader["hydro"] = {"variables": {}, "offsets": {}, "bytes": {}}
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
                    "read": read,
                    "type": descriptor[i, 2].strip(),
                    "buffer": None,
                    "pieces": {},
                    "unit": get_unit(key, units["ud"], units["ul"], units["ut"])}
        # data.meta["nvar_hydro"] = len(variables_hydro)
        # return hyd
