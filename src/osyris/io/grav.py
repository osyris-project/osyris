from .loader import Loader
from .units import get_unit
from . import utils

class GravLoader(Loader):

    def __init__(self, nout, path, select, units):

        super().__init__()

        # # Check if self-gravity files exist
        fname = utils.generate_fname(nout,path,ftype="grav",cpuid=1)
        gravity = os.path.exists(fname_grav)
        # Add gravity fields
        if gravity:
            loaders["grav"] = {"variables": {}, "offsets": {}, "bytes": {}}
            descriptor = {"grav_potential": "d"}
            for n in range(data.meta["ndim"]):
                descriptor["grav_acceleration_" + "xyz"[n]]: "d"
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
                loaders["grav"]["variables"][key] = {
                    "read": read,
                    "type": descriptor[key],
                    "buffer": None,
                    "pieces": {},
                    "unit": get_unit(key, units["ud"], units["ul"], units["ut"])}


        # # Read the number of variables from the hydro_file_descriptor.txt
        # # and select the ones to be read if specified by user
        # self.initialized = True
        # fname = infile+"/hydro_file_descriptor.txt"
        # try:
        #     descriptor = np.loadtxt(fname, dtype=str, delimiter=",")
        # except IOError:
        #     self.initialized = False

        # if self.initialized:
        #     # loader["hydro"] = {"variables": {}, "offsets": {}, "bytes": {}}
        #     for i in range(len(descriptor)):
        #         key = descriptor[i, 1].strip()
        #         read = True
        #         if key in select:
        #             if select[key] is False:
        #                 read = False
        #         if "hydro" in select:
        #             if select["hydro"] is False:
        #                 read = False
        #         self.variables[key] = {
        #             "read": read,
        #             "type": descriptor[i, 2].strip(),
        #             "buffer": None,
        #             "pieces": {},
        #             "unit": get_unit(key, units["ud"], units["ul"], units["ut"])}
        # # data.meta["nvar_hydro"] = len(variables_hydro)
        # # return hyd
