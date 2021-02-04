

class HydroLoader(Loader):

    def __init__(self, select):

        super().__init__()

        # Read the number of variables from the hydro_file_descriptor.txt
        # and select the ones to be read if specified by user
        hydro = True
        fname_hydro = infile+"/hydro_file_descriptor.txt"
        try:
            descriptor = np.loadtxt(fname_hydro, dtype=str, delimiter=",")
        except IOError:
            hydro = False

        if hydro:
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
                    "unit": get_unit(key, data.meta["unit_d"], data.meta["unit_l"], data.meta["unit_t"])}
        # data.meta["nvar_hydro"] = len(variables_hydro)
        # return hyd
