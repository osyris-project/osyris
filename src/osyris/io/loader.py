from . import utils

class Loader():

    def __init__(self):
        self.variables = {}
        self.offsets = {}
        self.meta = {}
        self.bytes = None
        return

    def allocate_buffers(self, ngridmax, twotondim):
        for item in self.variables.values():
            item["buffer"] = np.zeros([ngridmax,twotondim],
                dtype=np.dtype(item["type"]))


    def read_header(self, *args, **kwargs):
        return

    def read_level_header(self, *args, **kwargs):
        return

    def read_domain_header(self, *args, **kwargs):
        return

    def read_variables(self, ncache, ind):
        for item in self.variables.values():
            if item["read"]:
                item["buffer"][:ncache,ind] = np.array(utils.read_binary_data(
                    fmt="{}{}".format(ncache, item["type"]),
                    content=loaders[group].bytes, offsets=loaders[group].offsets)) * item["unit"].magnitude
            else:
                loaders[group].offsets[item["type"]] += ncache
                loaders[group].offsets["n"] += ncache

