# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)
import numpy as np
import os
from .. import config
from ..io import Loader
from .datagroup import Datagroup
from .tools import bytes_to_human_readable
# from ..utils import value_to_string


class Dataset:
    def __init__(self, nout=1, scale=None, path=""):
        self.groups = {}
        self.meta = {}

        if scale is None:
            scale = config.parameters["scale"]

        self.loader = Loader(nout=nout, scale=scale, path=path)

        self.meta.update(self.loader.load_meta_info())

        # code_units = {
        #     "ud": self.meta["unit_d"],
        #     "ul": self.meta["unit_l"],
        #     "ut": self.meta["unit_t"]
        # }

        # self.readers = {
        #     "amr": AmrReader(),
        #     "hydro": HydroReader(),
        #     "grav": GravReader(),
        #     "rt": RtReader()
        # }
        # "sinks":
        # SinksReader(nout=nout, path=path, code_units=code_units)
        # }

        # # Generate directory name from output number
        # infile = utils.generate_fname(nout, path)

        # # Read info file and create info dictionary
        # infofile = os.path.join(infile, "info_" + infile.split("_")[-1] + ".txt")
        # self.meta.update(utils.read_parameter_file(fname=infofile))

        # # Add additional information
        # self.meta["scale"] = scale
        # self.meta["infile"] = infile
        # self.meta["nout"] = nout
        # self.meta["path"] = path
        # self.meta["time"] *= config.get_unit("time", self.meta["unit_d"],
        #                                      self.meta["unit_l"], self.meta["unit_t"])

        # self.readers = {
        #     "amr":
        #     AmrReader(
        #         scale=scale,
        #         # select=select,
        #         code_units=code_units,
        #         meta=self.data.meta,
        #         infofile=infofile),
        #     "hydro":
        #     HydroReader(infile=infile, code_units=code_units),
        #     "grav":
        #     GravReader(
        #         nout=nout,
        #         path=path,
        #         # select=select,
        #         code_units=code_units,
        #         ndim=self.data.meta["ndim"]),
        #     "rt":
        #     RtReader(infile=infile, code_units=code_units),
        #     # "sinks":
        #     # SinksReader(nout=nout, path=path, code_units=code_units)
        # }

    def __iter__(self):
        return self.groups.__iter__()

    def __len__(self):
        return self.groups.__len__()

    def __getitem__(self, key):
        return self.groups[key]
        # if isinstance(key, str):
        #     return self.groups[key]
        # else:
        #     d = Dataset()
        #     for name, val in self.items():
        #         d[name] = val[key]
        #     d.meta.update(self.meta)
        #     return d

    def __setitem__(self, key, value):
        if not isinstance(value, Datagroup):
            raise TypeError("Only objects of type Datagroup can be inserted into a "
                            "Dataset.")
        self.groups.__setitem__(key, value)
        # if len(value.shape) > 0:
        #     shape = value.shape[0]
        # else:
        #     shape = 1
        # if self.shape is None:
        #     self.shape = shape
        # else:
        #     if self.shape != shape:
        #         raise RuntimeError(
        #             "Size mismatch on element insertion. Item "
        #             "shape is {} while container accepts shape {}.".format(
        #                 shape, self.shape))
        # if isinstance(value, Array):
        value.name = key
        value.parent = self
        # print(key, parent)
        #     self.groups[key] = value
        # else:
        #     self.groups[key] = Array(values=value, name=key, parent=self)

    def __delitem__(self, key):
        return self.groups.__delitem__(key)

    def __repr__(self):
        return str(self)

    def __str__(self):
        output = "Dataset: "
        if "infile" in self.meta:
            output += "{}: ".format(self.meta["infile"])
        output += "{}\n".format(self.print_size())
        for key, item in self.items():
            output += str(item) + "\n"
        return output

    def keys(self):
        return self.groups.keys()

    def items(self):
        return self.groups.items()

    def values(self):
        return self.groups.values()

    def load(self, *args, **kwargs):
        groups = self.loader.load(*args, meta=self.meta, **kwargs)
        for name, group in groups.items():
            self[name] = group
        config.additional_variables(self)
        return self

    # def set_scale(self, scale):
    #     for key in ["x", "y", "z", "dx"]:
    #         if key in self:
    #             self[key].to(scale)

    def nbytes(self):
        return np.sum([item.nbytes() for item in self.groups.values()])

    def print_size(self):
        return bytes_to_human_readable(self.nbytes())

    # def print_size(self):
    #     total_size = np.sum([item._array.nbytes for item in self.values()])
    #     multipliers = {"G": 1.0e9, "M": 1.0e6, "K": 1.0e3, "B": 1.0}
    #     size = "0B"
    #     for m, mult in multipliers.items():
    #         if total_size >= mult:
    #             size = "{:.2f} {}B".format(total_size / mult, m)
    #             break
    #     return size
