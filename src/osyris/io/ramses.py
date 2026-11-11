# SPDX-License-Identifier: BSD-3-Clause
from ..core import Dataset
from .loader import Loader

_default_additional_variables = {
    "mesh": {
        "B_field": lambda B_left, B_right: 0.5 * (B_left + B_right),
        "mass": lambda density, dx: density * dx**3,
    }
}


class RamsesDataset(Dataset):
    def __init__(self, nout, path=""):
        super().__init__()
        self.loader = Loader(nout=nout, path=path)
        self.meta.update(self.loader.load_metadata())
        self.set_units()
        self.meta["time"] *= self.units["time"]

    def load(self, *args, additional_variables=None, **kwargs):
        """
        Load the dataset and compute any additional variables.

        Parameters
        ----------
        additional_variables : dict, optional
            A dictionary specifying additional variables to compute.
            Example: to compute the mass inside each cell:
                additional_variables = {
                    "mesh": {
                        "mass": lambda density, dx: density * dx**3
                    }
                }
            where "mesh" refers to the mesh group in the dataset, "mass" is the new
            variable to compute, and "density" and "dx" are both variables found in the
            mesh group.
        """
        groups = self.loader.load(*args, meta=self.meta, units=self.units, **kwargs)
        for name, group in groups.items():
            self[name] = group
        self.compute_additional_variables(additional_variables=additional_variables)
        return self

    def compute_additional_variables(self, additional_variables=None):
        to_add = _default_additional_variables
        if additional_variables is not None:
            for key, value in additional_variables.items():
                if key in to_add and isinstance(to_add[key], dict):
                    to_add[key].update(value)
                else:
                    to_add[key] = value

        for group_name, variables in to_add.items():
            if group_name in self.groups:
                group = self.groups[group_name]
                for var_name, func in variables.items():
                    if var_name not in group:
                        required_vars = func.__code__.co_varnames[
                            : func.__code__.co_argcount
                        ]
                        if all(rv in group for rv in required_vars):
                            group[var_name] = func(*(group[rv] for rv in required_vars))

    def copy(self):
        nout = self.loader.nout
        path = self.loader.path
        out = self.__class__(nout=nout, path=path)
        for key, group in self.groups.items():
            out[key] = group.copy()
        out.meta = self.meta.copy()
        return out
