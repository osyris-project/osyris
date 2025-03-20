# SPDX-License-Identifier: BSD-3-Clause
from .. import config
from ..core import Dataset
from .loader import Loader


class RamsesDataset(Dataset):
    def __init__(self, nout, path=""):
        super().__init__()
        self.loader = Loader(nout=nout, path=path)
        self.meta.update(self.loader.load_metadata())
        self.set_units()
        self.meta["time"] *= self.units["time"]

    def load(self, *args, **kwargs):
        groups = self.loader.load(*args, meta=self.meta, units=self.units, **kwargs)
        for name, group in groups.items():
            self[name] = group
        config.additional_variables(self)
        return self

    def copy(self):
        nout = self.loader.nout
        path = self.loader.path
        out = self.__class__(nout=nout, path=path)
        for key, group in self.groups.items():
            out[key] = group.copy()
        out.meta = self.meta.copy()
        return out
