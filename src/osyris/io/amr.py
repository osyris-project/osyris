import numpy as np
from .loader import Loader
from .units import get_unit
from .. import units
from . import utils


class AmrLoader(Loader):
    def __init__(self, scale=None, code_units=None, ndim=1):

        super().__init__()

        self.initialized = True
        # AMR grid variables
        length_unit = get_unit("x", code_units["ud"], code_units["ul"],
                               code_units["ut"])
        if scale is not None:
            scale = units(scale)
            scaling = (length_unit.to(scale) / scale).magnitude * scale
        else:
            scaling = length_unit

        self.variables.update({
            "level": {
                "read": True,
                "type": "i",
                "buffer": None,
                "pieces": {},
                "unit": 1.0 * units.dimensionless
            },
            "cpu": {
                "read": True,
                "type": "i",
                "buffer": None,
                "pieces": {},
                "unit": 1.0 * units.dimensionless
            },
            "dx": {
                "read": True,
                "type": "d",
                "buffer": None,
                "pieces": {},
                "unit": scaling
            }
        })
        self.variables.update({
            "xyz_{}".format(c): {
                "read": True,
                "type": "d",
                "buffer": None,
                "pieces": {},
                "unit": scaling
            }
            for c in "xyz"[:ndim]
        })

    def allocate_buffers(self, ngridmax, twotondim):
        super().allocate_buffers(ngridmax, twotondim)
        self.xcent = np.zeros([8, 3], dtype=np.float64)
        self.xg = np.zeros([ngridmax, 3], dtype=np.float64)
        self.son = np.zeros([ngridmax, twotondim], dtype=np.int32)
        self.ref = np.zeros([ngridmax, twotondim], dtype=np.bool)

    def read_header(self, info):
        # nx,ny,nz
        self.offsets["i"] += 2
        self.offsets["n"] += 2
        [nx, ny, nz] = utils.read_binary_data(fmt="3i",
                                              content=self.bytes,
                                              offsets=self.offsets)
        ncoarse = nx * ny * nz
        self.meta["xbound"] = [
            float(int(nx / 2)),
            float(int(ny / 2)),
            float(int(nz / 2))
        ]

        # nboundary
        self.offsets["i"] += 2
        self.offsets["n"] += 2
        [self.meta["nboundary"]] = utils.read_binary_data(fmt="i",
                                                          content=self.bytes,
                                                          offsets=self.offsets)
        self.meta["ngridlevel"] = np.zeros(
            [info["ncpu"] + self.meta["nboundary"], info["levelmax"]], dtype=np.int32)

        # noutput
        self.offsets["i"] += 1
        self.offsets["n"] += 2
        self.offsets["d"] += 1
        [noutput] = utils.read_binary_data(fmt="i",
                                           content=self.bytes,
                                           offsets=self.offsets)
        # dtold, dtnew
        self.offsets["i"] += 2
        self.offsets["n"] += 3
        self.offsets["d"] += 1 + 2 * noutput
        info["dtold"] = utils.read_binary_data(fmt="{}d".format(info["levelmax"]),
                                               content=self.bytes,
                                               offsets=self.offsets)
        info["dtnew"] = utils.read_binary_data(fmt="{}d".format(info["levelmax"]),
                                               content=self.bytes,
                                               offsets=self.offsets)

        # Read the number of grids
        self.offsets["i"] += 2 + (2 * info["ncpu"] * info["levelmax"])
        self.offsets["n"] += 7
        self.offsets["d"] += 16
        self.meta["ngridlevel"][:info["ncpu"], :] = np.array(
            utils.read_binary_data(fmt="{}i".format(info["ncpu"] * info["levelmax"]),
                                   content=self.bytes,
                                   offsets=self.offsets)).reshape(
                                       info["levelmax"], info["ncpu"]).T

        # Read boundary grids if any
        self.offsets["i"] += 10 * info["levelmax"]
        self.offsets["n"] += 3
        if self.meta["nboundary"] > 0:
            self.offsets["i"] += (2 * self.meta["nboundary"] * info["levelmax"])
            # self.offsets["n"] += 4
            self.meta["ngridlevel"][info["ncpu"]:info["ncpu"] +
                                    self.meta["nboundary"], :] = np.array(
                                        utils.read_binary_data(
                                            fmt="{}i".format(self.meta["nboundary"] *
                                                             info["levelmax"]),
                                            content=self.bytes,
                                            offsets=self.offsets)).reshape(
                                                info["levelmax"],
                                                self.meta["nboundary"]).T
            self.offsets["n"] += 2

        # Determine bound key precision
        self.offsets["i"] += 5
        self.offsets["s"] += 128
        [key_size] = utils.read_binary_data(fmt="i",
                                            content=self.bytes,
                                            offsets=self.offsets,
                                            skip_head=False,
                                            increment=False)

        # Offset for AMR
        self.offsets["i"] += 3 * ncoarse
        self.offsets["n"] += 3
        self.offsets["s"] += key_size

    def read_level_header(self, ilevel, twotondim):
        # Geometry
        self.dxcell = 0.5**(ilevel + 1)
        for ind in range(twotondim):
            iz = int((ind) / 4)
            iy = int((ind - 4 * iz) / 2)
            ix = int((ind - 2 * iy - 4 * iz))
            self.xcent[ind, 0] = (float(ix) - 0.5) * self.dxcell
            self.xcent[ind, 1] = (float(iy) - 0.5) * self.dxcell
            self.xcent[ind, 2] = (float(iz) - 0.5) * self.dxcell

    def read_cacheline_header(self, ncache, ndim):
        # xg: grid coordinates
        self.offsets['i'] += ncache * 3
        self.offsets['n'] += 3
        for n in range(ndim):
            self.xg[:ncache, n] = utils.read_binary_data(fmt="{}d".format(ncache),
                                                         content=self.bytes,
                                                         offsets=self.offsets)

        # son indices
        self.offsets['i'] += ncache * (1 + 2 * ndim)
        self.offsets['n'] += 1 + 2 * ndim

    def read_variables(self, ncache, ind, ilevel, cpuid, info):

        self.son[:ncache, ind] = utils.read_binary_data(fmt="{}i".format(ncache),
                                                        content=self.bytes,
                                                        offsets=self.offsets)

        self.variables["level"]["buffer"][:ncache, ind] = ilevel + 1
        for n in range(info["ndim"]):
            key = "xyz_" + "xyz"[n]
            self.variables[key]["buffer"][:ncache, ind] = (
                self.xg[:ncache, n] + self.xcent[ind, n] - self.meta["xbound"][n]
            ) * info["boxlen"] * self.variables[key]["unit"].magnitude
        self.variables["dx"]["buffer"][:ncache, ind] = self.dxcell * info[
            "boxlen"] * self.variables["dx"]["unit"].magnitude
        self.variables["cpu"]["buffer"][:ncache, ind] = cpuid + 1

        # Note: use lmax here instead of levelmax because the user might not
        # want to load all levels. levelmax is always the max level in the
        # entire simulation.
        self.ref[:ncache, ind] = np.logical_not(
            np.logical_and(self.son[:ncache, ind] > 0, ilevel < info["lmax"] - 1))

    def make_conditions(self, select, ncache):
        conditions = super().make_conditions(select, ncache)
        conditions.update({"leaf": self.ref[:ncache, :]})
        return conditions

    def read_footer(self, ncache, twotondim):
        # Increment offsets with remainder of the file
        self.offsets['i'] += ncache * 2 * twotondim
        self.offsets['n'] += 2 * twotondim

    def step_over(self, ncache, twotondim, ndim):
        self.offsets['i'] += ncache * (4 + 3 * twotondim + 2 * ndim)
        self.offsets['d'] += ncache * ndim
        self.offsets['n'] += 4 + 3 * twotondim + 3 * ndim
