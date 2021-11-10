# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)
import numpy as np
from .hilbert import hilbert_cpu_list
from .reader import Reader, ReaderKind
from .. import config
from .. import units
from . import utils


class AmrReader(Reader):
    def __init__(self):
        super().__init__(kind=ReaderKind.AMR)
        self.cpu_list = None

    def initialize(self, meta, select):
        length_unit = config.get_unit("x", meta["unit_d"], meta["unit_l"],
                                      meta["unit_t"])
        if meta["scale"] is not None:
            scale = units(meta["scale"])
            scaling = (length_unit.to(scale) / scale).magnitude * scale
        else:
            scaling = length_unit

        scaling = utils.get_spatial_scaling(meta["unit_d"], meta["unit_l"],
                                            meta["unit_t"], meta["scale"])

        if select is False:
            meta["lmax"] = 0
            return

        # AMR grid variables
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
            for c in "xyz"[:meta["ndim"]]
        })

        self.cpu_list = hilbert_cpu_list(meta=meta,
                                         scaling=scaling,
                                         select=select,
                                         infofile=meta["infofile"])
        self.initialized = True

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
        info["dtold"] = np.array(
            utils.read_binary_data(fmt="{}d".format(info["levelmax"]),
                                   content=self.bytes,
                                   offsets=self.offsets))
        info["dtnew"] = np.array(
            utils.read_binary_data(fmt="{}d".format(info["levelmax"]),
                                   content=self.bytes,
                                   offsets=self.offsets))

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

        self.variables["level"]["buffer"]._array[:ncache, ind] = ilevel + 1
        for n in range(info["ndim"]):
            key = "xyz_" + "xyz"[n]
            self.variables[key]["buffer"]._array[:ncache, ind] = (
                self.xg[:ncache, n] + self.xcent[ind, n] - self.meta["xbound"][n]
            ) * info["boxlen"] * self.variables[key]["unit"].magnitude
        self.variables["dx"]["buffer"]._array[:ncache, ind] = self.dxcell * info[
            "boxlen"] * self.variables["dx"]["unit"].magnitude
        self.variables["cpu"]["buffer"]._array[:ncache, ind] = cpuid + 1

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
