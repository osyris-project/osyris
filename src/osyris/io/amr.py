import numpy as np
from .loader import Loader
from .units import get_unit
from .. import units
from . import utils


class AmrLoader(Loader):
    def __init__(self, scale, code_units):

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
            "x": {
                "read": True,
                "type": "d",
                "buffer": None,
                "pieces": {},
                "unit": scaling
            },
            "y": {
                "read": True,
                "type": "d",
                "buffer": None,
                "pieces": {},
                "unit": scaling
            },
            "z": {
                "read": True,
                "type": "d",
                "buffer": None,
                "pieces": {},
                "unit": scaling
            },
            "dx": {
                "read": True,
                "type": "d",
                "buffer": None,
                "pieces": {},
                "unit": scaling
            }
        })

    def allocate_buffers(self, ngridmax, twotondim):
        super().allocate_buffers(ngridmax, twotondim)
        xg = np.zeros([ngridmax,3],dtype=np.float64)


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

        # self.meta["nboundary"]
        self.offsets["i"] += 2
        self.offsets["n"] += 2
        [self.meta["nboundary"]] = utils.read_binary_data(fmt="i",
                                                          content=self.bytes,
                                                          offsets=self.offsets)
        self.meta["ngridlevel"] = np.zeros(
            [info["ncpu"] + self.meta["nboundary"], info["levelmax"]],
            dtype=np.int32)
        # print(self.meta["nboundary"])
        # # self.meta["nboundary"]
        # ninteg = 7
        # nlines = 5
        # [self.meta["nboundary"]] = utils.get_binary_data(fmt="i",content=self.bytes,ninteg=ninteg,nlines=nlines)
        # ngridlevel = np.zeros([info["ncpu"]+self.meta["nboundary"],info["levelmax"]],dtype=np.int32)
        # print(self.meta["nboundary"])
        # return

        # noutput
        self.offsets["i"] += 1
        self.offsets["n"] += 2
        self.offsets["d"] += 1
        [noutput] = utils.read_binary_data(fmt="i",
                                           content=self.bytes,
                                           offsets=self.offsets)
        # print(noutput)
        # # noutput
        # ninteg = 9
        # nfloat = 1
        # nlines = 8
        # [noutput] = utils.read_binary_data(fmt="i",content=self.bytes,ninteg=ninteg,nlines=nlines,nfloat=nfloat)

        # dtold, dtnew
        self.offsets["i"] += 2
        self.offsets["n"] += 3
        self.offsets["d"] += 1 + 2 * noutput
        info["dtold"] = utils.read_binary_data(fmt="{}d".format(
            info["levelmax"]),
                                                    content=self.bytes,
                                                    offsets=self.offsets)
        # nfloat += 1+info["levelmax"]
        # nlines += 1
        info["dtnew"] = utils.read_binary_data(fmt="{}d".format(
            info["levelmax"]),
                                                    content=self.bytes,
                                                    offsets=self.offsets)
        # print(info)
        # return
        # # dtold, dtnew
        # ninteg = 12
        # nfloat = 2+2*noutput
        # nlines = 12
        # info["dtold"] = utils.read_binary_data(fmt="%id"%(info["levelmax"]),\
        #                      content=self.bytes,ninteg=ninteg,nlines=nlines,nfloat=nfloat)
        # nfloat += 1+info["levelmax"]
        # nlines += 1
        # info["dtnew"] = utils.read_binary_data(fmt="%id"%(info["levelmax"]),\
        #                      content=self.bytes,ninteg=ninteg,nlines=nlines,nfloat=nfloat)

        # # hydro gamma
        # offsets_hydro["i"] += 5
        # offsets_hydro["n"] += 5
        # [info["gamma"]] = utils.read_binary_data(fmt="d",content=loaders["hydro"].bytes,offsets=offsets_hydro)
        # # hydro gamma
        # ninteg = 5
        # nfloat = 0
        # nlines = 5
        # [info["gamma"]] = utils.read_binary_data(fmt="d",content=loaders["hydro"].bytes,ninteg=ninteg,nlines=nlines)

        # Read the number of grids
        self.offsets["i"] += 2 + (2 * info["ncpu"] *
                                  info["levelmax"])
        self.offsets["n"] += 7
        self.offsets["d"] += 16
        # print(self.offsets["i"], self.offsets["d"], self.offsets["n"])
        self.meta["ngridlevel"][:info["ncpu"], :] = np.array(
            utils.read_binary_data(
                fmt="{}i".format(info["ncpu"] * info["levelmax"]),
                content=self.bytes,
                offsets=self.offsets)).reshape(info["levelmax"],
                                               info["ncpu"]).T
        # print("ngridlevel1", ngridlevel)
        # # # Read the number of grids
        # ninteg = 14+(2*info["ncpu"]*info["levelmax"])
        # nfloat = 18+(2*noutput)+(2*info["levelmax"])
        # nlines = 21
        # print(ninteg, nfloat, nlines)
        # ngridlevel[:info["ncpu"],:] = np.asarray(utils.get_binary_data(fmt="%ii"%(info["ncpu"]*info["levelmax"]),\
        #      content=self.bytes,ninteg=ninteg,nlines=nlines,nfloat=nfloat)).reshape(info["levelmax"],info["ncpu"]).T
        # print("ngridlevel2", ngridlevel)
        # return

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
                                   info["levelmax"], self.meta["nboundary"]).T
            # print("111", ngridlevel[info["ncpu"]:info["ncpu"]+self.meta["nboundary"],:])
            self.offsets["n"] += 2
        # else:
        #     self.offsets["n"] += 3

        # ninteg = 14+(3*info["ncpu"]*info["levelmax"])+(10*info["levelmax"])+(2*self.meta["nboundary"]*info["levelmax"])
        # nfloat = 18+(2*noutput)+(2*info["levelmax"])
        # nlines = 25
        # ngridlevel[info["ncpu"]:info["ncpu"]+self.meta["nboundary"],:] = np.asarray(utils.get_binary_data(fmt="%ii"%(self.meta["nboundary"]*info["levelmax"]),\
        #                                 content=self.bytes,ninteg=ninteg,nlines=nlines,nfloat=nfloat)).reshape(info["levelmax"],self.meta["nboundary"]).T
        # print("222", ngridlevel[info["ncpu"]:info["ncpu"]+self.meta["nboundary"],:])

        # Determine bound key precision
        self.offsets["i"] += 5
        # self.offsets["n"] += 2*min(1,self.meta["nboundary"])
        self.offsets["s"] += 128
        # print(self.offsets["i"], self.offsets["d"], self.offsets["n"], self.offsets["s"])
        [key_size] = utils.read_binary_data(fmt="i",
                                            content=self.bytes,
                                            offsets=self.offsets,
                                            skip_head=False,
                                            increment=False)
        # print("key_size", key_size)
        # # Determine bound key precision
        # ninteg = 14+(3*info["ncpu"]*info["levelmax"])+(10*info["levelmax"])+(3*self.meta["nboundary"]*info["levelmax"])+5
        # nfloat = 18+(2*noutput)+(2*info["levelmax"])
        # nlines = 21+2+3*min(1,self.meta["nboundary"])+1+1
        # nstrin = 128
        # print(ninteg, nfloat, nlines, nstrin)
        # [key_size] = utils.get_binary_data(fmt="i",content=self.bytes,ninteg=ninteg,nlines=nlines,nfloat=nfloat,nstrin=nstrin,correction=-4)
        # print("key_size2", key_size)

        # Offset for AMR
        self.offsets["i"] += 3 * ncoarse
        self.offsets["n"] += 3
        self.offsets["s"] += key_size

    def read_domain_header(self, ncache, ndim):
        # xg: grid coordinates
        self.offsets['i'] += ncache*3
        self.offsets['n'] +=  3
        # ninteg = ninteg_amr + ncache*3
        # nfloat = nfloat_amr
        # nlines = nlines_amr + 3
        # nstrin = nstrin_amr
        for n in range(ndim):
            xg[:ncache,n] = utils.read_binary_data(fmt="{}d".format(ncache),content=self.bytes,offsets=self.offsets)
            # offset = 4*ninteg + 8*(nlines+nfloat+n*(ncache+1)) + nstrin + 4
            # xg[:ncache,n] = struct.unpack("%id"%(ncache), self.bytes[offset:offset+8*ncache])

        # son indices
        self.offsets['i'] += ncache*(1 + 2*ndim)
        self.offsets['n'] += 1 + 2*ndim
        # ninteg = ninteg_amr + ncache*(4+2*data.meta["ndim"])
        # nfloat = nfloat_amr + ncache*data.meta["ndim"]
        # nlines = nlines_amr + 4 + 3*data.meta["ndim"]
        # nstrin = nstrin_amr

    def read_variables(self):
    	loaders["amr"].variables["level"]["buffer"][:ncache,ind] = ilevel + 1
        # scaling = get_unit(
        #             "x", data.meta["unit_d"], data.meta["unit_l"], data.meta["unit_t"]).magnitude
        for n in range(data.meta["ndim"]):
            # xyz[:ncache,ind,n] = xg[:ncache,n] + xcent[ind,n]-xbound[n]
            # var[:ncache,ind,-5+n] = xyz[:ncache,ind,n]*data.meta["boxlen"]
            key = "xyz"[n]
            loaders["amr"].variables[key]["buffer"][:ncache,ind] = (
                xg[:ncache,n] + xcent[ind,n]-loaders["amr"].meta["xbound"][n])*data.meta["boxlen"] * loaders["amr"].variables[key]["unit"].magnitude
        loaders["amr"].variables["dx"]["buffer"][:ncache,ind] = dxcell*data.meta["boxlen"] * loaders["amr"].variables["dx"]["unit"].magnitude
        loaders["amr"].variables["cpu"]["buffer"][:ncache,ind] = cpuid+1
        