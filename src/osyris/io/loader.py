# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np
import os
from .. import config  # import parameters, additional_variables
from . import utils
from .. import config
from ..core import Datagroup
from .amr import AmrReader
from .grav import GravReader
from .hydro import HydroReader
from .rt import RtReader
from .sinks import SinkReader


class Loader:
    def __init__(self, nout, scale, path):

        # self.data = Dataset()

        # if scale is None:
        #     scale = config.parameters["scale"]

        # Generate directory name from output number
        self.nout = nout
        self.path = path
        self.scale = scale
        self.infile = utils.generate_fname(nout, path)
        self.readers = {
            "amr": AmrReader(),
            "hydro": HydroReader(),
            # "grav": GravReader(),
            # "rt": RtReader()
            "sinks": SinkReader()
        }

    def load_metadata(self):

        # Read info file and create info dictionary
        infofile = os.path.join(self.infile,
                                "info_" + self.infile.split("_")[-1] + ".txt")
        meta = utils.read_parameter_file(fname=infofile)
        # Add additional information
        meta["infofile"] = infofile
        meta["scale"] = self.scale
        meta["infile"] = self.infile
        meta["nout"] = self.nout
        meta["path"] = self.path
        meta["time"] *= config.get_unit("time", meta["unit_d"], meta["unit_l"],
                                        meta["unit_t"])

        # Read namelist.txt file
        meta["namelist"] = {}
        namelist_file = os.path.join(self.infile, "namelist.txt")
        with open(namelist_file, 'r') as f:
            content = f.read()
        groups = content.split('&')
        for group in groups:
            if "=" in group:
                split = group.split('\n')
                trimmed = [ele for ele in split if ele.strip()]
                if trimmed[-1] == '/':
                    key = trimmed[0].strip()
                    meta["namelist"][key] = {}
                    for item in trimmed[1:-1]:
                        [left, right] = item.split('=')
                        try:
                            right = eval(right)
                        except (NameError, SyntaxError):
                            pass
                        meta["namelist"][key][left] = right

        return meta

    #     # # Take into account user specified lmax
    #     # if "level" in select:
    #     #     meta["lmax"] = utils.find_max_amr_level(
    #     #         levelmax=meta["levelmax"], select=select)
    #     # else:
    #     #     meta["lmax"] = meta["levelmax"]

    #     code_units = {"ud": meta["unit_d"], "ul": meta["unit_l"], "ut": meta["unit_t"]}

    #     self.reader_list = {
    #         "amr":
    #         AmrReader(
    #             scale=scale,
    #             # select=select,
    #             code_units=code_units,
    #             meta=meta,
    #             infofile=infofile),
    #         "hydro":
    #         HydroReader(infile=infile, code_units=code_units),
    #         "grav":
    #         GravReader(
    #             nout=nout,
    #             path=path,
    #             # select=select,
    #             code_units=code_units,
    #             ndim=meta["ndim"]),
    #         "rt":
    #         RtReader(infile=infile, code_units=code_units)
    #     }

    #     self.sinks = read_sinks(nout=nout, path=path, code_units=code_units)

    # # @property
    # # def meta(self):
    # #     return meta

    def load(self, groups=None, select=None, cpu_list=None, meta=None):

        out = {}
        meta = meta.copy()

        if groups is None:
            groups = list(self.readers.keys())
        if "amr" not in groups:
            groups.append("amr")

        if select is None:
            select = {}

        # Take into account user specified lmax
        if "level" in select:
            meta["lmax"] = utils.find_max_amr_level(levelmax=meta["levelmax"],
                                                    select=select)
        else:
            meta["lmax"] = meta["levelmax"]

        # Initialize readers
        readers = {}
        # for group, reader in self.reader_list.items():
        for group in groups:
            if not self.readers[group].initialized:
                first_load = self.readers[group].initialize(meta=meta, select=select)
                if first_load is not None:
                    out[group] = first_load
            if self.readers[group].initialized:
                readers[group] = self.readers[group]

        # Take into account user specified cpu list
        if cpu_list is None:
            cpu_list = self.readers["amr"].cpu_list if self.readers[
                "amr"].cpu_list is not None else range(1, meta["ncpu"] + 1)

        print("Processing {} files in {}".format(len(cpu_list), meta["infile"]))

        # Allocate work arrays
        twotondim = 2**meta["ndim"]
        for reader in readers.values():
            reader.allocate_buffers(ngridmax=meta["ngridmax"], twotondim=twotondim)

        iprog = 1
        istep = 10
        ncells_tot = 0
        npieces = 0

        # integer, double, line, string, quad, long
        null_offsets = {key: 0 for key in "idnsql"}

        # Loop over the cpus and read the AMR and HYDRO files in binary format
        for cpu_ind, cpu_num in enumerate(cpu_list):

            # Print progress
            percentage = int(float(cpu_ind) * 100.0 / float(len(cpu_list)))
            if percentage >= iprog * istep:
                print("{:>3d}% : read {:>10d} cells".format(percentage, ncells_tot))
                iprog += 1

            # Read binary files
            for group, reader in readers.items():
                fname = utils.generate_fname(meta["nout"],
                                             meta["path"],
                                             ftype=group,
                                             cpuid=cpu_num)
                with open(fname, mode='rb') as f:
                    reader.bytes = f.read()

            # Read file headers
            for reader in readers.values():
                reader.offsets.update(null_offsets)
                reader.read_header(meta)

            # Loop over levels
            for ilevel in range(meta["lmax"]):

                for reader in readers.values():
                    reader.read_level_header(ilevel, twotondim)

                # Loop over domains
                for domain in range(readers["amr"].meta["nboundary"] + meta["ncpu"]):

                    ncache = readers["amr"].meta["ngridlevel"][domain, ilevel]

                    for reader in readers.values():
                        reader.read_domain_header()

                    if ncache > 0:

                        if domain == cpu_num - 1:

                            for reader in readers.values():
                                reader.read_cacheline_header(ncache, meta["ndim"])

                            for ind in range(twotondim):

                                # Read variables in cells
                                for reader in readers.values():
                                    reader.read_variables(ncache, ind, ilevel,
                                                          cpu_num - 1, meta)

                            # Apply selection criteria: select only leaf cells and
                            # add any criteria requested by the user via select.
                            conditions = {}
                            for reader in readers.values():
                                conditions.update(reader.make_conditions(
                                    select, ncache))
                            # Combine all selection criteria together with AND
                            # operation by using a product on bools
                            sel = np.where(
                                np.prod(np.array(list(conditions.values())), axis=0))

                            # Count the number of cells
                            ncells = np.shape(sel)[1]
                            if ncells > 0:
                                ncells_tot += ncells
                                npieces += 1
                                # Add the cells in the pieces dictionaries
                                for reader in readers.values():
                                    for item in reader.variables.values():
                                        if item["read"]:
                                            item["pieces"][npieces] = item["buffer"][
                                                sel]

                            # Increment offsets with remainder of the file
                            for reader in readers.values():
                                reader.read_footer(ncache, twotondim)

                        else:

                            for reader in readers.values():
                                reader.step_over(ncache, twotondim, meta["ndim"])

        # Store the number of cells
        meta["ncells"] = ncells_tot

        # Merge all the data pieces into the Arrays
        for group, reader in readers.items():
            out[group] = Datagroup()
            for key, item in reader.variables.items():
                if item["read"]:
                    out[group][key] = np.concatenate(list(item["pieces"].values()))
            # If vector quantities are found, make them into vector Arrays
            utils.make_vector_arrays(out[group], ndim=meta["ndim"])

        # # Create additional variables derived from the ones already loaded
        # config.additional_variables(self.data)

        print("Total number of cells loaded: {}".format(ncells_tot))
        # print("Memory used: {}".format(self.data.print_size()))

        return out
