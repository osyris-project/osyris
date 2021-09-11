# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np
import os
from ..config import parameters, additional_variables
from . import utils
from ..core import Dataset
from .amr import AmrReader
from .grav import GravReader
from .hydro import HydroReader
from .rt import RtReader
from .sinks import read_sinks
from .units import get_unit


class Loader:
    def __init__(self, nout=1, scale=None, path=""):

        self.data = Dataset()

        if scale is None:
            scale = parameters["scale"]

        # Generate directory name from output number
        infile = utils.generate_fname(nout, path)

        # Read info file and create info dictionary
        infofile = os.path.join(infile, "info_" + infile.split("_")[-1] + ".txt")
        self.data.meta.update(utils.read_parameter_file(fname=infofile))

        # Add additional information
        self.data.meta["scale"] = scale
        self.data.meta["infile"] = infile
        self.data.meta["nout"] = nout
        self.data.meta["path"] = path
        self.data.meta["time"] *= get_unit("time", self.data.meta["unit_d"],
                                           self.data.meta["unit_l"],
                                           self.data.meta["unit_t"])

        # # Take into account user specified lmax
        # if "level" in select:
        #     self.data.meta["lmax"] = utils.find_max_amr_level(
        #         levelmax=self.data.meta["levelmax"], select=select)
        # else:
        #     self.data.meta["lmax"] = self.data.meta["levelmax"]

        code_units = {
            "ud": self.data.meta["unit_d"],
            "ul": self.data.meta["unit_l"],
            "ut": self.data.meta["unit_t"]
        }

        self.reader_list = {
            "amr":
            AmrReader(
                scale=scale,
                # select=select,
                code_units=code_units,
                meta=self.data.meta,
                infofile=infofile),
            "hydro":
            HydroReader(infile=infile, code_units=code_units),
            "grav":
            GravReader(
                nout=nout,
                path=path,
                # select=select,
                code_units=code_units,
                ndim=self.data.meta["ndim"]),
            "rt":
            RtReader(infile=infile, code_units=code_units)
        }

        self.sinks = read_sinks(nout, path)

    @property
    def meta(self):
        return self.data.meta

    def load(self, select=None, cpu_list=None):

        if select is None:
            select = {}

        # Take into account user specified lmax
        if "level" in select:
            self.data.meta["lmax"] = utils.find_max_amr_level(
                levelmax=self.data.meta["levelmax"], select=select)
        else:
            self.data.meta["lmax"] = self.data.meta["levelmax"]

        # Initialize readers
        readers = {}
        for group, reader in self.reader_list.items():
            if reader.initialize(select):
                readers[group] = reader

        # Take into account user specified cpu list
        if cpu_list is None:
            cpu_list = self.reader_list["amr"].cpu_list if self.reader_list[
                "amr"].cpu_list is not None else range(1, self.data.meta["ncpu"] + 1)

        print("Processing {} files in {}".format(len(cpu_list),
                                                 self.data.meta["infile"]))

        # Allocate work arrays
        twotondim = 2**self.data.meta["ndim"]
        for reader in readers.values():
            reader.allocate_buffers(ngridmax=self.data.meta["ngridmax"],
                                    twotondim=twotondim)

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
                fname = utils.generate_fname(self.data.meta["nout"],
                                             self.data.meta["path"],
                                             ftype=group,
                                             cpuid=cpu_num)
                with open(fname, mode='rb') as f:
                    reader.bytes = f.read()

            # Read file headers
            for reader in readers.values():
                reader.offsets.update(null_offsets)
                reader.read_header(self.data.meta)

            # Loop over levels
            for ilevel in range(self.data.meta["lmax"]):

                for reader in readers.values():
                    reader.read_level_header(ilevel, twotondim)

                # Loop over domains
                for domain in range(readers["amr"].meta["nboundary"] +
                                    self.data.meta["ncpu"]):

                    ncache = readers["amr"].meta["ngridlevel"][domain, ilevel]

                    for reader in readers.values():
                        reader.read_domain_header()

                    if ncache > 0:

                        if domain == cpu_num - 1:

                            for reader in readers.values():
                                reader.read_cacheline_header(ncache,
                                                             self.data.meta["ndim"])

                            for ind in range(twotondim):

                                # Read variables in cells
                                for reader in readers.values():
                                    reader.read_variables(ncache, ind, ilevel,
                                                          cpu_num - 1, self.data.meta)

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
                                reader.step_over(ncache, twotondim,
                                                 self.data.meta["ndim"])

        # Store the number of cells
        self.data.meta["ncells"] = ncells_tot

        # Merge all the data pieces into the Arrays
        for group in readers.values():
            for key, item in group.variables.items():
                if item["read"]:
                    self.data[key] = np.concatenate(list(item["pieces"].values()))

        # If vector quantities are found, make them into vector Arrays
        utils.make_vector_arrays(self.data)

        # Create additional variables derived from the ones already loaded
        additional_variables(self.data)

        print("Total number of cells loaded: {}".format(ncells_tot))
        print("Memory used: {}".format(self.data.print_size()))

        return self.data
