# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np
import os
from . import utils
from .. import config
from ..core import Datagroup
from .amr import AmrReader
from .grav import GravReader
from .hydro import HydroReader
from .rt import RtReader
from .sink import SinkReader


class Loader:
    def __init__(self, nout, scale, path):
        # Generate directory name from output number
        self.nout = nout
        self.path = path
        self.scale = scale
        self.infile = utils.generate_fname(nout, path)
        self.readers = {
            "amr": AmrReader(),
            "hydro": HydroReader(),
            "grav": GravReader(),
            "rt": RtReader(),
            "sink": SinkReader()
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

        return meta

    def load(self, select=None, cpu_list=None, meta=None):

        out = {}
        groups = list(self.readers.keys())

        if select is None:
            select = {group: {} for group in self.readers}
        else:
            for key in select:
                if key not in self.readers:
                    print("Warning: {} found in select is not a valid "
                          "Datagroup.".format(key))
            for group in self.readers:
                if group not in select:
                    select[group] = {}

        # Take into account user specified lmax
        if "level" in select:
            meta["lmax"] = utils.find_max_amr_level(levelmax=meta["levelmax"],
                                                    select=select["amr"])
        else:
            meta["lmax"] = meta["levelmax"]

        # Initialize readers
        readers = {}
        for group in groups:
            if not self.readers[group].initialized:
                first_load = self.readers[group].initialize(meta=meta,
                                                            select=select[group])
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
                            for group, reader in readers.items():
                                conditions.update(
                                    reader.make_conditions(select[group], ncache))
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

        print("Total number of cells loaded: {}".format(ncells_tot))

        return out
