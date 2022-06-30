# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)

import numpy as np
import os
from . import utils
from ..core import Datagroup, Array
from .amr import AmrReader
from .grav import GravReader
from .hydro import HydroReader
from .part import PartReader
from .rt import RtReader
from .sink import SinkReader
from .reader import ReaderKind


class Loader:

    def __init__(self, nout, path):
        # Generate directory name from output number
        self.nout = nout
        self.path = path
        self.infile = utils.generate_fname(nout, path)
        self.readers = {
            "amr": AmrReader(),
            "hydro": HydroReader(),
            "grav": GravReader(),
            "part": PartReader(),
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
        meta["infile"] = self.infile
        meta["nout"] = self.nout
        meta["path"] = self.path
        meta["ncells"] = 0
        meta["nparticles"] = 0
        return meta

    def load(self, select=None, cpu_list=None, sortby=None, meta=None, units=None):

        out = {}
        groups = list(self.readers.keys())

        _select = {group: {} for group in self.readers}
        if isinstance(select, dict):
            for key in select:
                if key not in self.readers:
                    print("Warning: {} found in select is not a valid "
                          "Datagroup.".format(key))
                else:
                    _select[key] = select[key]
        elif isinstance(select, str):
            for key in _select:
                if key != select:
                    _select[key] = False
        elif isinstance(select, list) or isinstance(select, tuple):
            for key in _select:
                if key not in select:
                    _select[key] = False
        # Replace aliases for x,y,z in select: x,y,x -> position_x,y,z
        for group in _select.values():
            if isinstance(group, dict):
                for c in "xyz":
                    if c in group:
                        group[f"position_{c}"] = group.pop(c)

        # Take into account user specified lmax
        meta["lmax"] = meta["levelmax"]
        if "amr" in _select:
            if _select["amr"]:
                if "level" in _select["amr"]:
                    meta["lmax"] = utils.find_max_amr_level(levelmax=meta["levelmax"],
                                                            select=_select["amr"])

        # Initialize readers
        readers = {}
        do_not_load_amr = True
        do_not_load_cpus = True
        for group in groups:
            loaded_on_init = self.readers[group].initialize(meta=meta,
                                                            units=units,
                                                            select=_select[group])
            if loaded_on_init is not None:
                out[group] = loaded_on_init
            if self.readers[group].initialized:
                readers[group] = self.readers[group]
                if self.readers[group].kind == ReaderKind.AMR:
                    do_not_load_amr = False
                if self.readers[group].kind in (ReaderKind.AMR, ReaderKind.PART):
                    do_not_load_cpus = False
        # If no reader requires the AMR tree to be read, set lmax to zero
        if do_not_load_amr:
            lmax = 0
        else:
            meta["ncells"] = 0
            lmax = meta["lmax"]
            if "amr" not in readers:
                readers["amr"] = self.readers["amr"]

        # Take into account user specified cpu list
        if cpu_list is None:
            cpu_list = self.readers["amr"].cpu_list if self.readers[
                "amr"].cpu_list is not None else range(1, meta["ncpu"] + 1)

        # If no reader requires the CPUs (if loading only sinks), make empty cpu list
        if do_not_load_cpus:
            cpu_list = []
        else:
            meta["nparticles"] = 0
            print("Processing {} files in {}".format(len(cpu_list), meta["infile"]))

        # Allocate work arrays
        twotondim = 2**meta["ndim"]

        iprog = 1
        istep = 10
        npieces = 0

        # byte, integer, double, line, string, quad, long
        null_offsets = {key: 0 for key in "bidnsql"}

        # Loop over the cpus and read the AMR and HYDRO files in binary format
        for cpu_ind, cpu_num in enumerate(cpu_list):

            # Print progress
            percentage = int(float(cpu_ind) * 100.0 / float(len(cpu_list)))
            if percentage >= iprog * istep:
                print("{:>3d}% : read {:>10d} cells, {:>10d} particles".format(
                    percentage, meta["ncells"], meta["nparticles"]))
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
            for ilevel in range(lmax):

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
                                reader.allocate_buffers(ncache, twotondim)

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
                                conditions.update(reader.make_conditions(
                                    _select[group]))
                            # Combine all selection criteria together with AND
                            # operation by using a product on bools
                            sel = np.prod(np.array([
                                c.values if isinstance(c, Array) else c
                                for c in conditions.values()
                            ]),
                                          axis=0).astype(bool)

                            # Count the number of cells
                            ncells = np.sum(sel)
                            if ncells > 0:
                                meta["ncells"] += ncells
                                npieces += 1
                                # Add the cells in the pieces dictionaries
                                for reader in readers.values():
                                    if reader.kind == ReaderKind.AMR:
                                        for item in reader.variables.values():
                                            if item["read"]:
                                                item["pieces"][npieces] = item[
                                                    "buffer"][sel]

                            # Increment offsets with remainder of the file
                            for reader in readers.values():
                                reader.read_footer(ncache, twotondim)

                        else:

                            for reader in readers.values():
                                reader.step_over(ncache, twotondim, meta["ndim"])

        # Merge all the data pieces into the Arrays
        for group, reader in readers.items():
            out[group] = Datagroup()
            for key, item in reader.variables.items():
                if item["read"] and len(item["pieces"]) > 0:
                    out[group][key] = np.concatenate(list(item["pieces"].values()))
            # If vector quantities are found, make them into vector Arrays
            utils.make_vector_arrays(out[group], ndim=meta["ndim"])

        print("Loaded: {} cells, {} particles.".format(meta["ncells"],
                                                       meta["nparticles"]))

        # Apply sorting if any requested from args
        if sortby is not None:
            for group, key in sortby.items():
                if group in out:
                    out[group].sortby(key)

        return out
