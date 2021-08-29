# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np
from ..config import parameters, additional_variables
from . import utils
from ..core import Dataset
from .amr import AmrLoader
from .grav import GravLoader
from .hydro import HydroLoader
from .rt import RtLoader
from .units import get_unit


def load(nout=1, scale=None, path="", select=None, cpu_list=None, bounding_box=None):

    data = Dataset()

    if select is None:
        select = {}

    if scale is None:
        scale = parameters["scale"]

    # Generate directory name from output number
    infile = utils.generate_fname(nout, path)

    # Read info file and create info dictionary
    infofile = infile + "/info_" + infile.split("_")[-1] + ".txt"
    data.meta.update(utils.read_parameter_file(fname=infofile))

    # Add additional information
    data.meta["scale"] = scale
    data.meta["infile"] = infile
    data.meta["path"] = path
    data.meta["time"] = data.meta["time"] * get_unit(
        "time", data.meta["unit_d"], data.meta["unit_l"], data.meta["unit_t"])

    # Take into account user specified lmax
    if "level" in select:
        data.meta["lmax"] = utils.find_max_amr_level(levelmax=data.meta["levelmax"],
                                                     select=select)
    else:
        data.meta["lmax"] = data.meta["levelmax"]

    code_units = {
        "ud": data.meta["unit_d"],
        "ul": data.meta["unit_l"],
        "ut": data.meta["unit_t"]
    }

    loader_list = {
        "amr":
        AmrLoader(scale=scale,
                  select=select,
                  code_units=code_units,
                  meta=data.meta,
                  infofile=infofile),
        "hydro":
        HydroLoader(infile=infile, select=select, code_units=code_units),
        "grav":
        GravLoader(nout=nout,
                   path=path,
                   select=select,
                   units=code_units,
                   ndim=data.meta["ndim"]),
        "rt":
        RtLoader(infile=infile, select=select, units=code_units)
    }

    # Initialize loaders
    loaders = {}
    for group, loader in loader_list.items():
        if loader.initialized:
            loaders[group] = loader

    # Take into account user specified cpu list
    if cpu_list is None:
        cpu_list = loader_list["amr"].cpu_list if loader_list[
            "amr"].cpu_list is not None else range(1, data.meta["ncpu"] + 1)

    print("Processing {} files in {}".format(len(cpu_list), infile))

    # Allocate work arrays
    twotondim = 2**data.meta["ndim"]
    for loader in loaders.values():
        loader.allocate_buffers(ngridmax=data.meta["ngridmax"], twotondim=twotondim)

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
        for group, loader in loaders.items():
            fname = utils.generate_fname(nout, path, ftype=group, cpuid=cpu_num)
            with open(fname, mode='rb') as f:
                loader.bytes = f.read()

        # Read file headers
        for loader in loaders.values():
            loader.offsets.update(null_offsets)
            loader.read_header(data.meta)

        # Loop over levels
        for ilevel in range(data.meta["lmax"]):

            for loader in loaders.values():
                loader.read_level_header(ilevel, twotondim)

            # Loop over domains
            for domain in range(loaders["amr"].meta["nboundary"] + data.meta["ncpu"]):

                ncache = loaders["amr"].meta["ngridlevel"][domain, ilevel]

                for loader in loaders.values():
                    loader.read_domain_header()

                if ncache > 0:

                    if domain == cpu_num - 1:

                        for loader in loaders.values():
                            loader.read_cacheline_header(ncache, data.meta["ndim"])

                        for ind in range(twotondim):

                            # Read variables in cells
                            for loader in loaders.values():
                                loader.read_variables(ncache, ind, ilevel, cpu_num - 1,
                                                      data.meta)

                        # Apply selection criteria: select only leaf cells and
                        # add any criteria requested by the user via select.
                        conditions = {}
                        for loader in loaders.values():
                            conditions.update(loader.make_conditions(select, ncache))
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
                            for loader in loaders.values():
                                for item in loader.variables.values():
                                    item["pieces"][npieces] = item["buffer"][sel]

                        # Increment offsets with remainder of the file
                        for loader in loaders.values():
                            loader.read_footer(ncache, twotondim)

                    else:

                        for loader in loaders.values():
                            loader.step_over(ncache, twotondim, data.meta["ndim"])

    # Store the number of cells
    data.meta["ncells"] = ncells_tot

    # Merge all the data pieces into the Arrays
    for group in loaders.values():
        for key, item in group.variables.items():
            data[key] = np.concatenate(list(item["pieces"].values()))

    # If vector quantities are found, make them into vector Arrays
    utils.make_vector_arrays(data)

    # Create additional variables derived from the ones already loaded
    additional_variables(data)

    print("Total number of cells loaded: {}".format(ncells_tot))
    print("Memory used: {}".format(data.print_size()))

    return data
