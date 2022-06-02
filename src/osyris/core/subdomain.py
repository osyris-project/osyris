# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)

import numpy as np
from .dataset import Dataset
from .datagroup import Datagroup
from .array import Array
from .vector import Vector
from .. import spatial
from .. import units


class Subdomain(Dataset):
    def __init__(self,
                 dataset,
                 selection,
                 dr=None,
                 dx=None,
                 dy=None,
                 dz=None,
                 origin=None,
                 basis=None,
                 dr_L=None):

        self._parent = dataset  # store reference to parent dataset
        self.meta = dataset.meta
        self.units = dataset.units

        self._sink_particles = False
        self._particles = False

        if origin is None:
            self.origin = Vector(x=0, y=0, z=0, unit=dataset.units['length'].units)
        else:
            self.origin = origin

        # find subdomain and extract data
        valid_comp_amr, valid_comp_sinks, valid_comp_particles = self._determine_subdomain(
            selection, dr, dx, dy, dz)
        self._extract_data(valid_comp_amr, valid_comp_sinks, valid_comp_particles)

        # translate positions in subdomain with new origin
        spatial.change_origin(self, self.origin)

        if basis is None:
            self.basis = [0, 0, 1]  # assume it's the ramses grid
        elif isinstance(basis, str):
            if dr_L is None:
                raise ValueError(
                    "Please provide the radius size with which to compute angular momentum (dr_L)"
                )
            ang_mom = spatial.get_ang_mom(self, dr_L)
            if basis.lower() == "top":
                basis = ang_mom / ang_mom.norm
            elif basis.lower() == "side":
                perp_v = Vector(1.0,
                                1.0,
                                (-1.0 * (ang_mom.x + ang_mom.y) / ang_mom.z).values,
                                unit=ang_mom.unit)
                basis = perp_v / perp_v.norm
            spatial.change_basis(self, basis)
        else:
            spatial.change_basis(self, basis)

    def _determine_subdomain(self, selection, dr, dx, dy, dz):
        # find indices of data within subdomain for amr dependent groups, sinks and particles
        centered_pos = self._parent["amr"]["position"] - self.origin
        if "sink" in self._parent:
            if len(self._parent["sink"]) > 0:
                self._sink_particles = True
                centered_sinks_pos = self._parent["sink"]["position"] - self.origin
        if "part" in self._parent:
            if len(self._parent["part"]) > 0:
                self._particles = True
                centered_particles_pos = self._parent["part"]["position"] - self.origin
        if isinstance(selection, str):
            if selection.lower() == "spherical":
                if dr is None:
                    raise ValueError("Please specify a valid selection radius")
                valid_comp_amr = (centered_pos.norm <= dr).values
                valid_comp_sinks = (centered_sinks_pos.norm <= dr).values
                valid_comp_particles = (centered_particles_pos.norm <= dr).values
            elif selection.lower() == "cubic":
                if sum(v is not None for v in [dx, dy, dz]) == 3:
                    # find amr indices
                    valid_comp_x_amr = (centered_pos.x >= -dx * .5) & (centered_pos.x <=
                                                                       dx * .5)
                    valid_comp_y_amr = (centered_pos.y >= -dy * .5) & (centered_pos.y <=
                                                                       dy * .5)
                    valid_comp_z_amr = (centered_pos.z >= -dz * .5) & (centered_pos.z <=
                                                                       dz * .5)
                    valid_comp_amr = (valid_comp_x_amr & valid_comp_y_amr
                                      & valid_comp_z_amr).values
                    # find sink & particles indices
                    if self._sink_particles:
                        valid_comp_x_sinks = (centered_sinks_pos.x >= -dx * .5) & (
                            centered_sinks_pos.x <= dx * .5)
                        valid_comp_y_sinks = (centered_sinks_pos.y >= -dy * .5) & (
                            centered_sinks_pos.y <= dy * .5)
                        valid_comp_z_sinks = (centered_sinks_pos.z >= -dz * .5) & (
                            centered_sinks_pos.z <= dz * .5)
                        valid_comp_sinks = (valid_comp_x_sinks & valid_comp_y_sinks
                                            & valid_comp_z_sinks).values
                        if not np.any(valid_comp_sinks):
                            self._sink_particles = False
                    if self._particles:
                        valid_comp_x_particles = (centered_particles_pos.x >=
                                                  -dx * .5) & (centered_particles_pos.x
                                                               <= dx * .5)
                        valid_comp_y_particles = (centered_particles_pos.y >=
                                                  -dy * .5) & (centered_particles_pos.y
                                                               <= dy * .5)
                        valid_comp_z_particles = (centered_particles_pos.z >=
                                                  -dz * .5) & (centered_particles_pos.z
                                                               <= dz * .5)
                        valid_comp_particles = (valid_comp_x_particles
                                                & valid_comp_y_particles
                                                & valid_comp_z_particles).values
                        if not np.any(valid_comp_particles):
                            self._particles = False
                else:
                    raise ValueError("Please specify a valid dx, dy and dz")
            else:
                raise ValueError(
                    "Unrecognized string '{}', valid criterions are 'spherical' and 'cubic'"
                    .format(selection.lower()))
        elif isinstance(selection, np.array):
            raise NotImplementedError(
                "Numpy array indexing for subdomain extraction coming soon...")
        else:
            raise ValueError(
                "Please specify a valid selection criterion (either 'spherical' or 'cubic')"
            )

        if not np.any(valid_comp_amr):
            raise ValueError("Empty domain")

        return valid_comp_amr, valid_comp_sinks, valid_comp_particles

    def _extract_data(self, valid_comp_amr, valid_comp_sinks, valid_comp_particles):
        # proceed to data extraction from parent dataset to subdomain
        self.groups = {}
        # extract amr dependent groups
        for group in ["hydro", "amr", "grav"]:
            self.groups[group] = Datagroup(parent=self, name=group)
            for element in self._parent[group]:
                self.groups[group][element] = self._parent[group][element][
                    valid_comp_amr]
        # extract sinks & particles
        if self._sink_particles:
            self.groups["sink"] = Datagroup(parent=self, name="sink")
            for element in self._parent["sink"]:
                self.groups["sink"][element] = self._parent["sink"][element][
                    valid_comp_sinks]
        if self._particles:
            self.groups["part"] = Datagroup(parent=self, name="part")
            for element in self._parent["part"]:
                self.groups["part"][element] = self._parent["part"][element][
                    valid_comp_particles]
