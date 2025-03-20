# SPDX-License-Identifier: BSD-3-Clause

import numpy as np

from ..core.vector import Vector, VectorBasis


def _basis_with_names(basis):
    basis.n.name = "normal"
    basis.u.name = "pos_u"
    basis.v.name = "pos_v"
    return basis


def get_direction(direction, data=None, dx=None, dy=None, origin=None):
    """
    Find direction vectors for slice.

    Direction can be:

    - a single letter ``"x"``, ``"y"``, ``"z"`` representing the normal to the plane
    - three letters, e.g. ``"zyx"``, where the first letter is the normal to the plane,
      and the following two letters are the two other perpendicual vectors
    - a list of 3 numbers, representing the 3 components of the normal to the plane
    - ``"top"`` for a 'top' view according to local angular momentum
    - ``"side"`` for a 'side' view according to local angular momentum
    - a list of 2 lists of 3 numbers, representing the normal to the plane and the first
      of the two other perpendicular vectors (the third one will be computed from the
      other two)

    The origin is a vector of 3 numbers (xyz).
    """
    dir_list = {
        "x": Vector(1, 0, 0, name="x"),
        "y": Vector(0, 1, 0, name="y"),
        "z": Vector(0, 0, 1, name="z"),
    }

    if isinstance(direction, str):
        direction = direction.lower()
        if direction in ["top", "side"]:
            pos = data["position"]
            if dx is None:
                sphere_rad = (
                    0.5
                    * (
                        pos.x.max()
                        - pos.x.min()
                        + pos.y.max()
                        - pos.y.min()
                        + pos.z.max()
                        - pos.z.min()
                    )
                    / 3.0
                )
            else:
                sphere_rad = 0.25 * (dx + dy)

            if origin is not None:
                pos = pos - origin

            # Compute angular momentum vector
            sphere = (pos.norm < sphere_rad).values
            weighted_pos = pos[sphere] * data["mass"][sphere]
            vel = data["velocity"][sphere]
            ang_mom = np.sum(weighted_pos.cross(vel))
            ang_mom.unit = ""

            basis = VectorBasis(n=ang_mom)
            if direction == "side":
                basis = basis.roll()
            print(basis)
            return _basis_with_names(basis)

        if set(direction) == set("xyz"):  # case where direction = "xyz", "zyx" etc.
            return VectorBasis(
                n=dir_list[direction[0]],
                u=dir_list[direction[1]],
                v=dir_list[direction[2]],
            )
        elif direction == "x":
            return VectorBasis(n=dir_list["x"], u=dir_list["y"], v=dir_list["z"])
        elif direction == "y":
            return VectorBasis(n=dir_list["y"], u=dir_list["z"], v=dir_list["x"])
        elif direction == "z":
            return VectorBasis(n=dir_list["z"], u=dir_list["x"], v=dir_list["y"])
    elif isinstance(direction, Vector):
        return _basis_with_names(VectorBasis(n=direction))
    elif isinstance(direction, VectorBasis):
        return _basis_with_names(
            VectorBasis(n=direction.n, u=direction.u, v=direction.v)
        )
    else:
        raise ValueError(f"Bad direction for slice: {direction}.")
