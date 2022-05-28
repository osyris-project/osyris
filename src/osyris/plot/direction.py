# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)

from ..core import Vector

import numpy as np


def _perpendicular_vector(v):
    """
    Compute a vector perpendicular to the input vector
    """

    if v.z.values == 0:
        return Vector(-v.y.values, v.x.values, 0, unit=v.unit)
    else:
        return Vector(1.0, 1.0, (-1.0 * (v.x + v.y) / v.z).values, unit=v.unit)


def get_direction(direction=None, dataset=None, dx=None, dy=None, origin=None):
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
        "z": Vector(0, 0, 1, name="z")
    }
    dir_names = ("normal", "pos_u", "pos_v")

    if dataset.meta["ndim"] < 3:
        return {
            dir_names[0]: Vector(0, 0, name="z"),
            dir_names[1]: Vector(1, 0, name="x"),
            dir_names[2]: Vector(0, 1, name="y")
        }

    if isinstance(direction, str):
        direction = direction.lower()
        if direction in ["top", "side"]:
            if dx is None:
                pos = dataset["amr"]["position"]
                sphere_rad = (0.5 * (pos.x.max() - pos.x.min() + pos.y.max() -
                                     pos.y.min() + pos.z.max() - pos.z.min()) /
                              3.) * dataset["amr"]["position"].unit
            else:
                sphere_rad = 0.25 * (dx + dy)

            xyz = dataset["amr"]["position"]
            if origin is not None:
                xyz = xyz - origin
            # Compute angular momentum vector
            sphere = (xyz.norm < sphere_rad).values
            pos = xyz * dataset["hydro"]["mass"]
            vel = dataset["hydro"]["velocity"]

            ang_mom = np.sum(pos[sphere].cross(vel[sphere]))
            if direction == "side":
                # Choose a vector perpendicular to the angular momentum vector
                dir3 = ang_mom
                dir1 = _perpendicular_vector(dir3)
                dir2 = dir3.cross(dir1)
            else:
                dir1 = ang_mom
                dir2 = _perpendicular_vector(dir1)
                dir3 = dir1.cross(dir2)
            dir_vecs = {}
            print("Basis vectors:")
            for key, vec in zip(dir_names, (dir1, dir2, dir3)):
                vec.name = key
                dir_vecs[key] = vec
                print(vec)

        if set(direction) == set("xyz"):  # case where direction = "xyz", "zyx" etc.
            return {
                "normal": dir_list[direction[0]],
                "pos_u": dir_list[direction[1]],
                "pos_v": dir_list[direction[2]]
            }
        elif direction == "x":
            return {
                "normal": dir_list["x"],
                "pos_u": dir_list["y"],
                "pos_v": dir_list["z"]
            }
        elif direction == "y":
            return {
                "normal": dir_list["y"],
                "pos_u": dir_list["z"],
                "pos_v": dir_list["x"]
            }
        elif direction == "z":
            return {
                "normal": dir_list["z"],
                "pos_u": dir_list["x"],
                "pos_v": dir_list["y"]
            }
    elif isinstance(direction, Vector):
        dir_vecs = {"normal": direction}
        u = _perpendicular_vector(dir_vecs["normal"])
        u.name = "pos_u"
        dir_vecs["pos_u"] = u
        v = direction.cross(dir_vecs["pos_u"])
        v.name = "pos_v"
        dir_vecs["pos_v"] = v
    else:
        raise ValueError(f"Bad direction for slice: {direction}.")

    for key, vec in dir_vecs.items():
        v = vec / vec.norm
        v.name = vec.name
        dir_vecs[key] = v

    return dir_vecs
