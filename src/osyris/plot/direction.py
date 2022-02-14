# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np
from ..core import Array


def _perpendicular_vector(v):
    """
    Compute a vector perpendicular to the input vector
    """

    # x = y = z = 0 is not an acceptable solution
    if v[0] == v[1] == v[2] == 0:
        raise ValueError("zero-vector")

    if v[2] == 0:
        return [-v[1], v[0], 0]
    else:
        return [1.0, 1.0, -1.0 * (v[0] + v[1]) / v[2]]


def get_direction(direction=None, dataset=None, dx=None, dy=None, origin=None):
    """
    Find direction vectors for slice.

    Direction can be:


    The origin can be either a vector of 3 numbers (xyz), or it can be "sink17"
    for sink particles.
    """

    # # Transform origin to coordinates if sink is requested
    # try:
    #     if origin.startswith("sink"):
    #         isink = np.where(holder.sinks["id"] ==
    #          int(origin.split(":")[1]))[0][0]
    #         origin = [holder.sinks["x"][isink], holder.sinks["y"]
    #                   [isink], holder.sinks["z"][isink]]
    # except AttributeError:
    #     pass

    ndim = dataset.meta["ndim"]

    dir_list = {"x": [1, 0, 0], "y": [0, 1, 0], "z": [0, 0, 1]}
    lab_list = {"x": "pos_x", "y": "pos_y", "z": "pos_z"}
    dir_labs = {"x": "pos_u", "y": "pos_v", "z": "pos_w"}

    if ndim < 3:
        dir_vecs = np.array([[0., 0.], [1., 0.], [0., 1.]])
        dir_labs = {"x": "pos_x", "y": "pos_y", "z": "pos_z"}

    elif direction in ["top", "side"]:
        # sphere_rad = 0.0 * dataset["amr"]["position"].unit.units
        if dx is None:
            pos = dataset["amr"]["position"].values
            sphere_rad = (0.5 * (pos[:, 0].max() - pos[:, 0].min() + pos[:, 1].max() -
                                 pos[:, 1].min() + pos[:, 2].max() - pos[:, 2].min()) /
                          3.) * dataset["amr"]["position"].unit.units
        else:
            sphere_rad = 0.25 * (dx + dy)

        #    sphere_rad
        # or dy is None:
        #    raise RuntimeError("When using automatic slice orientation, "
        #                       "dx cannot be None.")
        # sphere_rad = 0.25 * (dx + dy)
        xyz = dataset["amr"]["position"] - origin
        # if origin is not None:
        #     xyz -= origin
        # Compute angular momentum vector
        # sphere = np.where(xyz.norm < sphere_rad.magnitude)
        sphere = xyz.norm < sphere_rad.magnitude
        pos = xyz * dataset["hydro"]["mass"]
        vel = dataset["hydro"]["velocity"]

        AngMom = np.sum(np.cross(pos.array[sphere], vel.array[sphere]), axis=0)
        if direction == "side":
            # Choose a vector perpendicular to the angular momentum vector
            dir3 = AngMom
            dir1 = _perpendicular_vector(dir3)
            dir2 = np.cross(dir1, dir3)
        else:
            dir1 = AngMom
            dir2 = _perpendicular_vector(dir1)
            dir3 = np.cross(dir1, dir2)
        norm1 = np.linalg.norm(dir1)
        print("Normal slice vector: [%.5e,%.5e,%.5e]" %
              (dir1[0] / norm1, dir1[1] / norm1, dir1[2] / norm1))
        dir_vecs = np.array([dir1, dir2, dir3])

    elif isinstance(direction, str):
        if len(direction) == 3:  # This is the case where direction = "xyz"
            dir_vecs = np.array([
                dir_list[direction[0]], dir_list[direction[1]], dir_list[direction[2]]
            ])
            dir_labs = {
                "x": f"pos_{direction[1]}",
                "y": f"pos_{direction[2]}",
                "z": f"pos_{direction[0]}"
            }
        elif direction == "x":
            dir_vecs = np.array([dir_list["x"], dir_list["y"], dir_list["z"]])
            dir_labs = {"x": "pos_y", "y": "pos_z", "z": "pos_x"}
            # dir_labs = [dir_labs["x"], dir_labs["y"], dir_labs["z"]]
        elif direction == "y":
            dir_vecs = np.array([dir_list["y"], dir_list["z"], dir_list["x"]])
            # dir_labs = [dir_labs["y"], dir_labs["z"], dir_labs["x"]]
            dir_labs = {"x": "pos_z", "y": "pos_x", "z": "pos_y"}
        elif direction == "z":
            dir_vecs = np.array([dir_list["z"], dir_list["x"], dir_list["y"]])
            dir_labs = {"x": "pos_x", "y": "pos_y", "z": "pos_z"}
            # dir_labs = [dir_labs["z"], dir_labs["x"], dir_labs["y"]]
    # This is the case where direction = [1,1,2]
    # (i.e. is a vector with 3 numbers)
    elif len(direction) == 3:
        dir1 = direction
        dir2 = _perpendicular_vector(dir1)
        dir3 = np.cross(dir1, dir2).tolist()
        dir_vecs = np.array([dir1, dir2, dir3])
    # This is the case where two vectors are specified:
    # direction = [[1,0,1],[0,1,0]]
    elif len(direction) == 2:
        dir_vecs = np.array(
            [direction[0], direction[1],
             np.cross(direction[0], direction[1])])
    else:
        raise ValueError(f"Bad direction for slice: {direction}.")

    # Avoid division by zero in norm
    norm = np.linalg.norm(dir_vecs, axis=1).reshape(3, 1)
    norm[norm == 0.] = 1.0
    dir_vecs = dir_vecs / norm

    # if origin is None:
    #     origin = Array(values=np.zeros([1, ndim]), unit=dataset["amr"]["position"].unit)

    return dir_vecs, dir_labs
