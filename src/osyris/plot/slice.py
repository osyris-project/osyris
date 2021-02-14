# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2019 Osyris contributors (https://github.com/nvaytet/osyris)
# @author Neil Vaytet

import numpy as np
from ..core.tools import perpendicular_vector


def get_slice_direction(direction=None,
                        parent=None,
                        dx=None,
                        origin=[0, 0, 0]):
    """
    Find direction vectors for slice.

    Direction can be:


    The origin can be either a vector of 3 numbers (xyz), or it can be "sink17"
    for sink particles.
    """

    # # Transform origin to coordinates if sink is requested
    # try:
    #     if origin.startswith("sink"):
    #         isink = np.where(holder.sinks["id"] == int(origin.split(":")[1]))[0][0]
    #         origin = [holder.sinks["x"][isink], holder.sinks["y"]
    #                   [isink], holder.sinks["z"][isink]]
    # except AttributeError:
    #     pass

    dir_list = {"x": [1, 0, 0], "y": [0, 1, 0], "z": [0, 0, 1]}
    dir_type = len(np.shape(direction))

    if direction in ["auto", "top", "side"]:
        sphere_rad = 0.5 * dx
        xyz = parent["xyz"] - origin
        # x_loc = parent["x"] - origin[0]
        # y_loc = parent["y"] - origin[1]
        # z_loc = parent["z"] - origin[2]
        # r_loc = np.linalg.norm([x_loc, y_loc, z_loc], axis=0)
        # print("r_loc", r_loc)
        # Compute angular momentum vector
        sphere = np.where(xyz.norm < sphere_rad.magnitude)
        # sphere = np.where(xyz.values < sphere_rad)
        # print(x_loc[sphere])
        # print(sphere)
        # print(sphere.shape)
        # print("===")
        # print((x_loc[sphere], y_loc[sphere],
        #                  z_loc[sphere])*parent["mass"][sphere])
        # pos = np.array([x_loc[sphere], y_loc[sphere],
        #                  z_loc[sphere]]).T * parent["mass"].values[sphere]
        pos = xyz * parent["mass"]  #.values[sphere]
        print(pos.shape)
        vel = parent["velocity"]  # [sphere]
        print(vel.shape)

        # vel = np.vstack((parent["velocity_x"][sphere],
        #                  parent["velocity_y"][sphere], holder.get("velocity_z")[sphere])).T
        AngMom = np.sum(np.cross(pos.array[sphere], vel.array[sphere]), axis=0)
        if direction == "side":
            # Choose a vector perpendicular to the angular momentum vector
            dir3 = AngMom
            dir1 = perpendicular_vector(dir3)
            dir2 = np.cross(dir1, dir3)
        else:
            dir1 = AngMom
            dir2 = perpendicular_vector(dir1)
            dir3 = np.cross(dir1, dir2)
        # elif view ==
        # else:
        #     raise ValueError("Unknown view direction.")
        norm1 = np.linalg.norm(dir1)
        print("Normal slice vector: [%.5e,%.5e,%.5e]" %
              (dir1[0] / norm1, dir1[1] / norm1, dir1[2] / norm1))
        # dir_vecs = [["z", dir1], ["x", dir2], ["y", dir3]]
        dir_vecs = np.array([dir1, dir2, dir3])

    elif isinstance(direction, str):
        if len(direction) == 3:  # This is the case where direction = "xyz"
            # dir_vecs = [[direction[0], dir_list[direction[0]]],
            #             [direction[1], dir_list[direction[1]]],
            #             [direction[2], dir_list[direction[2]]]]
            dir_vecs = np.array([
                dir_list[direction[0]], dir_list[direction[1]],
                dir_list[direction[2]]
            ])
        elif direction == "x":
            # dir_vecs = [["x", dir_list["x"]], [
            #     "y", dir_list["y"]], ["z", dir_list["z"]]]
            dir_vecs = np.array([dir_list["x"], dir_list["y"], dir_list["z"]])
        elif direction == "y":
            dir_vecs = np.array([dir_list["y"], dir_list["z"], dir_list["x"]])
            # dir_vecs = [["y", dir_list["y"]], [
            #     "z", dir_list["z"]], ["x", dir_list["x"]]]
        elif direction == "z":
            dir_vecs = np.array([dir_list["z"], dir_list["x"], dir_list["y"]])
            # dir_vecs = [["z", dir_list["z"]], [
            #     "x", dir_list["x"]], ["y", dir_list["y"]]]
    # This is the case where direction = [1,1,2] (i.e. is a vector with 3 numbers)
    elif len(direction) == 3:
        dir1 = direction
        dir2 = perpendicular_vector(dir1)
        dir3 = np.cross(dir1, dir2).tolist()
        dir_vecs = np.array([dir1, dir2, dir3])
        # dir_vecs = [["z", dir1], ["x", dir2], ["y", dir3]]
    # This is the case where two vectors are specified: direction = [[1,0,1],[0,1,0]]
    elif len(direction) == 2:
        dir_vecs = np.array(
            [direction[0], direction[1],
             np.cross(direction[0], direction[1])])
        # dir_vecs = [["z", direction[0]],
        #             ["x", direction[1]],
        #             ["y", np.cross(direction[0], direction[1]).tolist()]]
    else:
        print("Bad direction for slice: ", direction)
        return

    dir_vecs = dir_vecs / np.linalg.norm(dir_vecs, axis=1).reshape(3, 1)
    return dir_vecs, origin
