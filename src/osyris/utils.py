import numpy as np

def get_slice_direction(direction, origin=[0, 0, 0]):
    """
    Find direction vectors for slice
    """

    # # Make it possible to call with only one size in the arguments
    # if dy == 0.0:
    #     dy = dx

    # Transform origin to coordinates if sink is requested
    try:
        if origin.startswith("sink"):
            isink = np.where(holder.sinks["id"] == int(origin.split(":")[1]))[0][0]
            origin = [holder.sinks["x"][isink], holder.sinks["y"]
                      [isink], holder.sinks["z"][isink]]
    except AttributeError:
        pass

    dir_list = {"x": [1, 0, 0], "y": [0, 1, 0], "z": [0, 0, 1]}
    dir_type = len(np.shape(direction))

    if dir_type == 0:  # This is the case where direction contains just one character "x", "y" or "z"
        # This is the case where direction = "auto"
        if direction.startswith("auto"):
            params = direction.split(":")
            if len(params) == 1:
                view = "top"
            else:
                view = params[1]
            if len(params) < 3:
                sphere_rad = 0.5 * \
                    ((np.nanmax(holder.get("x"))-np.nanmin(holder.get("x")))
                     if dx == 0.0 else dx)
            else:
                sphere_rad = float(params[2])
            x_loc = holder.get("x") - origin[0]
            y_loc = holder.get("y") - origin[1]
            z_loc = holder.get("z") - origin[2]
            r_loc = np.linalg.norm([x_loc, y_loc, z_loc], axis=0)
            # Compute angular momentum vector
            sphere = np.where(r_loc < sphere_rad)
            pos = np.vstack((x_loc[sphere], y_loc[sphere],
                             z_loc[sphere])*holder.get("mass")[sphere]).T
            vel = np.vstack((holder.get("velocity_x")[sphere], holder.get(
                "velocity_y")[sphere], holder.get("velocity_z")[sphere])).T
            AngMom = np.sum(np.cross(pos, vel), axis=0)
            if view == "top":
                dir1 = AngMom
                dir2 = perpendicular_vector(dir1)
                dir3 = np.cross(dir1, dir2)
            elif view == "side":
                # Choose a vector perpendicular to the angular momentum vector
                dir3 = AngMom
                dir1 = perpendicular_vector(dir3)
                dir2 = np.cross(dir1, dir3)
            else:
                print("Unknown view direction")
                return
            norm1 = np.linalg.norm(dir1)
            print("Normal slice vector: [%.5e,%.5e,%.5e]" % (
                dir1[0]/norm1, dir1[1]/norm1, dir1[2]/norm1))
            dir_vecs = [["z", dir1], ["x", dir2], ["y", dir3]]
        elif len(direction) == 3:  # This is the case where direction = "xyz"
            dir_vecs = [[direction[0], dir_list[direction[0]]],
                        [direction[1], dir_list[direction[1]]],
                        [direction[2], dir_list[direction[2]]]]
        elif direction == "x":
            dir_vecs = [["x", dir_list["x"]], [
                "y", dir_list["y"]], ["z", dir_list["z"]]]
        elif direction == "y":
            dir_vecs = [["y", dir_list["y"]], [
                "z", dir_list["z"]], ["x", dir_list["x"]]]
        elif direction == "z":
            dir_vecs = [["z", dir_list["z"]], [
                "x", dir_list["x"]], ["y", dir_list["y"]]]
    # This is the case where direction = [1,1,2] (i.e. is a vector with 3 numbers)
    elif dir_type == 1:
        dir1 = direction
        dir2 = perpendicular_vector(dir1)
        dir3 = np.cross(dir1, dir2).tolist()
        dir_vecs = [["z", dir1], ["x", dir2], ["y", dir3]]
    # This is the case where two vectors are specified: direction = [[1,0,1],[0,1,0]]
    elif dir_type == 2:
        dir_vecs = [["z", direction[0]],
                    ["x", direction[1]],
                    ["y", np.cross(direction[0], direction[1]).tolist()]]
    else:
        print("Bad direction for slice: ", direction)
        return

    # boxmin_x = np.nanmin(holder.get(dir_vecs[1][0]))
    # boxmax_x = np.nanmax(holder.get(dir_vecs[1][0]))
    # boxmin_y = np.nanmin(holder.get(dir_vecs[2][0]))
    # boxmax_y = np.nanmax(holder.get(dir_vecs[2][0]))
    # if dx+dy == 0.0:
    #     dx = boxmax_x - boxmin_x
    #     dy = boxmax_y - boxmin_y
    # elif dx == 0.0:
    #     dx = dy

    for i in range(3):
        dir_vecs[i][1] /= np.linalg.norm(dir_vecs[i][1])

    # box = [boxmin_x, boxmax_x, boxmin_y, boxmax_y]

    # return dx, dy, box, dir_vecs, origin
    return dir_vecs, origin


def create_vector_containers(df):
    """
    Create dummy variables containing the components of the vectors
    """
    if df.attrs["ndim"] > 1:
        for key in df:
            if key.endswith("_x"):
                rawkey = key[:-2]
                ok = rawkey+"_y" in df
                # try:
                #     k1 = len(self.get(rawkey+"_y"))
                # except AttributeError:
                #     ok = False
                if df.attrs["ndim"] > 2:
                    ok = ok and rawkey+"_z" in df
                    # try:
                    #     k2 = len(self.get(rawkey+"_z"))
                    # except AttributeError:
                    #     ok = False

                if ok:
                    # self.vector_field(name=rawkey,label=rawkey)
                    df[rawkey] = np.linalg.norm(
                        np.array([df[rawkey+'_x'].values,
                                  df[rawkey+'_y'].values,
                                  df[rawkey+'_z'].values]).T, axis=1)
                    df[rawkey].attrs["vector"] = {'x': df[rawkey+'_x'],
                                                  'y': df[rawkey+'_y'],
                                                  'z': df[rawkey+'_z']}

    return
