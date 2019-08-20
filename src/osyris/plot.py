# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2019 Osyris contributors (https://github.com/nvaytet/osyris)
# @author Neil Vaytet

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections
from matplotlib.colors import LogNorm, Normalize
from . import config as conf


def perpendicular_vector(v):
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


def parse_arguments(args, args_osyris, args_plot):
    """
    Separate arguments to plotting functions between osyris specific arguments and
    Matplotlib arguments.
    """

    # Save a copy so we can exclude these parameters specific to osyris from the ones
    # to be sent to the matplotlib routines.
    args_osyris["cb_format"] = None
    ignore = set(args_osyris.keys())
    # Now we go through the arguments taken from the function call - scalar_args - and
    # add them to the osyris arguments.
    for key in args.keys():
        args_osyris[key] = args[key]
    # Define colorbar scaling
    cmap = args_osyris["cmap"]
    # norm = None
    # Default number of contours
    try:
        nc = args_osyris["nc"]
    except KeyError:
        nc = 21
    if (args_osyris["vmin"] is not None) and (args_osyris["vmin"] == args_osyris["vmax"]):
        args_osyris["vmin"] /= 1.1
        args_osyris["vmax"] *= 1.1
    norm = Normalize(vmin=args_osyris["vmin"], vmax=args_osyris["vmax"])
    # Default contour levels
    try:
        clevels = np.linspace(args_osyris["vmin"], args_osyris["vmax"], nc)
    except TypeError:
        clevels = [0.0, 1.0]
    # Perform log normalization if it is required
    try:
        if cmap.startswith("log") or cmap.endswith("log"):
            #cmap_save = cmap
            cmap = cmap.replace("log", "")
            chars = [",", ":", ";", " "]
            for ch in chars:
                cmap = cmap.replace(ch, "")
            if len(cmap) == 0:
                cmap = conf.default_values["colormap"]  # cmap_save
            #norm = LogNorm()
            norm = LogNorm(vmin=args_osyris["vmin"], vmax=args_osyris["vmax"])
            clevels = np.logspace(
                np.log10(args_osyris["vmin"]), np.log10(args_osyris["vmax"]), nc)
            args_osyris["cb_format"] = "%.1e"
    except AttributeError:
        pass

    # We now define default parameters for the matplotlib function
    keylist = args_plot.keys()
    if "levels" in keylist:
        args_plot["levels"] = clevels
    if "cmap" in keylist:
        args_plot["cmap"] = cmap
    if "norm" in keylist:
        args_plot["norm"] = norm
    # Then run through the vec_args, adding them to the plotting arguments, but
    # ignoring all osyris specific arguments
    keys = set(args.keys())
    for key in keys.difference(ignore):
        args_plot[key] = args[key]

    return


def get_slice_direction(holder, direction, dx=0, dy=0, origin=[0, 0, 0]):
    """
    Find direction vectors for slice
    """

    # Make it possible to call with only one size in the arguments
    if dy == 0.0:
        dy = dx

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
            #vel    = holder.get("velocity")[sphere]
            AngMom = np.sum(np.cross(pos, vel), axis=0)
            if view == "top":
                dir1 = AngMom
                # [1.0, 1.0, -1.0 * (dir1[0] + dir1[1]) / dir1[2]]
                dir2 = perpendicular_vector(dir1)
                dir3 = np.cross(dir1, dir2)
            elif view == "side":
                # Choose a vector perpendicular to the angular momentum vector
                dir3 = AngMom
                # [1.0, 1.0, -1.0 * (dir3[0] + dir3[1]) / dir3[2]]
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

    boxmin_x = np.nanmin(holder.get(dir_vecs[1][0]))
    boxmax_x = np.nanmax(holder.get(dir_vecs[1][0]))
    boxmin_y = np.nanmin(holder.get(dir_vecs[2][0]))
    boxmax_y = np.nanmax(holder.get(dir_vecs[2][0]))
    if dx+dy == 0.0:
        dx = boxmax_x - boxmin_x
        dy = boxmax_y - boxmin_y
    elif dx == 0.0:
        dx = dy

    for i in range(3):
        dir_vecs[i][1] /= np.linalg.norm(dir_vecs[i][1])

    box = [boxmin_x, boxmax_x, boxmin_y, boxmax_y]

    return dx, dy, box, dir_vecs, origin


def render_map(scalar=False, image=False, contour=False, scatter=False, vec=False, stream=False, outline=False, x=0, y=0,
               z_scal=0, z_imag=0, z_cont=0, z_outl=0, u_vect=0, v_vect=0, w_vect=0, u_strm=0, v_strm=0,
               w_strm=0, fname=None, axes=None, title=None, sinks=True, new_window=False,
               clear=True, block=False, xmin=0, xmax=0, ymin=0, ymax=0,
               resolution=128, scalar_args={}, image_args={}, contour_args={}, vec_args={},
               stream_args={}, scatter_args={}, outline_args={}, sink_args={}, dz=0, holder=None,
               thePlane=0, origin=[0, 0, 0], dir_vecs=[["z", [0, 0, 0]], ["x", [0, 0, 0]], ["y", [0, 0, 0]]], x_raw=None, y_raw=None, equal_axes=True,
               only_leafs=True):
    """
    Use matplotlib to plot histogram, slice or column density maps
    """

    if axes:
        theAxes = axes
    elif new_window:
        plt.figure()
        plt.subplot(111)
        theAxes = plt.gca()
    else:
        if clear:
            plt.clf()
        plt.subplot(111)
        theAxes = plt.gca()

    x += np.inner(origin, dir_vecs[1][1])
    y += np.inner(origin, dir_vecs[2][1])

    # Plot scalar field
    if scalar:

        # Round off AMR levels to integers
        if scalar.label == "level" or scalar.label == "cpu":
            z_scal = np.around(z_scal)
        # Parse scalar plot arguments
        scalar_args_osyris = {"vmin": np.nanmin(z_scal), "vmax": np.nanmax(
            z_scal), "cbar": True, "cbax": None, "cmap": conf.default_values["colormap"], "nc": 21}
        scalar_args_plot = {"levels": 1, "cmap": 1, "norm": 1}
        parse_arguments(scalar_args, scalar_args_osyris, scalar_args_plot)
        contf = theAxes.contourf(x, y, z_scal, **scalar_args_plot)
        if scalar_args_osyris["cbar"]:
            scb = plt.colorbar(
                contf, ax=theAxes, cax=scalar_args_osyris["cbax"], format=scalar_args_osyris["cb_format"])
            scb.ax.set_ylabel(
                scalar.label+(" ["+scalar.unit+"]" if len(scalar.unit) > 0 else ""))
            scb.ax.yaxis.set_label_coords(-1.1, 0.5)

    # Plot image
    if image:

        # Round off AMR levels to integers
        if image.label == "level" or image.label == "cpu":
            z_imag = np.around(z_imag)
        # Here we define a set of default parameters
        image_args_osyris = {"vmin": np.nanmin(z_imag), "vmax": np.nanmax(
            z_imag), "cbar": True, "cbax": None, "cmap": conf.default_values["colormap"], "nc": 21}
        # cmap and norm are just dummy arguments to tell the parsing function that they are required
        image_args_plot = {"cmap": 1, "norm": 1,
                           "interpolation": "none", "origin": "lower"}
        parse_arguments(image_args, image_args_osyris, image_args_plot)
        img = theAxes.imshow(
            z_imag, extent=[xmin, xmax, ymin, ymax], **image_args_plot)
        if image_args_osyris["cbar"]:
            icb = plt.colorbar(
                img, ax=theAxes, cax=image_args_osyris["cbax"], format=image_args_osyris["cb_format"])
            icb.ax.set_ylabel(
                image.label+(" ["+image.unit+"]" if len(image.unit) > 0 else ""))
            icb.ax.yaxis.set_label_coords(-1.1, 0.5)

    # Plot contours
    if contour:

        # Round off AMR levels to integers
        if contour.label == "level" or contour.label == "cpu":
            z_cont = np.around(z_cont)
        # Here we define a set of default parameters
        contour_args_osyris = {"vmin": np.nanmin(z_cont), "vmax": np.nanmax(z_cont), "cbar": False, "cbax": None,
                               "cmap": conf.default_values["colormap"], "nc": 21, "label": False, "fmt": "%1.3f"}
        # levels, cmap and norm are just dummy arguments to tell the parsing function that they are required
        contour_args_plot = {"levels": 1, "cmap": 1,
                             "norm": 1, "zorder": 10, "linestyles": "solid"}
        parse_arguments(contour_args, contour_args_osyris, contour_args_plot)
        cont = theAxes.contour(x, y, z_cont, **contour_args_plot)
        if contour_args_osyris["label"]:
            theAxes.clabel(cont, inline=1, fmt=contour_args_osyris["fmt"])
        if contour_args_osyris["cbar"]:
            ccb = plt.colorbar(
                cont, ax=theAxes, cax=contour_args_osyris["cbax"], format=contour_args_osyris["cb_format"])
            ccb.ax.set_ylabel(
                contour.label+(" ["+contour.unit+"]" if len(contour.unit) > 0 else ""))
            ccb.ax.yaxis.set_label_coords(-1.1, 0.5)

    # Plot scatter points
    if scatter:
        cube = np.where(np.logical_and(x_raw >= xmin, np.logical_and(x_raw <= xmax,
                                                                     np.logical_and(y_raw >= ymin, y_raw <= ymax))))
        try:
            vmin = np.nanmin(holder.get(
                scatter.name, only_leafs=only_leafs)[cube])
            vmax = np.nanmax(holder.get(
                scatter.name, only_leafs=only_leafs)[cube])
            scbar = True
            scmap = conf.default_values["colormap"]
        except AttributeError:
            vmin = vmax = None
            scbar = False
            scmap = None
        scatter_args_osyris = {"iskip": 1, "cmap": scmap, "vmin": vmin,
                               "vmax": vmax, "cbar": scbar, "cbax": None, "nc": 21}
        scatter_args_plot = {"cmap": 1, "marker": ".",
                             "c": "b", "edgecolor": "None", "s": 20, "norm": 1}
        parse_arguments(scatter_args, scatter_args_osyris, scatter_args_plot)
        # Check if a variable is given as a color
        try:
            scatter_args_plot["c"] = holder.get(scatter.name, only_leafs=only_leafs)[
                cube][::scatter_args_osyris["iskip"]]
        except AttributeError:
            pass
        scat = theAxes.scatter(x_raw[cube][::scatter_args_osyris["iskip"]],
                               y_raw[cube][::scatter_args_osyris["iskip"]], **scatter_args_plot)
        if scatter_args_osyris["cbar"]:
            rcb = plt.colorbar(
                scat, ax=theAxes, cax=scatter_args_osyris["cbax"], format=scatter_args_osyris["cb_format"])
            rcb.ax.set_ylabel(
                scatter.label+(" ["+scatter.unit+"]" if len(scatter.unit) > 0 else ""))
            rcb.ax.yaxis.set_label_coords(-1.1, 0.5)

    # Draw outline
    if outline:
        outline_args_plot = {"levels": [1.0], "colors": "grey"}
        for key in outline_args.keys():
            outline_args_plot[key] = outline_args[key]
        outl = theAxes.contour(x, y, z_outl, **outline_args_plot)

    # Plot vector field
    if vec:

        vec_args_osyris = {"vskip": int(0.047*resolution), "vscale": np.nanmax(w_vect), "vsize": 15.0, "vkey": True, "vkey_pos": [0.70, -0.08],
                           "cbar": False, "cbax": None, "vmin": np.nanmin(w_vect), "vmax": np.nanmax(w_vect), "nc": 21, "cmap": None, "colors": None,
                           "normalize_arrows": False}
        vec_args_plot = {"cmap": 1, "pivot": "mid", "color": "w", "norm": None}
        parse_arguments(vec_args, vec_args_osyris, vec_args_plot)
        vskip = vec_args_osyris["vskip"]
        if not "scale" in vec_args_plot.keys():
            vec_args_plot["scale"] = vec_args_osyris["vsize"] * \
                vec_args_osyris["vscale"]
        if vec_args_osyris["normalize_arrows"]:
            arrow_norm = np.sqrt(u_vect**2+v_vect**2)
            u_vect = u_vect/arrow_norm
            v_vect = v_vect/arrow_norm
        if vec_args_osyris["colors"]:
            vect = theAxes.quiver(x[::vskip], y[::vskip], u_vect[::vskip, ::vskip], v_vect[::vskip, ::vskip],
                                  w_vect[::vskip, ::vskip], **vec_args_plot)
            if vec_args_osyris["cbar"]:
                vcb = plt.colorbar(
                    vect, ax=theAxes, cax=vec_args_osyris["cbax"], orientation="horizontal", format=vec_args_osyris["cb_format"])
                vcb.ax.set_xlabel(vec_args_osyris["colors"].label+(
                    " ["+vec_args_osyris["colors"].unit+"]" if len(vec_args_osyris["colors"].unit) > 0 else ""))
        elif vec_args_plot["cmap"]:
            vect = theAxes.quiver(x[::vskip], y[::vskip], u_vect[::vskip, ::vskip], v_vect[::vskip, ::vskip],
                                  w_vect[::vskip, ::vskip], **vec_args_plot)
            if vec_args_osyris["cbar"]:
                vcb = plt.colorbar(
                    vect, ax=theAxes, cax=vec_args_osyris["cbax"], orientation="horizontal", format=vec_args_osyris["cb_format"])
                vcb.ax.set_xlabel(
                    vec.label+(" ["+vec.unit+"]" if len(vec.unit) > 0 else ""))
        else:
            vect = theAxes.quiver(x[::vskip], y[::vskip], u_vect[::vskip, ::vskip], v_vect[::vskip, ::vskip],
                                  **vec_args_plot)

        # Plot the scale of the vectors under the axes
        unit_u = vec.unit
        if vec_args_osyris["vkey"]:
            theAxes.quiverkey(vect, vec_args_osyris["vkey_pos"][0], vec_args_osyris["vkey_pos"][1],
                              vec_args_osyris["vscale"], "%.2f [%s]" % (vec_args_osyris["vscale"],
                                                                        unit_u), labelpos="E", labelcolor="k", coordinates="axes", color="k",
                              zorder=100)

    # Plot streamlines
    if stream:

        # Here we define a set of default parameters
        stream_args_osyris = {"cbar": False, "cbax": None, "sskip": 1, "vmin": np.nanmin(
            w_strm), "vmax": np.nanmax(w_strm), "nc": 21, "cmap": None}
        stream_args_plot = {"cmap": 1, "color": "w", "norm": None}
        parse_arguments(stream_args, stream_args_osyris, stream_args_plot)
        sskip = stream_args_osyris["sskip"]
        if stream_args_plot["cmap"]:
            stream_args_plot["color"] = w_strm[::sskip, ::sskip]
        strm = theAxes.streamplot(
            x[::sskip], y[::sskip], u_strm[::sskip, ::sskip], v_strm[::sskip, ::sskip], **stream_args_plot)
        if stream_args_osyris["cbar"]:
            scb = plt.colorbar(
                strm.lines, ax=theAxes, cax=stream_args_osyris["cbax"], orientation="horizontal", format=stream_args_osyris["cb_format"])
            scb.ax.set_xlabel(
                stream.label+(" ["+stream.unit+"]" if len(stream.unit) > 0 else ""))

    # Plot sink particles
    if holder.info["nsinks"] > 0 and sinks:
        dx = xmax-xmin
        if dz == 0:
            thickness = 0.05*dx
        else:
            thickness = 0.5*dz
        dist = (thePlane[0]*holder.sinks["x"]+thePlane[1]*holder.sinks["y"]+thePlane[2]*holder.sinks["z"]+thePlane[3]) \
            / np.sqrt(thePlane[0]**2 + thePlane[1]**2 + thePlane[2]**2)
        sinkcoords = np.transpose(
            [holder.sinks["x"], holder.sinks["y"], holder.sinks["z"]])
        sink_x = np.inner(sinkcoords, dir_vecs[1][1])
        sink_y = np.inner(sinkcoords, dir_vecs[2][1])
        subset = np.where(np.logical_and(dist <= thickness, np.logical_and(
            np.absolute(sink_x) <= 0.5*dx, np.absolute(sink_y) <= 0.5*dx)))
        srad = np.maximum(holder.sinks["radius"]
                          [subset], np.full(len(subset), dx*0.01))
        xy = np.array([sink_x[subset], sink_y[subset]]).T
        patches = [plt.Circle(cent, size) for cent, size in zip(xy, srad)]

        if "colors" in sink_args.keys():
            try:
                # This is a custom array of numbers
                sink_colors = np.array(sink_args["colors"])[subset] + 0.0
            except (TypeError, IndexError):
                # Go and find values in sinks dict
                sink_colors = holder.sinks[sink_args["colors"]][subset]
            sk_vmin = np.nanmin(sink_colors)
            sk_vmax = np.nanmax(sink_colors)
            sk_cmap = conf.default_values["colormap"]
            sk_cbar = True
        else:
            sk_vmin = 0
            sk_vmax = 1
            sk_cmap = None
            sk_cbar = False
        sink_args_osyris = {"cbar": sk_cbar, "cbax": None, "vmin": sk_vmin,
                            "vmax": sk_vmax, "cmap": sk_cmap, "colors": None}
        sink_args_plot = {"facecolors": "w", "edgecolors": "k",
                          "linewidths": 2, "alpha": 0.7, "cmap": 1, "norm": 1}
        parse_arguments(sink_args, sink_args_osyris, sink_args_plot)
        coll = matplotlib.collections.PatchCollection(
            patches, **sink_args_plot)
        if "colors" in sink_args.keys():
            coll.set_array(np.array(sink_colors))
            coll.set_clim([sink_args_osyris["vmin"], sink_args_osyris["vmax"]])
        theAxes.add_collection(coll)
        if sink_args_osyris["cbar"]:
            skcb = plt.colorbar(coll, ax=theAxes, cax=sink_args_osyris["cbax"])
            skcb.ax.set_ylabel("Sink "+sink_args["colors"])
            skcb.ax.yaxis.set_label_coords(-1.1, 0.5)

    try:
        title += ""
        theAxes.set_title(title)
    except TypeError:
        theAxes.set_title("Time = %.3f %s" % (
            holder.info["time"]/conf.constants[conf.default_values["time_unit"]], conf.default_values["time_unit"]))

    if axes:
        pass
    elif new_window:
        theAxes.set_xlim([xmin, xmax])
        theAxes.set_ylim([ymin, ymax])
    else:
        if clear:
            theAxes.set_xlim([xmin, xmax])
            theAxes.set_ylim([ymin, ymax])
        else:
            theAxes.set_xlim([min(theAxes.get_xlim()[0], xmin),
                              max(theAxes.get_xlim()[1], xmax)])
            theAxes.set_ylim([min(theAxes.get_ylim()[0], ymin),
                              max(theAxes.get_ylim()[1], ymax)])

    # Define axes labels
    xlab = getattr(holder, dir_vecs[1][0]).label
    if len(getattr(holder, dir_vecs[1][0]).unit) > 0:
        xlab += " ["+getattr(holder, dir_vecs[1][0]).unit+"]"
    ylab = getattr(holder, dir_vecs[2][0]).label
    if len(getattr(holder, dir_vecs[2][0]).unit) > 0:
        ylab += " ["+getattr(holder, dir_vecs[2][0]).unit+"]"
    theAxes.set_xlabel(xlab)
    theAxes.set_ylabel(ylab)

    if equal_axes:
        theAxes.set_aspect("equal")

    if fname:
        plt.savefig(fname, bbox_inches="tight")
    elif axes:
        pass
    else:
        plt.show(block=block)

    return
