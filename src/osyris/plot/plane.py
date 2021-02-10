# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2019 Osyris contributors (https://github.com/nvaytet/osyris)
# @author Neil Vaytet

import numpy as np
from .slice import get_slice_direction
from .render import render
from .parser import parse_layer
from ..core import Plot, Array



def plane(*layers,
    # scalar=False, image=False, contour=False, vec=False, stream=False, axes=None,
               direction="z", dx=0.0, dy=0.0, fname=None, title=None, sinks=True, copy=False,
               mode=None, norm=None, vmin=None, vmax=None,operation="mean",
               origin=[0, 0, 0],
               resolution=256,
               # new_window=False, update=None,
               # clear=True, plot=True, block=False, interpolation="linear", sink_args={},
               # scalar_args={}, image_args={}, contour_args={}, vec_args={}, stream_args={},
               # outline=False, outline_args={}, lmax=0, slice_direction=None):
               ax=None,
               **kwargs):
    """
    Plot a 2D slice through the data cube. The arguments are:
    - scalar     : the scalar field to be plotted, e.g. mydata.density
    - image      : the scalar field to be plotted with an image
    - contour    : the scalar field to be plotted with contours
    - vec        : the vector field to be overplotted, e.g. mydata.velocity
    - stream     : the field for streamlines to be overplotted, e.g. mydata.B
    - direction  : the direction normal to the plane of the slice
    - dx         : the x extent of the slice, in units of scale (see data loader)
    - dy         : the y extent of the slice, in units of scale. If not specified, dy = dx
    - dz         : the thickness of the slice
    - axes       : if specified, the data is plotted on the specified axes (see demo).
    - resolution : number of pixels in the slice.
    - fname      : if specified, the figure is saved to file.
    """
    # print(layers)

    if isinstance(layers, Array):
        layers = [layers]

    to_process = []
    to_render = []
    for layer in layers:
        data, settings, params = parse_layer(layer, mode=mode, norm=norm,
            vmin=vmin, vmax=vmax, operation=operation, **kwargs)
        # print(data, settings, params)
        # if settings["mode"] in ["vec", "stream"]:
        #     to_process[data.name] = data._array
        # else:



# if holder.info["ndim"] < 3:
    #         datau1 = vec.x.values[cube]
    #         datav1 = vec.y.values[cube]
    #     else:
    #         vectors = np.transpose(
    #             [vec.x.values[cube], vec.y.values[cube], vec.z.values[cube]])
    #         datau1 = np.inner(vectors, dir_vecs[1][1])
    #         datav1 = np.inner(vectors, dir_vecs[2][1])

    #     to_process["u_vect"] = datau1
    #     to_process["v_vect"] = datav1
    #     to_process["w_vect"] = np.sqrt(datau1**2+datav1**2)


        to_process.append(data)
        to_render.append({"mode": settings["mode"],
                                "params": params,
                                "unit": data.unit.units,
                                "name": data.name})
        # operations[data.name] = settings["operation"]


    # # Find parent container of object to plot
    # print(to_process)
    # key = list(to_process.keys())[0]
    parent = to_process[0].parent
    # parent = layers[0].parent



    # if scalar:
    #     holder = scalar.parent
    # elif image:
    #     holder = image.parent
    # elif contour:
    #     holder = contour.parent
    # elif vec:
    #     holder = vec.parent
    # elif stream:
    #     holder = stream.parent
    # else:
    #     print("Nothing to plot.")
    #     return

    # if holder.info["ndim"] < 2:
    #     print("plot_slice Error: Cannot plot slice from 1D data. Exiting...")
    #     return

    # # Possibility of updating the data from inside the plotting routines
    # try:
    #     update += 0
    #     holder.update_values(nout=update)
    # except TypeError:
    #     pass

    # Make it possible to call with only one size in the arguments
    if dy == 0.0:
        dy = dx


    # # Determine if interpolation is to be done in 2d or 3d
    # if (interpolation.count("3d") > 0) or (interpolation.count("3D") > 0):
    #     inter_3d = True
    #     # size_fact = np.sqrt(3.0)
    #     chars = [",", ":", ";", " ", "3d", "3D"]
    #     for ch in chars:
    #         interpolation = interpolation.replace(ch, "")
    #     if len(interpolation) == 0:
    #         interpolation = "linear"
    # else:
    # inter_3d = False
    # size_fact = 1.0
    # sqrt3 = np.sqrt(3.0)

    # # Get slice extent and direction vectors
    # if slice_direction is not None:
    #     [dx, dy, box, dir_vecs, origin] = slice_direction
    # else:
    #     dx, dy, box, dir_vecs, origin = get_slice_direction(
    #         holder, direction, dx, dy, origin=origin)

    dir_vecs, origin = get_slice_direction(direction=direction,
    	parent=parent, dx=0.5*(dx+dy), origin=origin)


    # # Define equation of a plane
    # a_plane = dir_vecs[0][1][0]
    # b_plane = dir_vecs[0][1][1]
    # c_plane = dir_vecs[0][1][2]
    # d_plane = -dir_vecs[0][1][0]*origin[0] - \
    #     dir_vecs[0][1][1]*origin[1]-dir_vecs[0][1][2]*origin[2]

    # # Distance to the plane
    # dist1 = (a_plane*parent["x"].values +
    #          b_plane*parent["y"].values +
    #          c_plane*parent["z"].values +
    #          d_plane) / np.sqrt(a_plane**2 + b_plane**2 + c_plane**2)
    # print('1', parent["xyz"] - origin)
    # print('2', (parent["xyz"] - origin) * dir_vecs[0][1])
    # print('3', np.sum((parent["xyz"] - origin) * dir_vecs[0][1]))
    xyz = parent["xyz"] - origin
    diagonal = parent["dx"] * np.sqrt(3.0) * 0.5

    dist1 = np.sum(xyz * dir_vecs[0][1], axis=1) / np.linalg.norm(dir_vecs[0][1])

    # Distance from center
    # dist2 = xyz.values - sqrt3*0.5*parent["dx"].values
    dist2 = xyz - diagonal
    # print(dist2)
    # dist2 = np.sqrt((parent["x"].values-origin[0])**2 +
    #                 (parent["y"].values-origin[1])**2 +
    #                 (parent["z"].values-origin[2])**2) - \
    #     sqrt3*0.5*parent["dx"].values

    # Select only the cells in contact with the slice., i.e. at a distance less than dx/2
    # cube = np.ravel(np.where(np.logical_and(np.abs(dist1) <= 0.5001*parent["dx"].values*sqrt3,
    #                                np.abs(dist2) <= max(dx, dy)*0.5*np.sqrt(2.0))))
    cube = np.ravel(np.where(np.logical_and(np.abs(dist1) <= 1.0001*diagonal,
                                   np.abs(dist2.values) <= max(dx, dy)*0.5*np.sqrt(2.0))))
    # print(cube)

    ncells = len(cube)

    # Project coordinates onto the plane by taking dot product with axes vectors
    # coords = np.transpose([parent["x"].values.take(cube)-origin[0],
    #                        parent["y"].values.take(cube)-origin[1],
    #                        parent["z"].values.take(cube)-origin[2]])
    coords = xyz._array[cube]
    # print("coords.shape", coords.shape)
    # print(coords._array)
    # print(dir_vecs)
    # print(dir_vecs[1][1])
    # print(dir_vecs[2][1])
    datax = np.inner(coords, dir_vecs[1][1])
    datay = np.inner(coords, dir_vecs[2][1])
    # datadx = sqrt3*0.5*parent["dx"].values.take(cube)
    datadx = diagonal[cube]

    # Define slice extent and resolution
    # xmin = max(-0.5*dx, box[0])
    # xmax = min(xmin+dx, box[1])
    # ymin = max(-0.5*dy, box[2])
    # ymax = min(ymin+dy, box[3])
    xmin = -0.5 * dx
    xmax = xmin + dx
    ymin = -0.5 * dy
    ymax = ymin + dy
    nx = resolution
    ny = resolution
    dpx = (xmax-xmin)/float(nx)
    dpy = (ymax-ymin)/float(ny)
    x = np.linspace(xmin+0.5*dpx, xmax-0.5*dpx, nx)
    y = np.linspace(ymin+0.5*dpy, ymax-0.5*dpy, ny)


    # to_render = {"scalar": None, "image": None, "contour": None,
    #              "u_vect": None, "v_vect": None, "w_vect": None,
    #              "u_strm": None, "v_strm": None, "w_strm": None}

    # to_process = {}

    # if scalar:
    #     to_process["scalar"] = scalar.values[cube]

    # if image:
    #     to_process["image"] = image.values[cube]

    # if contour:
    #     to_process["contour"] = contour.values[cube]

    # if vec:
    #     if holder.info["ndim"] < 3:
    #         datau1 = vec.x.values[cube]
    #         datav1 = vec.y.values[cube]
    #     else:
    #         vectors = np.transpose(
    #             [vec.x.values[cube], vec.y.values[cube], vec.z.values[cube]])
    #         datau1 = np.inner(vectors, dir_vecs[1][1])
    #         datav1 = np.inner(vectors, dir_vecs[2][1])

    #     to_process["u_vect"] = datau1
    #     to_process["v_vect"] = datav1
    #     to_process["w_vect"] = np.sqrt(datau1**2+datav1**2)

    # if stream:
    #     if holder.info["ndim"] < 3:
    #         datau2 = stream.x.values[cube]
    #         datav2 = stream.y.values[cube]
    #     else:
    #         streams = np.transpose(
    #             [stream.x.values[cube], stream.y.values[cube], stream.z.values[cube]])
    #         datau2 = np.inner(streams, dir_vecs[1][1])
    #         datav2 = np.inner(streams, dir_vecs[2][1])

    #     to_process["u_strm"] = datau2
    #     to_process["v_strm"] = datav2
    #     to_process["w_strm"] = np.sqrt(datau2**2+datav2**2)

    counts = np.zeros([ny, nx])
    for ind in range(len(to_process)):
        # to_render[data.name] = {"mode": settings["mode"]
        if to_render[ind]["mode"] in ["vec", "stream"]:
            if to_process[ind].ndim < 3:
                to_process[ind] = to_process[ind]._array[cube]
            else:
                # vectors = np.transpose(
                #     [vec.x.values[cube], vec.y.values[cube], vec.z.values[cube]])
                # u = np.inner(to_process[key]._array.take(cube, axis=0), dir_vecs[1][1])
                # v = np.inner(to_process[key]._array.take(cube, axis=0), dir_vecs[2][1])
                uv = np.inner(to_process[ind]._array.take(cube, axis=0), dir_vecs[1:3][1])
                # w = np.sqrt(u**2+v**2)
                w = None
                if "color" in to_render[ind]["params"]:
                    if not isinstance(to_render[ind]["params"]["color"], str):
                        w = to_render[ind]["params"]["color"]
                if w is None:
                    w = np.linalg.norm(uv, axis=1)
                to_process[ind] = np.concatenate((uv, w.reshape(ncells, 1)), axis=1)
            to_render[ind]["data"] = np.zeros([ny, nx, to_process[ind].shape[1]])
        else:
            to_process[ind] = to_process[ind].values.take(cube)
            to_render[ind]["data"] = np.zeros([ny, nx])

    datax -= xmin
    datay -= ymin
    istart = ((datax - datadx) / dpx).astype(np.int)
    iend = ((datax + datadx) / dpx + 1).astype(np.int)
    jstart = ((datay - datadx) / dpy).astype(np.int)
    jend = ((datay + datadx) / dpy + 1).astype(np.int)

    # print("=============")
    # print(istart)
    # print(type(istart))
    # print(istart[0])
    # print("=============")

    for i in range(len(istart)):
        i0 = istart[i]
        i1 = iend[i]
        j0 = jstart[i]
        j1 = jend[i]
        if i0 <= nx and j0 <= ny and i1 > 0 and j1 > 0:
            i0 = max(i0, 0)
            i1 = min(i1, nx)
            j0 = max(j0, 0)
            j1 = min(j1, ny)
            for ind in range(len(to_process)):
                # print(to_render[key]["data"].shape)
                to_render[ind]["data"][j0:j1, i0:i1] += to_process[ind][i]
            counts[j0:j1, i0:i1] += 1.0

    # Normalize by counts
    counts = np.ma.masked_where(counts == 0.0, counts)
    for ind in range(len(to_process)):
        if to_render[ind]["data"].ndim > counts.ndim:
            to_render[ind]["data"] /= counts.reshape(ny, nx, 1)
        else:
            to_render[ind]["data"] /= counts

    # Render the map

    figure = render(x=x, y=y, data=to_render, ax=ax)

    figure["ax"].set_xlabel(parent["xyz"].x.label)
    figure["ax"].set_ylabel(parent["xyz"].y.label)

    return Plot(x=x, y=y, layers=to_render, fig=figure["fig"], ax=figure["ax"])


    # figure = render(scalar=scalar, image=image, contour=contour, vec=vec, stream=stream, x=x, y=y, z_scal=to_render["scalar"],
    #                z_imag=to_render["image"], z_cont=to_render["contour"], u_vect=to_render["u_vect"], v_vect=to_render["v_vect"], w_vect=to_render["w_vect"], u_strm=to_render["u_strm"],
    #                v_strm=to_render["v_strm"], w_strm=to_render["w_strm"], xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fname=fname,
    #                axes=axes, title=title, sinks=sinks, new_window=new_window, clear=clear, block=block,
    #                resolution=resolution, thePlane=[
    #                    a_plane, b_plane, c_plane, d_plane],
    #                origin=origin, dir_vecs=dir_vecs, scalar_args=scalar_args, image_args=image_args,
    #                contour_args=contour_args, vec_args=vec_args, stream_args=stream_args, outline=outline,
    #                outline_args=outline_args, sink_args=sink_args, x_raw=datax, y_raw=datay, holder=holder)

    # return OsyrisPlot(x=x, y=y, layers=to_render, fig=figure["fig"], ax=figure["ax"])

    # # if copy:
    # #     return x, y, to_render
    # # else:
    # #     return
