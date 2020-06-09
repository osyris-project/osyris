# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2019 Osyris contributors (https://github.com/nvaytet/osyris)
# @author Neil Vaytet

import numpy as np
from scipy.interpolate import griddata
from scipy.stats import binned_statistic_2d
from .plot import get_slice_direction, render_map

from timeit import default_timer as timer


def plot_slice(scalar=False, image=False, contour=False, vec=False, stream=False, axes=None,
               direction="z", dx=0.0, dy=0.0, fname=None, title=None, sinks=True, copy=False,
               origin=[0, 0, 0], resolution=128, new_window=False, update=None,
               clear=True, plot=True, block=False, interpolation="linear", sink_args={},
               scalar_args={}, image_args={}, contour_args={}, vec_args={}, stream_args={},
               outline=False, outline_args={}, lmax=0, slice_direction=None):
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

    # Find parent container of object to plot
    if scalar:
        holder = scalar.parent
    elif image:
        holder = image.parent
    elif contour:
        holder = contour.parent
    elif vec:
        holder = vec.parent
    elif stream:
        holder = stream.parent
    else:
        print("Nothing to plot.")
        return

    if holder.info["ndim"] < 2:
        print("plot_slice Error: Cannot plot slice from 1D data. Exiting...")
        return

    # Possibility of updating the data from inside the plotting routines
    try:
        update += 0
        holder.update_values(nout=update)
    except TypeError:
        pass

    # Determine if interpolation is to be done in 2d or 3d
    if (interpolation.count("3d") > 0) or (interpolation.count("3D") > 0):
        inter_3d = True
        # size_fact = np.sqrt(3.0)
        chars = [",", ":", ";", " ", "3d", "3D"]
        for ch in chars:
            interpolation = interpolation.replace(ch, "")
        if len(interpolation) == 0:
            interpolation = "linear"
    else:
        inter_3d = False
        size_fact = 1.0

    # Get slice extent and direction vectors
    if slice_direction is not None:
        [dx, dy, box, dir_vecs, origin] = slice_direction
    else:
        dx, dy, box, dir_vecs, origin = get_slice_direction(
            holder, direction, dx, dy, origin=origin)

    # # Try to automatically determine lmax to speedup process
    # if lmax == 0:
    #     # dxlmax = holder.info["boxsize_scaled"] * \
    #     #     (0.5**holder.info["levelmax_active"])
    #     # target = min(dx, dy)/float(resolution)
    #     lmax = int(round(np.log((min(dx, dy)/float(resolution)) /
    #                             holder.info["boxsize_scaled"])/(np.log(0.5)))) + 2
    #     # subset = np.where(holder.get("level") <= lmax)
    #     # print(lmax)
    #     # nmax = holder.info["level_indices"][int(lmax)]
    #     # print(nmax)
    # # else:
    # #     lmax = holder.info["levelmax_active"]
    # # print(lmax, int(lmax))
    # if lmax >= int(holder.info["levelmax_active"]):
    #     nmax = holder.info["ncells"] + 1
    # else:
    #     nmax = holder.info["level_indices"][int(lmax)]
    # print(lmax, nmax)
    sqrt3 = np.sqrt(3.0)

    # Define equation of a plane
    a_plane = dir_vecs[0][1][0]
    b_plane = dir_vecs[0][1][1]
    c_plane = dir_vecs[0][1][2]
    d_plane = -dir_vecs[0][1][0]*origin[0] - \
        dir_vecs[0][1][1]*origin[1]-dir_vecs[0][1][2]*origin[2]

    # Distance to the plane
    dist1 = (a_plane*holder.get("x") +
             b_plane*holder.get("y") +
             c_plane*holder.get("z") +
             d_plane) / np.sqrt(a_plane**2 + b_plane**2 + c_plane**2)
    # Distance from center
    dist2 = np.sqrt((holder.get("x")-origin[0])**2 +
                    (holder.get("y")-origin[1])**2 +
                    (holder.get("z")-origin[2])**2) - \
        sqrt3*0.5*holder.get("dx")

    # Select only the cells in contact with the slice., i.e. at a distance less than dx/2
    # cube = np.where(np.logical_and(np.abs(dist1) <= 0.5*holder.get("dx")*size_fact,
    #                                np.abs(dist2) <= max(dx, dy)*0.5*np.sqrt(2.0)))
    cube = np.where(np.logical_and(np.abs(dist1) <= 0.5001*holder.get("dx")*sqrt3,
                                   np.abs(dist2) <= max(dx, dy)*0.5*np.sqrt(2.0)))
    # cube = np.where(np.abs(dist1) <= 0.5001*holder.get("dx")*sqrt3)

    # Project coordinates onto the plane by taking dot product with axes vectors
    coords = np.transpose([holder.get("x")[cube]-origin[0],
                           holder.get("y")[cube]-origin[1],
                           holder.get("z")[cube]-origin[2]])
    datax = np.inner(coords, dir_vecs[1][1])
    datay = np.inner(coords, dir_vecs[2][1])
    print(len(holder.get("x")[cube]))
    datadx = sqrt3*0.5*holder.get("dx")[cube]

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
    # grid_x, grid_y = np.meshgrid(x, y)

    # # xe = np.linspace(xmin, xmax, nx*2 + 1)
    # # ye = np.linspace(ymin, ymax, ny*2 + 1)
    # # xe = np.linspace(xmin, xmax, nx + 1)
    # # ye = np.linspace(ymin, ymax, ny + 1)
    # xe = np.linspace(xmin-0.1*dx, xmax+0.1*dx, 1.5*nx)
    # ye = np.linspace(ymin-0.1*dy, ymax+0.1*dy, 1.5*ny)

    # xc = 0.5 * (xe[:-1] + xe[1:])
    # yc = 0.5 * (ye[:-1] + ye[1:])
    # grid_x2, grid_y2 = np.meshgrid(xc, yc)

    # # Doing different things depending on whether interpolation is 2D or 3D
    # if inter_3d:
    #     dataz = np.inner(coords, dir_vecs[0][1])
    #     grid_z = np.zeros([nx, ny])
    #     points = np.transpose([datax, datay, dataz])
    #     grids = (grid_x, grid_y, grid_z)
    # else:
    #     points = np.transpose([datax, datay])
    #     # points = np.transpose([grid_x2.flatten(), grid_y2.flatten()])
    #     grids = (grid_x, grid_y)

    # # nums, yedg, xedg = np.histogram2d(datay, datax, bins=(ye, xe))
    # # nums = nums.flatten()
    # # subs = np.where(nums > 0)
    # # points = np.transpose([grid_x2.flatten()[subs], grid_y2.flatten()[subs]])

    # # Use scipy interpolation function to make image
    # z_scal = z_imag = z_cont = u_vect = v_vect = w_vect = u_strm = v_strm = w_strm = 0

    to_render = {"scalar": None, "image": None, "contour": None,
                 "u_vect": None, "v_vect": None, "w_vect": None,
                 "u_strm": None, "v_strm": None, "w_strm": None}


    to_process = {}



    if scalar:
        to_process["scalar"] = scalar.values[cube]
        # # vals, xedg, yedg = np.histogram2d(datax, datay, bins=(xe, ye), weights=scalar.values[cube])
        # # nums, xedg, yedg = np.histogram2d(datax, datay, bins=(xe, ye))
        # vals, yedg, xedg = np.histogram2d(datay, datax, bins=(ye, xe), weights=scalar.values[cube])
        # print(len(points), len(vals.flatten()))
        # # vals = vals.flatten()
        # # subs = np.where(nums > 0)
        # print(len(nums[subs]))
        # z_scal = griddata(
        #     points, vals.flatten()[subs]/nums[subs], grids, method=interpolation)

        # # z_scal = griddata(
        # #     np.transpose([datax, datay]), scalar.values[cube], grids, method=interpolation)
        

        # axes.scatter(points[:, 0], points[:, 1])
        # return

    if image:
        # z_imag = griddata(
        #     points, image.values[cube], grids, method=interpolation)
        to_process["image"] = image.values[cube]
    if contour:
        to_process["contour"] = contour.values[cube]
        # z_cont = griddata(
        #     points, contour.values[cube], grids, method=interpolation)
    if vec:
        if holder.info["ndim"] < 3:
            datau1 = vec.x.values[cube]
            datav1 = vec.y.values[cube]
        else:
            vectors = np.transpose(
                [vec.x.values[cube], vec.y.values[cube], vec.z.values[cube]])
            datau1 = np.inner(vectors, dir_vecs[1][1])
            datav1 = np.inner(vectors, dir_vecs[2][1])

        to_process["u_vect"] = datau1
        to_process["v_vect"] = datav1
        to_process["w_vect"] = np.sqrt(datau1**2+datav1**2)

        # # u_vect = griddata(points, datau1, grids, method=interpolation)
        # # v_vect = griddata(points, datav1, grids, method=interpolation)
        # # if "colors" in vec_args.keys():
        # #     w_vect = griddata(
        # #         points, vec_args["colors"].values[cube], grids, method=interpolation)
        # # else:
        # #     w_vect = griddata(points, np.sqrt(
        # #         datau1**2+datav1**2), grids, method=interpolation)

        # ### should use binned_statistic !!!####


        # vals_u, yedg, xedg = np.histogram2d(datay, datax, bins=(ye, xe), weights=datau1)
        # vals_v, yedg, xedg = np.histogram2d(datay, datax, bins=(ye, xe), weights=datav1)
        # vals_w, yedg, xedg = np.histogram2d(datay, datax, bins=(ye, xe), weights=np.sqrt(datau1**2+datav1**2))
        # # vals_u = vals_u.flatten()
        # # vals_v = vals_v.flatten()
        # # subs = np.where(nums > 0)
        # print(len(nums[subs]))
        # u_vect = griddata(
        #     points, vals_u.flatten()[subs]/nums[subs], grids, method=interpolation)
        # v_vect = griddata(
        #     points, vals_v.flatten()[subs]/nums[subs], grids, method=interpolation)
        # w_vect = griddata(
        #     points, vals_w.flatten()[subs]/nums[subs], grids, method=interpolation)



    if stream:
        if holder.info["ndim"] < 3:
            datau2 = stream.x.values[cube]
            datav2 = stream.y.values[cube]
        else:
            streams = np.transpose(
                [stream.x.values[cube], stream.y.values[cube], stream.z.values[cube]])
            datau2 = np.inner(streams, dir_vecs[1][1])
            datav2 = np.inner(streams, dir_vecs[2][1])

        to_process["u_strm"] = datau2
        to_process["v_strm"] = datav2
        to_process["w_strm"] = np.sqrt(datau2**2+datav2**2)

        # u_strm = griddata(points, datau2, grids, method=interpolation)
        # v_strm = griddata(points, datav2, grids, method=interpolation)
        # w_strm = griddata(points, np.sqrt(datau2**2+datav2**2),
        #                   grids, method=interpolation)
    # print(list(to_process.values()))


    # results, y_edges, x_edges, bin_number = binned_statistic_2d(x=datay, y=datax, values=list(to_process.values()), statistic='mean', bins=[ye, xe])
    # # print(results)
    # # print(len(results))

    # # nums, yedg, xedg = np.histogram2d(datay, datax, bins=(ye, xe))
    # # nums = results[0].flatten()
    # subs = np.where(np.isfinite(results[0].flatten()))
    # points = np.transpose([grid_x2.flatten()[subs], grid_y2.flatten()[subs]])


    # start = timer()
    # for i, key in enumerate(to_process.keys()):
        # to_render[key] = griddata(
        #     points, results[i].flatten()[subs], grids, method=interpolation)
    start = timer()
    counts = np.zeros([ny, nx])
    for key in to_process.keys():
        to_render[key] = np.zeros([ny, nx])
        # counts[key] = np.zeros([ny, nx])
    print("init", timer() - start)
    start = timer()

    datax -= xmin
    datay -= ymin
    istart = ((datax - datadx) / dpx).astype(np.int)
    iend = ((datax + datadx) / dpx + 1).astype(np.int)
    jstart = ((datay - datadx) / dpy).astype(np.int)
    jend = ((datay + datadx) / dpy + 1).astype(np.int)
    print(istart, type(istart))
    # to_render["scalar"]
    print("ijstart", timer() - start)
    start = timer()

    # sel = np.where(np.logical_and(istart >= 0, np.logical_and(jstart >= 0,
    #     np.logical_and(iend <= nx+1, jend <= ny+1))))

    # # np.add.at(to_render["scalar"],([jstart[sel]:jend[sel]],[istart[sel]:iend[sel]]), 1)
    # istart = istart[sel]
    # iend = iend[sel]
    # jstart = jstart[sel]
    # jend = jend[sel]
    # print("sel", timer() - start)
    # start = timer()
    # # n = len(datax[sel])
    # # keys = to_process.keys()
    # # # key = "scalar"
    # # print('n', n)
    # # print(istart)
    # # # to_render = np.zeros([ny, nx])
    # for key in to_process.keys():
    #     to_process[key] = to_process[key][sel]

    for i in range(len(istart)):
        # i0 = int((datax[i] - datadx[i] - xmin) / dpx)
        # i1 = int((datax[i] + datadx[i] - xmin) / dpx) + 1
        # j0 = int((datay[i] - datadx[i] - ymin) / dpy)
        # j1 = int((datay[i] + datadx[i] - ymin) / dpy) + 1
        i0 = istart[i]
        i1 = iend[i]
        j0 = jstart[i]
        j1 = jend[i]
        if i0 <= nx and j0 <= ny and i1 > 0 and j1 > 0:
            i0 = max(i0, 0)
            i1 = min(i1, nx)
            j0 = max(j0, 0)
            j1 = min(j1, ny)
            for key in to_process.keys():
                to_render[key][j0:j1, i0:i1] += to_process[key][i]
            counts[j0:j1, i0:i1] += 1.0
        # to_render[j0:j1, i0:i1] += 1.0
        # # print(key, timer() - start)
    print("loop1", timer() - start)
    start = timer()

    # Normalize by counts
    for key in to_process.keys():
        to_render[key] /= counts

    # Render the map
    if plot:
        render_map(scalar=scalar, image=image, contour=contour, vec=vec, stream=stream, x=x, y=y, z_scal=to_render["scalar"],
                   z_imag=to_render["image"], z_cont=to_render["contour"], u_vect=to_render["u_vect"], v_vect=to_render["v_vect"], w_vect=to_render["w_vect"], u_strm=to_render["u_strm"],
                   v_strm=to_render["v_strm"], w_strm=to_render["w_strm"], xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fname=fname,
                   axes=axes, title=title, sinks=sinks, new_window=new_window, clear=clear, block=block,
                   resolution=resolution, thePlane=[
                       a_plane, b_plane, c_plane, d_plane],
                   origin=origin, dir_vecs=dir_vecs, scalar_args=scalar_args, image_args=image_args,
                   contour_args=contour_args, vec_args=vec_args, stream_args=stream_args, outline=outline,
                   outline_args=outline_args, sink_args=sink_args, x_raw=datax, y_raw=datay, holder=holder)

    if copy:
        return x, y, to_render
    else:
        return
