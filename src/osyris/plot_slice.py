# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2019 Osyris contributors (https://github.com/nvaytet/osyris)
# @author Neil Vaytet

import numpy as np
from scipy.interpolate import griddata
from .plot import get_slice_direction, render_map


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
        size_fact = np.sqrt(3.0)
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

    # Try to automatically determine lmax to speedup process
    if lmax == 0:
        dxlmax = holder.info["boxsize_scaled"] * \
            (0.5**holder.info["levelmax_active"])
        target = min(dx, dy)/float(resolution)
        lmax = round(np.log((min(dx, dy)/float(resolution)) /
                            holder.info["boxsize_scaled"])/(np.log(0.5)))
    subset = np.where(np.logical_or(np.logical_and(holder.get("level", only_leafs=False) < lmax,
                                                   holder.get("leaf", only_leafs=False) > 0.0), holder.get("level", only_leafs=False) == lmax))

    # Define equation of a plane
    a_plane = dir_vecs[0][1][0]
    b_plane = dir_vecs[0][1][1]
    c_plane = dir_vecs[0][1][2]
    d_plane = -dir_vecs[0][1][0]*origin[0] - \
        dir_vecs[0][1][1]*origin[1]-dir_vecs[0][1][2]*origin[2]

    # Distance to the plane
    dist1 = (a_plane*holder.get("x", only_leafs=False)[subset] +
             b_plane*holder.get("y", only_leafs=False)[subset] +
             c_plane*holder.get("z", only_leafs=False)[subset] +
             d_plane) / np.sqrt(a_plane**2 + b_plane**2 + c_plane**2)
    # Distance from center
    dist2 = np.sqrt((holder.get("x", only_leafs=False)[subset]-origin[0])**2 +
                    (holder.get("y", only_leafs=False)[subset]-origin[1])**2 +
                    (holder.get("z", only_leafs=False)[subset]-origin[2])**2) - \
        np.sqrt(3.0)*0.5*holder.get("dx", only_leafs=False)[subset]

    # Select only the cells in contact with the slice., i.e. at a distance less than dx/2
    cube = np.where(np.logical_and(np.abs(dist1) <= 0.5*holder.get("dx", only_leafs=False)[subset]*size_fact,
                                   np.abs(dist2) <= max(dx, dy)*0.5*np.sqrt(2.0)))

    # Project coordinates onto the plane by taking dot product with axes vectors
    coords = np.transpose([holder.get("x", only_leafs=False)[subset][cube]-origin[0],
                           holder.get("y", only_leafs=False)[
        subset][cube]-origin[1],
        holder.get("z", only_leafs=False)[subset][cube]-origin[2]])
    datax = np.inner(coords, dir_vecs[1][1])
    datay = np.inner(coords, dir_vecs[2][1])

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
    grid_x, grid_y = np.meshgrid(x, y)

    # Doing different things depending on whether interpolation is 2D or 3D
    if inter_3d:
        dataz = np.inner(coords, dir_vecs[0][1])
        grid_z = np.zeros([nx, ny])
        points = np.transpose([datax, datay, dataz])
        grids = (grid_x, grid_y, grid_z)
    else:
        points = np.transpose([datax, datay])
        grids = (grid_x, grid_y)

    # Use scipy interpolation function to make image
    z_scal = z_imag = z_cont = u_vect = v_vect = w_vect = u_strm = v_strm = w_strm = 0
    if scalar:
        z_scal = griddata(
            points, scalar.values[subset][cube], grids, method=interpolation)
    if image:
        z_imag = griddata(
            points, image.values[subset][cube], grids, method=interpolation)
    if contour:
        z_cont = griddata(
            points, contour.values[subset][cube], grids, method=interpolation)
    if vec:
        if holder.info["ndim"] < 3:
            datau1 = vec.x.values[subset][cube]
            datav1 = vec.y.values[subset][cube]
        else:
            vectors = np.transpose(
                [vec.x.values[subset][cube], vec.y.values[subset][cube], vec.z.values[subset][cube]])
            datau1 = np.inner(vectors, dir_vecs[1][1])
            datav1 = np.inner(vectors, dir_vecs[2][1])
        u_vect = griddata(points, datau1, grids, method=interpolation)
        v_vect = griddata(points, datav1, grids, method=interpolation)
        if "colors" in vec_args.keys():
            w_vect = griddata(
                points, vec_args["colors"].values[subset][cube], grids, method=interpolation)
        else:
            w_vect = griddata(points, np.sqrt(
                datau1**2+datav1**2), grids, method=interpolation)
    if stream:
        if holder.info["ndim"] < 3:
            datau2 = stream.x.values[subset][cube]
            datav2 = stream.y.values[subset][cube]
        else:
            streams = np.transpose(
                [stream.x.values[subset][cube], stream.y.values[subset][cube], stream.z.values[subset][cube]])
            datau2 = np.inner(streams, dir_vecs[1][1])
            datav2 = np.inner(streams, dir_vecs[2][1])
        u_strm = griddata(points, datau2, grids, method=interpolation)
        v_strm = griddata(points, datav2, grids, method=interpolation)
        w_strm = griddata(points, np.sqrt(datau2**2+datav2**2),
                          grids, method=interpolation)

    # Render the map
    if plot:
        render_map(scalar=scalar, image=image, contour=contour, vec=vec, stream=stream, x=x, y=y, z_scal=z_scal,
                   z_imag=z_imag, z_cont=z_cont, u_vect=u_vect, v_vect=v_vect, w_vect=w_vect, u_strm=u_strm,
                   v_strm=v_strm, w_strm=w_strm, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fname=fname,
                   axes=axes, title=title, sinks=sinks, new_window=new_window, clear=clear, block=block,
                   resolution=resolution, thePlane=[
                       a_plane, b_plane, c_plane, d_plane],
                   origin=origin, dir_vecs=dir_vecs, scalar_args=scalar_args, image_args=image_args,
                   contour_args=contour_args, vec_args=vec_args, stream_args=stream_args, outline=outline,
                   outline_args=outline_args, sink_args=sink_args, x_raw=datax, y_raw=datay, holder=holder)

    if copy:
        return x, y, z_scal, z_imag, z_cont, u_vect, v_vect, w_vect, u_strm, v_strm, w_strm
    else:
        return
