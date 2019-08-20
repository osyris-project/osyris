# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2019 Osyris contributors (https://github.com/nvaytet/osyris)
# @author Neil Vaytet

import numpy as np
from scipy.interpolate import griddata
from .plot import get_slice_direction, render_map


def plot_column_density(scalar=False, image=False, contour=False, vec=False, stream=False,
                        direction="z", dx=0.0, dy=0.0, dz=0.0, fname=None, axes=None, title=None,
                        origin=[0, 0, 0], resolution=128, sinks=True, summed=True, copy=False,
                        new_window=False, update=None, clear=True, plot=True, block=False, nz=0,
                        interpolation="linear", verbose=False, outline=False, outline_args={},
                        scalar_args={}, image_args={}, contour_args={}, vec_args={}, stream_args={},
                        sink_args={}, lmax=0):
    """
    Plot a column density through the data cube. The arguments are:
    - scalar     : the scalar field to be plotted, e.g. mydata.density
    - image      : the scalar field to be plotted with an image
    - contour    : the scalar field to be plotted with contours
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
        print("plot_column_density Error: Cannot plot slice from 1D data. Exiting...")
        return

    # Possibility of updating the data from inside the plotting routines
    try:
        update += 0
        holder.update_values(nout=update)
    except TypeError:
        pass

    # Get direction vectors once and for all for the column_density.
    # This should be computed here and not inside the plot_slice routine as the origin
    # changes along the z direction inside the loop below.
    dx, dy, box, dir_vecs, origin = get_slice_direction(
        holder, direction, dx, dy, origin=origin)

    # Compute domain dimension for integration
    if dz == 0.0:
        dz = max(dx, dy)
    zmin = -0.5*dz
    zmax = 0.5*dz
    nx = resolution
    ny = resolution
    if nz == 0:
        nz = resolution
    dpz = (zmax-zmin)/float(nz)
    z = np.linspace(zmin+0.5*dpz, zmax-0.5*dpz, nz)

    # We now create empty data arrays that will be filled by the cell data
    z_scal = z_imag = z_cont = u_vect = v_vect = w_vect = u_strm = v_strm = w_strm = 0
    if scalar:
        z_scal = np.zeros([ny, nx])
    if image:
        z_imag = np.zeros([ny, nx])
    if contour:
        z_cont = np.zeros([ny, nx])
    if vec:
        u_vect = np.zeros([ny, nx])
        v_vect = np.zeros([ny, nx])
        w_vect = np.zeros([ny, nx])
    if stream:
        u_strm = np.zeros([ny, nx])
        v_strm = np.zeros([ny, nx])
        w_strm = np.zeros([ny, nx])

    # Define equation of a plane
    a_plane = dir_vecs[0][1][0]
    b_plane = dir_vecs[0][1][1]
    c_plane = dir_vecs[0][1][2]
    d_plane = -dir_vecs[0][1][0]*origin[0] - \
        dir_vecs[0][1][1]*origin[1]-dir_vecs[0][1][2]*origin[2]

    iprog = 1
    istep = 10

    # Begin loop over vertical direction, calling plot_slice with plot=False and copy=True
    # to retrieve the slice data and create a sum
    for iz in range(nz):

        # Print progress
        if verbose:
            percentage = int(float(iz)*100.0/float(nz))
            if percentage >= iprog*istep:
                print("Column density: %3i%% done" % percentage)
                iprog += 1

        [x, y, z_scal_slice, z_imag_slice, z_cont_slice, u_vect_slice, v_vect_slice,
            w_vect_slice, u_strm_slice, v_strm_slice, w_strm_slice] = \
            plot_slice(scalar=scalar, image=image, contour=contour, vec=vec, stream=stream,
                       direction=direction, dx=dx, dy=dy, sinks=sinks, copy=True, resolution=resolution,
                       origin=[origin[0], origin[1], origin[2]+z[iz]], plot=False, interpolation=interpolation, lmax=lmax,
                       slice_direction=[dx, dy, box, dir_vecs, [origin[0], origin[1], origin[2]+z[iz]]])

        # Increment the sum
        if scalar:
            z_scal += z_scal_slice
        if image:
            z_imag += z_imag_slice
        if contour:
            z_cont += z_cont_slice
        if vec:
            u_vect += u_vect_slice
            v_vect += v_vect_slice
            w_vect += w_vect_slice
        if stream:
            u_strm += u_strm_slice
            v_strm += v_strm_slice
            w_strm += w_strm_slice

    # If summed=True, this is a real column density.
    # Else, only the average of the quantity is requested
    if summed:
        column = (zmax-zmin)*conf.constants[holder.info["scale"]]/float(nz)
    else:
        column = 1.0/float(nz)

    if scalar:
        z_scal *= column
    if image:
        z_imag *= column
    if contour:
        z_cont *= column
    if vec:
        u_vect *= column
        v_vect *= column
        w_vect *= column
    if stream:
        u_strm *= column
        v_strm *= column
        w_strm *= column

    # Render the map
    if plot:
        dpx = x[1] - x[0]
        dpy = y[1] - y[0]
        xmin = x[0] - 0.5*dpx
        xmax = x[-1] + 0.5*dpx
        ymin = y[0] - 0.5*dpy
        ymax = y[-1] + 0.5*dpy
        render_map(scalar=scalar, image=image, contour=contour, vec=vec, stream=stream, x=x, y=y, z_scal=z_scal,
                   z_imag=z_imag, z_cont=z_cont, u_vect=u_vect, v_vect=v_vect, w_vect=w_vect, u_strm=u_strm,
                   v_strm=v_strm, w_strm=w_strm, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fname=fname, dz=dz,
                   axes=axes, title=title, sinks=sinks, new_window=new_window, clear=clear, block=block,
                   resolution=resolution, thePlane=[
                       a_plane, b_plane, c_plane, d_plane],
                   origin=origin, dir_vecs=dir_vecs, scalar_args=scalar_args, image_args=image_args,
                   contour_args=contour_args, vec_args=vec_args, stream_args=stream_args, outline=outline,
                   outline_args=outline_args, sink_args=sink_args, holder=holder)

    if copy:
        return x, y, z_scal, z_imag, z_cont, u_vect, v_vect, w_vect, u_strm, v_strm, w_strm
    else:
        return
