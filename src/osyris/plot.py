# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2019 Osyris contributors (https://github.com/nvaytet/osyris)
# @author Neil Vaytet

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections
from matplotlib.colors import LogNorm, Normalize
from scipy.interpolate import griddata
from . import config as conf


def plot_histogram(var_x,
                   var_y,
                   axes=None,
                   block=False,
                   cbar=True,
                   clear=True,
                   contour=False,
                   contour_args={},
                   copy=False,
                   equal_axes=False,
                   fname=None,
                   image=False,
                   image_args={},
                   new_window=False,
                   only_leafs=True,
                   outline=False,
                   outline_args={},
                   plot=True,
                   resolution=256,
                   scalar=False,
                   scalar_args={},
                   scatter=False,
                   scatter_args={},
                   summed=False,
                   title=None,
                   update=None,
                   xmin=None,
                   xmax=None,
                   ymin=None,
                   ymax=None):
    """
    Plot a 2D histogram with two variables as input.
    This is used for instance to plot the temperature as a function of density
    for every cell in the mesh. You can supply a third variable for the
    colours. By default, if no variable for colour is supplied for either
    scalar, scatter, contour or image, the colour is used for the number of
    cells in each pixels.

    :param OsyrisField var_x: The variable for the x axis, e.g. mydata.log_rho.
    :param OsyrisField var_y: The variable for the x axis, e.g. mydata.log_T.
    :param MPL_axes axes: If specified, the data is plotted on the specified
        axes. The default is `None`.
    :param bool block: If ``True``, plotting the figure holds the terminal
        console until it is closed. The default is `False`.
    :param bool cbar: Make colorbar next to the plot if ``True``. The default
        is ``True``.
    :param bool clear: If ``True`` the figure is cleared before plotting. The
        default is `True`.
    :param OsyrisField contour: An osyris data field containing the variable
        to be represented using line contours. The default is ``False``.
    :param dict contour_args: A dictionary to hold additional matplotlib
        keyword arguments to be passed to the ``contour`` function. The default
        is empty.
    :param bool copy: If `True`, the function returns the data that was used to
        plot the histogram. Default is ``False``.

    """

    # :param equal_axes=False,
    # :param fname=None,
    # :param image=False,
    # :param image_args={},
    # :param new_window=False,
    # :param only_leafs=True,
    # :param outline=False,
    # :param outline_args={},
    # :param plot=True,
    # :param resolution=256,
    # :param scalar=False,
    # :param scalar_args={},
    # :param scatter=False,
    # :param scatter_args={},
    # :param summed=False,
    # :param title=None,
    # :param update=None,
    # :param xmin=None,
    # :param xmax=None,
    # :param ymin=None,
    # :param ymax=None):

    # :param str sender: The person sending the message
    # :param str recipient: The recipient of the message
    # :param str message_body: The body of the message
    # :param priority: The priority of the message, can be a number 1-5
    # :type priority: integer or None
    # :return: the message id
    # :rtype: int
    # :raises ValueError: if the message_body exceeds 160 characters
    # :raises TypeError: if the message_body is not a basestring

    # Parameters
    # ----------
    # x : type
    #     Description of parameter `x`.

    # var_x : osyrisField
    #     An osyris data field containing the variable for the x axis,
    #     e.g. mydata.log_rho.

    # var_y`: (*osyrisField*) An osyris data field containing the variable for the y
    # axis, e.g. mydata.log_T. There is no default.

    # scalar`: (*osyrisField*) An osyris data field containing the variable for the
    # colors to be used in the filled contours. Default is `False`.

    # image`: (*osyrisField*) An osyris data field containing the variable for the
    # colors to be used by the imshow function. Default is `False`.

    # contour`: (*osyrisField*) An osyris data field containing the variable for the
    # colors to be used in the filled contours. Default is `False`.

    # scatter`: (*osyrisField*) An osyris data field containing the variable for the
    # colors to be used in the scatter function. Default is `False`. **Note:**

    # fname`: (*string*) If specified, the figure is saved to file. Default is `None`.

    # axes`: (*MPL axes*) If specified, the data is plotted on the specified axes (see
    # demo). Default is `None`.

    # resolution`: (*integer*) The data is binned in a 2D matrix of size `resolution`.
    # Default is 256.

    # copy`: (*logical*) If `True`, the function returns the data that was used to plot
    # the histogram. Default is `False`.

    # title`: (*string*) The title for the figure/axes. Default is `None`. If no title
    # is specified, the simulation time will be used as a title. If you want no title on
    # the figure, use `title=""`.

    # xmin`: (*float*) Lower limit for the x axis. Default is `None`, implying an
    # automatic search for the lowest x value in the supplied var_x data.

    # xmax`: (*float*) Upper limit for the x axis. Default is `None`, implying an
    # automatic search for the highest x value in the supplied var_x data.

    # ymin`: (*float*) Lower limit for the y axis. Default is `None`, implying an
    # automatic search for the lowest y value in the supplied var_x data.

    # ymax`: (*float*) Upper limit for the y axis. Default is `None`, implying an
    # automatic search for the highest y value in the supplied var_x data.

    # plot`: (*logical*) Do not plot the data on a figure if `False`. Default is
    # `True`.

    # new_window`: (*logical*) Create a new plotting window if `True`. Default is
    # `False`.

    # update`: (*integer*) Update the values of the current snapshot by reading a
    # new simulation output, specified by the output number `update=new_output_number`.
    # This saves the user from having to first run `mydata.update_values(new_output_number)`
    # before calling `plot_histogram` again, it can all be done in one line. Default
    # is `None`.

    # cbar`: (*logical*) Make colorbar next to the plot if `True`. Default is `True`.

    # outline`: (*logical*) Print a contour as an outline around all the data points
    # in the histogram. Default is `False`.

    # summed`: (*logical*) If `True`, the data are summed in the histogram pixels,
    # rather than averaged. If no variable for colour is supplied for either
    # scalar, scatter, contour or image, the colour is used for the number of cells in
    # each pixels and `summed=True` is used. Otherwise, default is `False`.

    # clear`: (*logical*) If `True` the figure is cleared. Default is `True`.

    # block`: (*logical*) If `True` plotting the figure holds the terminal console
    # until it is closed. Default is `False`.

    # equal_axes`: (*logical*) If `True` the aspect ratio is conserved between x and y
    # axes. Default is `False`.

    # scalar_args`: (*dict*) A python dictionary to hold additional matplotlib arguments
    # to be passed to the `contourf` function. Default is empty.

    # image_args`: (*dict*) A python dictionary to hold additional matplotlib arguments
    # to be passed to the `imshow` function. Default is empty.

    # contour_args`: (*dict*) A python dictionary to hold additional matplotlib arguments
    # to be passed to the `contour` function. Default is empty.

    # scatter_args`: (*dict*) A python dictionary to hold additional matplotlib arguments
    # to be passed to the `scatter` function. Default is empty.

    # outline_args`: (*dict*) A python dictionary to hold additional matplotlib arguments
    # to be passed to the `contour` function. Default is empty.


    # Find parent container of object to plot
    holder = var_x.parent

    # Possibility of updating the data from inside the plotting routines
    try:
        update += 0
        holder.update_values(nout=update)
    except TypeError:
        pass

    # Parameters
    nx = resolution+1
    ny = resolution+1

    # Get the data values and units
    datax = holder.get(var_x.name, only_leafs=only_leafs)
    datay = holder.get(var_y.name, only_leafs=only_leafs)
    xlabel = var_x.label+" ["+var_x.unit+"]"
    ylabel = var_y.label+" ["+var_y.unit+"]"
    default_var = "histo_cell_density"

    # Define plotting range
    autoxmin = False
    autoxmax = False
    autoymin = False
    autoymax = False
    try:
        xmin += 0
    except TypeError:
        datax[np.isneginf(datax)] = np.inf
        xmin = np.nanmin(datax)
        autoxmin = True
    try:
        xmax += 0
    except TypeError:
        datax[np.isinf(datax)] = np.NINF
        xmax = np.nanmax(datax)
        autoxmax = True
    try:
        ymin += 0
    except TypeError:
        datay[np.isneginf(datay)] = np.inf
        ymin = np.nanmin(datay)
        autoymin = True
    try:
        ymax += 0
    except TypeError:
        datay[np.isinf(datay)] = np.NINF
        ymax = np.nanmax(datay)
        autoymax = True

    # Protect against empty plots if xmin==xmax or ymin==ymax
    if xmin == xmax:
        if xmin == 0.0:
            xmin = -0.1
            xmax = 0.1
        else:
            xmin = xmin - 0.05*abs(xmin)
            xmax = xmax + 0.05*abs(xmax)
    if ymin == ymax:
        if ymin == 0.0:
            ymin = -0.1
            ymax = 0.1
        else:
            ymin = ymin - 0.05*abs(ymin)
            ymax = ymax + 0.05*abs(ymax)

    dx = xmax-xmin
    dy = ymax-ymin
    if autoxmin:
        xmin = xmin - 0.05*dx
    if autoxmax:
        xmax = xmax + 0.05*dx
    if autoymin:
        ymin = ymin - 0.05*dy
    if autoymax:
        ymax = ymax + 0.05*dy

    # Construct some edge specifiers for the histogram2d function call
    xe = np.linspace(xmin, xmax, nx)
    ye = np.linspace(ymin, ymax, ny)
    # Call the numpy histogram2d function
    z0, yedges1, xedges1 = np.histogram2d(datay, datax, bins=(ye, xe))
    # In the contour plots, x and y are the centers of the cells, instead of the edges.
    x = np.zeros([nx-1])
    y = np.zeros([ny-1])
    for i in range(nx-1):
        x[i] = 0.5*(xe[i]+xe[i+1])
    for j in range(ny-1):
        y[j] = 0.5*(ye[j]+ye[j+1])

    # Use numpy histogram2d function to make image
    z_scal = z_imag = z_cont = z_outl = False
    empty = True
    if scalar:
        try:
            z1, yedges1, xedges1 = np.histogram2d(datay, datax, bins=(
                ye, xe), weights=holder.get(scalar.name, only_leafs=only_leafs))
            if summed:
                z_scal = np.ma.masked_where(z0 == 0.0, z1)
            else:
                with np.errstate(divide="ignore", invalid="ignore"):
                    z_scal = np.ma.masked_where(z0 == 0.0, z1/z0)
        except AttributeError:
            holder.new_field(name=default_var, unit="",
                             label="Number of cells", verbose=True)
            scalar = getattr(holder, default_var)
            z_scal = np.ma.masked_where(z0 == 0.0, z0)
        empty = False
    if image:
        try:
            z1, yedges1, xedges1 = np.histogram2d(datay, datax, bins=(
                ye, xe), weights=holder.get(image.name, only_leafs=only_leafs))
            if summed:
                z_imag = np.ma.masked_where(z0 == 0.0, z1)
            else:
                with np.errstate(divide="ignore", invalid="ignore"):
                    z_imag = np.ma.masked_where(z0 == 0.0, z1/z0)
        except AttributeError:
            holder.new_field(name=default_var, unit="",
                             label="Number of cells", verbose=True)
            image = getattr(holder, default_var)
            z_imag = np.ma.masked_where(z0 == 0.0, z0)
        empty = False
    if contour:
        try:
            z1, yedges1, xedges1 = np.histogram2d(datay, datax, bins=(
                ye, xe), weights=holder.get(contour.name, only_leafs=only_leafs))
            if summed:
                z_cont = np.ma.masked_where(z0 == 0.0, z1)
            else:
                with np.errstate(divide="ignore", invalid="ignore"):
                    z_cont = np.ma.masked_where(z0 == 0.0, z1/z0)
        except AttributeError:
            holder.new_field(name=default_var, unit="",
                             label="Number of cells", verbose=True)
            contour = getattr(holder, default_var)
            z_cont = np.ma.masked_where(z0 == 0.0, z0)
        empty = False
    if scatter:
        empty = False

    # If no variable is requested for z/color dimension, store number of cells by default
    if empty:
        holder.new_field(name=default_var, unit="",
                         label="Number of cells", verbose=True)
        scalar = getattr(holder, default_var)
        z_scal = np.ma.masked_where(z0 == 0.0, z0)

    if outline:
        z_outl = z0

    if plot:
        render_map(scalar=scalar, image=image, contour=contour, scatter=scatter, x=x, y=y, z_scal=z_scal,
                   z_imag=z_imag, z_cont=z_cont, z_outl=z_outl, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fname=fname,
                   axes=axes, title=title, new_window=new_window, clear=clear, block=block,
                   resolution=resolution, scalar_args=scalar_args, image_args=image_args, holder=holder,
                   contour_args=contour_args, scatter_args=scatter_args, equal_axes=equal_axes, x_raw=datax, y_raw=datay,
                   outline=outline, outline_args=outline_args, sinks=False, only_leafs=only_leafs,
                   dir_vecs=[["", [0, 0, 0]], [var_x.name, [0, 0, 0]], [var_y.name, [0, 0, 0]]])

    if hasattr(holder, default_var):
        holder.delete_field(default_var)

    if copy:
        return x, y, z_scal, z_imag, z_cont, z_outl
    else:
        return

# =======================================================================================
# Plot a 2D slice through the data cube. The arguments are:
# - scalar     : the scalar field to be plotted, e.g. mydata.density
# - image      : the scalar field to be plotted with an image
# - contour    : the scalar field to be plotted with contours
# - vec        : the vector field to be overplotted, e.g. mydata.velocity
# - stream     : the field for streamlines to be overplotted, e.g. mydata.B
# - direction  : the direction normal to the plane of the slice
# - dx         : the x extent of the slice, in units of scale (see data loader)
# - dy         : the y extent of the slice, in units of scale. If not specified, dy = dx
# - dz         : the thickness of the slice
# - axes       : if specified, the data is plotted on the specified axes (see demo).
# - resolution : number of pixels in the slice.
# - fname      : if specified, the figure is saved to file.
# =======================================================================================


def plot_slice(scalar=False, image=False, contour=False, vec=False, stream=False, axes=None,
               direction="z", dx=0.0, dy=0.0, fname=None, title=None, sinks=True, copy=False,
               origin=[0, 0, 0], resolution=128, new_window=False, update=None,
               clear=True, plot=True, block=False, interpolation="linear", sink_args={},
               scalar_args={}, image_args={}, contour_args={}, vec_args={}, stream_args={},
               outline=False, outline_args={}, lmax=0, slice_direction=None):

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

# =======================================================================================
# Plot a column density through the data cube. The arguments are:
# - scalar     : the scalar field to be plotted, e.g. mydata.density
# - image      : the scalar field to be plotted with an image
# - contour    : the scalar field to be plotted with contours
# - dx         : the x extent of the slice, in units of scale (see data loader)
# - dy         : the y extent of the slice, in units of scale. If not specified, dy = dx
# - dz         : the thickness of the slice
# - axes       : if specified, the data is plotted on the specified axes (see demo).
# - resolution : number of pixels in the slice.
# - fname      : if specified, the figure is saved to file.
# =======================================================================================


def plot_column_density(scalar=False, image=False, contour=False, vec=False, stream=False,
                        direction="z", dx=0.0, dy=0.0, dz=0.0, fname=None, axes=None, title=None,
                        origin=[0, 0, 0], resolution=128, sinks=True, summed=True, copy=False,
                        new_window=False, update=None, clear=True, plot=True, block=False, nz=0,
                        interpolation="linear", verbose=False, outline=False, outline_args={},
                        scalar_args={}, image_args={}, contour_args={}, vec_args={}, stream_args={},
                        sink_args={}, lmax=0):

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

# =======================================================================================
# Compute a vector perpendicular to the input vector
# =======================================================================================


def perpendicular_vector(v):

    # x = y = z = 0 is not an acceptable solution
    if v[0] == v[1] == v[2] == 0:
        raise ValueError("zero-vector")

    if v[2] == 0:
        return [-v[1], v[0], 0]
    else:
        return [1.0, 1.0, -1.0 * (v[0] + v[1]) / v[2]]

# =======================================================================================
# Separate arguments to plotting functions between osyris specific arguments and
# Matplotlib arguments.
# =======================================================================================


def parse_arguments(args, args_osyris, args_plot):

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
    if args_osyris["vmin"] == args_osyris["vmax"]:
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

# =======================================================================================
# Find direction vectors for slice
# =======================================================================================


def get_slice_direction(holder, direction, dx=0, dy=0, origin=[0, 0, 0]):

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

# =======================================================================================
# Use matplotlib to plot histogram, slice or column density maps
# =======================================================================================


def render_map(scalar=False, image=False, contour=False, scatter=False, vec=False, stream=False, outline=False, x=0, y=0,
               z_scal=0, z_imag=0, z_cont=0, z_outl=0, u_vect=0, v_vect=0, w_vect=0, u_strm=0, v_strm=0,
               w_strm=0, fname=None, axes=None, title=None, sinks=True, new_window=False,
               clear=True, block=False, xmin=0, xmax=0, ymin=0, ymax=0,
               resolution=128, scalar_args={}, image_args={}, contour_args={}, vec_args={},
               stream_args={}, scatter_args={}, outline_args={}, sink_args={}, dz=0, holder=None,
               thePlane=0, origin=[0, 0, 0], dir_vecs=[["z", [0, 0, 0]], ["x", [0, 0, 0]], ["y", [0, 0, 0]]], x_raw=None, y_raw=None, equal_axes=True,
               only_leafs=True):

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
        outline_args_plot = {"levels": [np.nanmin(z_outl)], "colors": "grey"}
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
