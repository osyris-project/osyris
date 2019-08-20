# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2019 Osyris contributors (https://github.com/nvaytet/osyris)
# @author Neil Vaytet

import numpy as np
from .plot import render_map


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
