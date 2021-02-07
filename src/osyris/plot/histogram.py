# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2019 Osyris contributors (https://github.com/nvaytet/osyris)
# @author Neil Vaytet

import numpy as np
from ..core import Plot
from .. import units
from .render import render
from ..core.tools import to_bin_centers, finmin, finmax
from .parser import parse_layer, get_norm
# from .engine import OsyrisField
from scipy.stats import binned_statistic_2d


def histogram(x, y, *layers,
    mode=None,
    # image=None,
    # contour=None,
    # scatter=None,
                   ax=None,
                   logx=False,
                   logy=False,
                   loglog=False,
                   norm=None,
                   # block=False,
                   cbar=True,
                   # clear=True,
                   # contour=False,
                   # contour_args={},
                   # copy=False,
                   # equal_axes=False,
                   filename=None,
                   # image=False,
                   # image_args={},
                   # new_window=False,
                   # outline=False,
                   # outline_args={},
                   # plot=True,
                   resolution=256,
                   # scalar_args={},
                   # scatter=False,
                   # scatter_args={},
                   operation="sum",
                   title=None,
                   # update=None,
                   xmin=None,
                   xmax=None,
                   ymin=None,
                   ymax=None,
                   vmin=None,
                   vmax=None,
                   **kwargs):
    """
    Plot a 2D histogram with two variables as input.
    """
    if loglog:
        logx = logy = True

    nx = resolution
    ny = resolution

    # # Get the data values and units
    # datax = holder.get(var_x.name)
    # datay = holder.get(var_y.name)
    xlabel = x.name # var_x.label+" ["+var_x.unit+"]"
    ylabel = y.name # var_y.label+" ["+var_y.unit+"]"
    default_var = "histo_cell_density"

    # Define plotting range
    autoxmin = False
    autoxmax = False
    autoymin = False
    autoymax = False

    if xmin is None:
        xmin = finmin(x.values)
        autoxmin = True
    if xmax is None:
        xmax = finmax(x.values)
        autoxmax = True
    if ymin is None:
        ymin = finmin(y.values)
        autoymin = True
    if ymax is None:
        ymax = finmax(y.values)
        autoymax = True

    if logx:
        [xmin, xmax] = np.log10([xmin, xmax])
    if logy:
        [ymin, ymax] = np.log10([ymin, ymax])
       #  xmax = np.log10(xmax)
    # ymin = np.log10(ymin)
    # ymax = np.log10(ymax)

    # try:
    #     xmin += 0
    # except TypeError:
    #     datax[np.isneginf(datax)] = np.inf
    #     xmin = np.nanmin(datax)
    #     autoxmin = True
    # try:
    #     xmax += 0
    # except TypeError:
    #     datax[np.isinf(datax)] = np.NINF
    #     xmax = np.nanmax(datax)
    #     autoxmax = True
    # try:
    #     ymin += 0
    # except TypeError:
    #     datay[np.isneginf(datay)] = np.inf
    #     ymin = np.nanmin(datay)
    #     autoymin = True
    # try:
    #     ymax += 0
    # except TypeError:
    #     datay[np.isinf(datay)] = np.NINF
    #     ymax = np.nanmax(datay)
    #     autoymax = True

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

    # Construct some bin edges
    # xedges = np.linspace(xmin, xmax, nx+1)
    # yedges = np.linspace(ymin, ymax, ny+1)
    if logx:
        xedges = np.logspace(xmin, xmax, nx+1)
    else:
        xedges = np.linspace(xmin, xmax, nx+1)
    if logy:
        yedges = np.logspace(ymin, ymax, ny+1)
    else:
        yedges = np.linspace(ymin, ymax, ny+1)


    # # Call the numpy histogram2d function
    # counts, _, _ = np.histogram2d(y, x, bins=(yedges, xedges))
    # masked = np.ma.masked_where(z0 == 0.0, z0)

    # In the contour plots, x and y are the centers of the cells, instead of the edges.
    xcenters = to_bin_centers(xedges)
    ycenters = to_bin_centers(yedges)

    # x = np.zeros([nx-1])
    # y = np.zeros([ny-1])
    # for i in range(nx-1):
    #     x[i] = 0.5*(xe[i]+xe[i+1])
    # for j in range(ny-1):
    #     y[j] = 0.5*(ye[j]+ye[j+1])

    # # Use numpy histogram2d function to make image
    # z_scal = z_imag = z_cont = z_outl = False

    to_render = {}
    # "contf": None, "image": None, "contour": None,
    #              "outline": None, "scatter": None}

    to_process = {"counts": np.ones_like(x.values)}

    empty = True

    # cell_count = OsyrisField(name=default_var,
    #                    unit="", label="Number of cells")

    operations = {}
    if layers is not None:
        for layer in layers:
            data, settings, params = parse_layer(layer, mode=mode, norm=norm,
                vmin=vmin, vmax=vmax, operation=operation, **kwargs)
            to_process[data.name] = data.values
            to_render[data.name] = {"mode": settings["mode"],
                                    "params": params,
                                    "unit": data.unit.units}
            operations[data.name] = settings["operation"]

    # print("========")
    # print(to_process)
    # print(to_render)
    # print("========")

    if (operation == "mean") and "sum" in operations.values():
        counts, _, _ = np.histogram2d(y.values, x.values, bins=(yedges, xedges))




    # if image is not None:
    #     if image is True:
    #         image = cell_count
    #         to_render["image"] = z1
    #     else:
    #         to_process["image"] = image.values
    #     empty = False
    # if contour:
    #     if contour is True:
    #         contour = cell_count
    #         to_render["contour"] = z1
    #     else:
    #         to_process["contour"] = contour.values
    #     empty = False
    # if scatter:
    #     if scatter is True:
    #         scatter = cell_count
    #     empty = False
    # if scalar or empty:
    #     if hasattr(scalar, "values"):
    #         to_process["scalar"] = scalar.values
    #         empty = False
    #     else:
    #         scalar = cell_count
    #         to_render["scalar"] = z1

    # # if outline:
    # #     to_render["outline"] = z0

    # if len(to_process) > 0:

    # print(xedges)
    # print(yedges)
    binned, _, _, _ = binned_statistic_2d(x=y.values, y=x.values, values=list(to_process.values()), statistic=operation, bins=[yedges, xedges])

    # Here we assume that dictionary retains order of insertion: counts
    # are the first key
    mask = binned[0] == 0.0
    for i, key in enumerate(set(to_process.keys()) - set(["counts"])):
        data = binned[i + 1]
        print(key, operations[key], operation)
        if operations[key] != operation:
            if operation == "sum":
                with np.errstate(invalid="ignore"):
                    data /= binned[0]
            else:
                data *= counts
        to_render[key]["data"] = np.ma.masked_where(mask, data)

    if len(to_render) == 0:
        to_render["counts"] = {"data": np.ma.masked_where(mask, binned[0]),
                               "mode": mode,
                               "params":{
                               "norm": get_norm(norm=norm,
                                                   vmin=vmin,
                                                   vmax=vmax),
                               "vmin": vmin,
                               "vmax": vmax},
                               "unit": 1.0 * units.dimensionless}


    figure = render(x=xcenters, y=ycenters,
        data=to_render, logx=logx, logy=logy)
        # norm=norm, vmin=vmin, vmax=vmax,
        # **kwargs)

    # figure = render(scalar=scalar, image=image, contour=contour, vec=vec, stream=stream, x=x, y=y, z_scal=to_render["scalar"],
    #                z_imag=to_render["image"], z_cont=to_render["contour"], u_vect=to_render["u_vect"], v_vect=to_render["v_vect"], w_vect=to_render["w_vect"], u_strm=to_render["u_strm"],
    #                v_strm=to_render["v_strm"], w_strm=to_render["w_strm"], xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fname=fname,
    #                axes=axes, title=title, sinks=sinks, new_window=new_window, clear=clear, block=block,
    #                resolution=resolution, thePlane=[
    #                    a_plane, b_plane, c_plane, d_plane],
    #                origin=origin, dir_vecs=dir_vecs, scalar_args=scalar_args, image_args=image_args,
    #                contour_args=contour_args, vec_args=vec_args, stream_args=stream_args, outline=outline,
    #                outline_args=outline_args, sink_args=sink_args, x_raw=datax, y_raw=datay, holder=holder)

    # figure["ax"].set_xscale("log")
    # figure["ax"].set_yscale("log")

    figure["ax"].set_xlabel(x.label)
    figure["ax"].set_ylabel(y.label)


    return Plot(x=xcenters, y=ycenters, layers=to_render, fig=figure["fig"], ax=figure["ax"])



    # if plot:
    #     render_map(scalar=scalar, image=image, contour=contour, scatter=scatter, x=x, y=y, z_scal=to_render["scalar"],
    #                z_imag=to_render["image"], z_cont=to_render["contour"], z_outl=to_render["outline"], xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fname=fname,
    #                axes=axes, title=title, new_window=new_window, clear=clear, block=block,
    #                resolution=resolution, scalar_args=scalar_args, image_args=image_args, holder=holder,
    #                contour_args=contour_args, scatter_args=scatter_args, equal_axes=equal_axes, x_raw=datax, y_raw=datay,
    #                outline=outline, outline_args=outline_args, sinks=False,
    #                dir_vecs=[["", [0, 0, 0]], [var_x.name, [0, 0, 0]], [var_y.name, [0, 0, 0]]])

    # if hasattr(holder, default_var):
    #     holder.delete_field(default_var)

    # if copy:
    #     return x, y, to_render
    # else:
    #     return
