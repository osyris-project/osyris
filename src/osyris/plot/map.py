# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)

import numpy as np
import numpy.ma as ma
from pint import Quantity
from typing import Union
from .direction import get_direction
from .render import render
from .scatter import scatter
from .parser import parse_layer
from ..core import Plot, Array, Vector
from ..core.tools import apply_mask
from .utils import evaluate_on_grid


def _add_scatter(to_scatter, origin, dir_vecs, dx, dy, ax, map_unit):
    xyz = to_scatter[0]["data"] - origin
    viewport = np.maximum(dx, dy)
    radius = None
    if "s" in to_scatter[0]["params"]:
        size = to_scatter[0]["params"]["s"]
        if isinstance(size, Array) or isinstance(size, Quantity):
            radius = size.to(dx.units)
            to_scatter[0]["params"]["s"] = radius
    if radius is None:
        # Fudge factor to select sinks close to the plane
        radius = Array(values=viewport * 0.05)
    dist_to_plane = xyz.dot(dir_vecs["normal"])
    global_selection = np.arange(len(to_scatter[0]["data"]))
    select = (np.abs(dist_to_plane) <= radius).values
    global_selection = global_selection[select]
    if len(select) > 0:
        # Project coordinates onto the plane by taking dot product with axes vectors
        coords = xyz[select]
        datax = coords.dot(dir_vecs["pos_u"])
        datay = coords.dot(dir_vecs["pos_v"])

        if dx is not None:
            # Limit selection further by using distance from center
            select2 = (np.abs(coords.norm) <= viewport * 0.6 * np.sqrt(2.0)).values
            datax = datax[select2]
            datay = datay[select2]
            global_selection = global_selection[select2]
        if "c" in to_scatter[0]["params"]:
            # TODO: also check that parents are the same to ensure size match?
            if isinstance(to_scatter[0]["params"]["c"], Array):
                to_scatter[0]["params"]["c"] = to_scatter[0]["params"]["c"][
                    global_selection
                ]
        datax.name = dir_vecs["pos_u"].name
        datay.name = dir_vecs["pos_v"].name
        scatter(
            x=datax.to(map_unit), y=datay.to(map_unit), ax=ax, **to_scatter[0]["params"]
        )


def map(
    *layers,
    direction: Union[str, list] = "z",
    dx: Quantity = None,
    dy: Quantity = None,
    dz: Quantity = None,
    filename: str = None,
    title: str = None,
    plot: bool = True,
    mode: str = None,
    norm: str = None,
    vmin: float = None,
    vmax: float = None,
    origin: Array = None,
    resolution: Union[int, dict] = None,
    operation: str = "sum",
    ax: object = None,
    **kwargs,
) -> Plot:
    """
    Create a 2D spatial map of a region inside a simulation domain.
    By default, the map represents a plane with zero thickness.
    A thick slab or cube can also be computed by specifying a thickness via the
    ``dz`` argument. In this case, the resulting 3D box is integrated along the ``z``
    direction before being sent to the image rendering.

    :param layers: Dicts or Arrays representing the quantities to be mapped onto the
        generated image.

    :param direction: The vector normal to the map. Possible choices are:

       * ``'x'``, ``'y'``, or ``'z'`` representing the cartesian axes
       * a list of 3 numbers representing the components of the vector,
         e.g. ``[1, 0.5, 2]``
       * ``'top'`` or ``'side'`` for automatic top or side view of a disk, according to
         the angular momentum computed around the center of the plotted region

    :param dx: The horizontal size of the plotted region. Default is ``None``,
        in which case the entire horizontal range of the simulation domain is plotted.

    :param dy: The vertical size of the plotted region. If not specified, it will
        either be equal to ``dx`` if ``dx`` is not ``None``, or the entire vertical
        range of the simulation domain if ``dx`` is ``None``. Default is ``None``.

    :param dz: The depth range over which the ``z`` dimension is to be integrated.
        Default is ``None``, in which case a plane with no thickness is plotted.

    :param filename: If specified, the returned figure is also saved to file.
        Default is ``None``.

    :param title: The title of the figure. Default is ``None``.

    :param plot: Make a plot if ``True``. If not, just return the ``Plot`` object
        containing the data that would be used to generate the plot.
        Default is ``True``.

    :param mode: The rendering mode for the map. Possible choices are ``'image'``,
        ``'contourf'``, and ``'contour'`` for scalar Arrays, ``'vec'`` and
        ``'stream'`` for vector quantities. Default is ``None``, which selects the
        ``render_mode`` set in the user configuration file (``'image'`` by default).

    :param norm: The colormap normalization. Possible values are ``'linear'`` and
        ``'log'``. Default is ``None`` (= ``'linear'``).

    :param vmin: Minimum value for colorbar range. Default is ``None``.

    :param vmax: Maximum value for colorbar range. Default is ``None``.

    :param origin: An Array describing the position of the center of the map
        (with 2 or 3 components depending on the dimensionality of the simulation).

    :param resolution: Resolution of the generated map. This can either be an
        integer or a dict. In the case of an integer, it represents the number of
        pixels used for the horizontal and vertical dimensions. For a dictionary,
        the following syntax should be used: ``resolution={'x': 128, 'y': 192}``.
        Default is ``256``.

    :param operation: Numpy operation to apply along the ``z`` dimension if ``dz`` is
        not ``None``. Example values are ``'sum'``, ``'nansum'``, ``'nanmin'``,
        ``'nanmax'``, ``'nanmean'``, ``'mean'``, ``'min'``, and ``'max'``.
        Default is ``'sum'``.

    :param ax: A matplotlib axes inside which the figure will be plotted.
        Default is ``None``, in which case some new axes a created.
    """

    if isinstance(layers, Array):
        layers = [layers]

    to_process = []
    to_render = []
    to_scatter = []
    for layer in layers:
        data, settings, params = parse_layer(
            layer=layer, mode=mode, norm=norm, vmin=vmin, vmax=vmax, **kwargs
        )
        if settings["mode"] == "scatter":
            to_scatter.append({"data": data, "params": params})
        else:
            to_process.append(data)
            to_render.append(
                {
                    "mode": settings["mode"],
                    "params": params,
                    "unit": data.unit,
                    "name": data.name,
                }
            )

    dataset = to_process[0].parent.parent
    ndim = dataset.meta["ndim"]

    thick = dz is not None

    spatial_unit = dataset["amr"]["position"].unit
    map_unit = spatial_unit

    # Set window size
    if dx is not None:
        map_unit = dx.units
        dx = dx.to(spatial_unit)
    dy = dx if dy is None else dy.to(spatial_unit)
    dz = dx if dz is None else dz.to(spatial_unit)

    if origin is None:
        origin = Vector(*[0 for n in range(ndim)], unit=spatial_unit)

    dir_vecs = get_direction(
        direction=direction, dataset=dataset, dx=dx, dy=dy, origin=origin
    )

    # Distance to the plane
    diagonal = np.sqrt(ndim)
    xyz = dataset["amr"]["position"] - origin
    selection_distance = 0.5 * diagonal * (dz if thick else dataset["amr"]["dx"])

    normal = dir_vecs["normal"]
    vec_u = dir_vecs["pos_u"]
    vec_v = dir_vecs["pos_v"]

    dist_to_plane = xyz.dot(normal)
    # Create an array of indices to allow further narrowing of the selection below
    global_indices = np.arange(len(dataset["amr"]["dx"]))
    # Select cells close to the plane, including factor of sqrt(ndim)
    close_to_plane = (np.abs(dist_to_plane) <= selection_distance).values
    indices_close_to_plane = global_indices[close_to_plane]

    if len(indices_close_to_plane) == 0:
        raise RuntimeError(
            "No cells were selected to construct the map. "
            "The resulting figure would be empty."
        )

    xmin = None
    if dx is not None:
        xmin = -0.5 * dx.magnitude
        xmax = xmin + dx.magnitude
        ymin = -0.5 * dy.magnitude
        ymax = ymin + dy.magnitude
        zmin = -0.5 * dz.magnitude
        zmax = zmin + dz.magnitude
        # Limit selection further by using distance from center
        radial_distance = (
            xyz[indices_close_to_plane]
            - 0.5 * dataset["amr"]["dx"][indices_close_to_plane] * diagonal
        )
        radial_selection = (
            np.abs(radial_distance.norm.values)
            <= max(dx.magnitude, dy.magnitude, dz.magnitude) * 0.6 * diagonal
        )
        indices_close_to_plane = indices_close_to_plane[radial_selection]

    # Project coordinates onto the plane by taking dot product with axes vectors
    coords = xyz[indices_close_to_plane]
    datax = coords.dot(vec_u)
    datay = coords.dot(vec_v)
    dataz = coords.dot(normal)
    datadx = dataset["amr"]["dx"][indices_close_to_plane] * 0.5

    if xmin is None:
        xmin = (datax - datadx).min().values
        xmax = (datax + datadx).max().values
        ymin = (datay - datadx).min().values
        ymax = (datay + datadx).max().values
        zmin = (dataz - datadx).min().values
        zmax = (dataz + datadx).max().values
        dx = (xmax - xmin) * datadx.unit
        dy = (ymax - ymin) * datadx.unit

    scalar_layer = []
    to_binning = []  # contains the variables in cells close to the plane
    for ind in range(len(to_process)):
        if to_render[ind]["mode"] in ["vec", "stream", "lic"]:
            uv = to_process[ind][indices_close_to_plane]
            if to_process[ind].z is None:
                u = uv.x.values
                v = uv.y.values
            else:
                u = uv.dot(vec_u).values
                v = uv.dot(vec_v).values

            w = None
            if isinstance(to_render[ind]["params"].get("color"), (Array, Vector)):
                w = to_render[ind]["params"]["color"].norm.values[
                    indices_close_to_plane
                ]
            else:
                w = u * u
                w += v * v
                w = np.sqrt(w)
            to_binning.append(apply_mask(u))
            to_binning.append(apply_mask(v))
            to_binning.append(w)
            scalar_layer.append(False)
        else:
            to_binning.append(
                apply_mask(to_process[ind].norm.values[indices_close_to_plane])
            )
            scalar_layer.append(True)

    # Create a grid of pixel centers
    default_resolution = 256
    if resolution is None:
        resolution = default_resolution
    if isinstance(resolution, int):
        resolution = {"x": resolution, "y": resolution}
    else:
        for xy in "xy":
            if xy not in resolution:
                resolution[xy] = default_resolution
    xspacing = (xmax - xmin) / resolution["x"]
    yspacing = (ymax - ymin) / resolution["y"]
    xcenters = np.linspace(
        xmin + 0.5 * xspacing, xmax - 0.5 * xspacing, resolution["x"]
    )
    ycenters = np.linspace(
        ymin + 0.5 * yspacing, ymax - 0.5 * yspacing, resolution["y"]
    )

    if thick:
        if "z" not in resolution:
            resolution["z"] = round((zmax - zmin) / (0.5 * (xspacing + yspacing)))
        zspacing = (zmax - zmin) / resolution["z"]
        zcenters = np.linspace(
            zmin + 0.5 * zspacing, zmax - 0.5 * zspacing, resolution["z"]
        )
    else:
        zmin = 0.0
        zspacing = zmax - zmin
        zcenters = [0.0]

    xg, yg, zg = np.meshgrid(xcenters, ycenters, zcenters, indexing="ij")
    xgrid = xg.T
    ygrid = yg.T
    zgrid = zg.T
    u_array = np.array(
        [
            vec_u.x.values,
            vec_u.y.values,
            vec_u.z.values if vec_u.z is not None else np.zeros_like(vec_u.x.values),
        ]
    )
    v_array = np.array(
        [
            vec_v.x.values,
            vec_v.y.values,
            vec_v.z.values if vec_v.z is not None else np.zeros_like(vec_v.x.values),
        ]
    )
    n_array = np.array(
        [
            normal.x.values,
            normal.y.values,
            normal.z.values if normal.z is not None else np.zeros_like(normal.x.values),
        ]
    )

    pixel_positions = (
        xgrid.reshape(xgrid.shape + (1,)) * u_array
        + ygrid.reshape(ygrid.shape + (1,)) * v_array
        + zgrid.reshape(zgrid.shape + (1,)) * n_array
    )

    # Evaluate the values of the data layers at the grid positions
    div = dx.to(datadx.unit).magnitude
    binned = evaluate_on_grid(
        cell_positions_in_new_basis_x=apply_mask(datax.values / div),
        cell_positions_in_new_basis_y=apply_mask(datay.values / div),
        cell_positions_in_new_basis_z=apply_mask(dataz.values / div),
        cell_positions_in_original_basis_x=coords.x.values / div,
        cell_positions_in_original_basis_y=(
            coords.y.values / div if coords.y is not None else None
        ),
        cell_positions_in_original_basis_z=(
            coords.z.values / div if coords.z is not None else None
        ),
        cell_values=np.array(to_binning),
        cell_sizes=datadx.values / div,
        grid_lower_edge_in_new_basis_x=xmin / div,
        grid_lower_edge_in_new_basis_y=ymin / div,
        grid_lower_edge_in_new_basis_z=zmin / div,
        grid_spacing_in_new_basis_x=xspacing / div,
        grid_spacing_in_new_basis_y=yspacing / div,
        grid_spacing_in_new_basis_z=zspacing / div,
        grid_positions_in_original_basis=pixel_positions / div,
        ndim=ndim,
    )

    # Apply operation along depth
    binned = getattr(np, operation)(binned, axis=1)

    # Handle thick maps
    if thick and ((operation == "sum") or (operation == "nansum")):
        binned *= zspacing
        for layer in to_render:
            layer["unit"] = layer["unit"] * dataz.unit

    # Mask NaN values
    mask = np.isnan(binned[-1, ...])
    mask_vec = np.broadcast_to(mask.reshape(*mask.shape, 1), mask.shape + (3,))

    # Now we fill the arrays to be sent to the renderer, also constructing vectors
    counter = 0
    for ind in range(len(to_render)):
        if scalar_layer[ind]:
            to_render[ind]["data"] = ma.masked_where(
                mask, binned[counter, ...], copy=False
            )
            counter += 1
        else:
            to_render[ind]["data"] = ma.masked_where(
                mask_vec,
                np.array(
                    [
                        binned[counter, ...].T,
                        binned[counter + 1, ...].T,
                        binned[counter + 2, ...].T,
                    ]
                ).T,
                copy=False,
            )
            counter += 3

    scale_ratio = (1.0 * spatial_unit).to(map_unit).magnitude
    xcenters *= scale_ratio
    ycenters *= scale_ratio

    to_return = {
        "x": xcenters,
        "y": ycenters,
        "layers": to_render,
        "filename": filename,
    }
    if plot:
        # Render the map
        figure = render(x=xcenters, y=ycenters, data=to_render, ax=ax)
        figure["ax"].set_xlabel(
            Array(values=0, unit=map_unit, name=dir_vecs["pos_u"].name).label
        )
        figure["ax"].set_ylabel(
            Array(values=0, unit=map_unit, name=dir_vecs["pos_v"].name).label
        )
        if ax is None:
            figure["ax"].set_aspect("equal")

        # Add scatter layer
        if len(to_scatter) > 0:
            _add_scatter(
                to_scatter=to_scatter,
                origin=origin,
                dir_vecs=dir_vecs,
                dx=dx,
                dy=dy,
                ax=figure["ax"],
                map_unit=map_unit,
            )

        xmin *= scale_ratio
        xmax *= scale_ratio
        ymin *= scale_ratio
        ymax *= scale_ratio
        figure["ax"].set_xlim(xmin, xmax)
        figure["ax"].set_ylim(ymin, ymax)

        to_return.update({"fig": figure["fig"], "ax": figure["ax"]})

    return Plot(**to_return)
