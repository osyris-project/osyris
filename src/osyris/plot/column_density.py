# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np
import numpy.ma as ma
from pint.quantity import Quantity
from .slice import get_slice_direction
from .render import render
from .scatter import scatter
from .parser import parse_layer
from ..core import Plot, Array
from ..core.tools import apply_mask
from .utils import evaluate_on_grid_3d


def _add_scatter(to_scatter, origin, dir_vecs, dx, dy, ax):
    xyz = to_scatter[0]["data"] - origin
    viewport = max(dx.magnitude, dy.magnitude)
    radius = None
    if "s" in to_scatter[0]["params"]:
        size = to_scatter[0]["params"]["s"]
        if isinstance(size, Array) or isinstance(size, Quantity):
            radius = size.to(dx.units)
            to_scatter[0]["params"]["s"] = radius
    if radius is None:
        # Fudge factor to select sinks close to the plane
        radius = Array(values=viewport * 0.05, unit=dx.units)
    dist1 = np.sum(xyz * dir_vecs[0], axis=1)
    global_selection = np.arange(len(to_scatter[0]["data"]))
    select = np.ravel(np.where(np.abs(dist1) <= radius))
    global_selection = global_selection[select]
    if len(select) > 0:
        # Project coordinates onto the plane by taking dot product with axes vectors
        coords = xyz[select]
        datax = np.inner(coords, dir_vecs[1])
        datay = np.inner(coords, dir_vecs[2])
        if dx is not None:
            # Limit selection further by using distance from center
            dist2 = coords
            select2 = np.ravel(
                np.where(np.abs(dist2.norm.values) <= viewport * 0.6 * np.sqrt(2.0)))
            datax = datax[select2]
            datay = datay[select2]
            global_selection = global_selection[select2]
        if "c" in to_scatter[0]["params"]:
            # TODO: also check that parents are the same to ensure size match?
            if isinstance(to_scatter[0]["params"]["c"], Array):
                to_scatter[0]["params"]["c"] = to_scatter[0]["params"]["c"][
                    global_selection]
        scatter(x=datax, y=datay, ax=ax, **to_scatter[0]["params"])


def column_density(*layers,
                   direction="z",
                   dx=None,
                   dy=None,
                   dz=None,
                   fname=None,
                   title=None,
                   plot=True,
                   mode=None,
                   norm=None,
                   vmin=None,
                   vmax=None,
                   origin=None,
                   resolution=256,
                   ax=None,
                   **kwargs):
    """
    Plot a 2D slice through the data domain.
    """

    if isinstance(layers, Array):
        layers = [layers]

    to_process = []
    to_render = []
    operations = []
    to_scatter = []
    for layer in layers:
        data, settings, params = parse_layer(layer=layer,
                                             mode=mode,
                                             norm=norm,
                                             vmin=vmin,
                                             vmax=vmax,
                                             **kwargs)
        if settings["mode"] == "scatter":
            to_scatter.append({"data": data, "params": params})
        else:
            to_process.append(data)
            to_render.append({
                "mode": settings["mode"],
                "params": params,
                "unit": data.unit.units,
                "name": data.name
            })
            operations.append(settings["operation"])

    dataset = to_process[0].parent.parent

    # Set window size
    if dy is None:
        dy = dx
    if dz is None:
        dz = dx
    if dx is not None and not isinstance(dx, Quantity):
        dx *= dataset["amr"]["xyz"].unit
    if dy is not None and not isinstance(dy, Quantity):
        dy *= dataset["amr"]["xyz"].unit
    if dz is not None and not isinstance(dz, Quantity):
        dz *= dataset["amr"]["xyz"].unit

    dir_vecs, origin = get_slice_direction(direction=direction,
                                           dataset=dataset,
                                           dx=dx,
                                           dy=dy,
                                           origin=origin)

    # Distance to the plane
    xyz = dataset["amr"]["xyz"] - origin
    diagonal = dataset["amr"]["dx"] * 0.5 * np.sqrt(dataset.meta["ndim"])
    dist_to_plane = np.sum(xyz * dir_vecs[0], axis=1)
    # Create an array of indices to allow further narrowing of the selection below
    global_indices = np.arange(len(dataset["amr"]["dx"]))
    # Select cells close to the plane, including factor of sqrt(ndim)
    close_to_plane = np.ravel(
        np.where(np.abs(dist_to_plane) <= 0.5 * dz * np.sqrt(dataset.meta["ndim"])))
    indices_close_to_plane = global_indices[close_to_plane]

    if len(indices_close_to_plane) == 0:
        raise RuntimeError("No cells were selected to construct the column density. "
                           "The resulting figure would be empty.")

    xmin = None
    if dx is not None:
        xmin = -0.5 * dx.magnitude
        xmax = xmin + dx.magnitude
        ymin = -0.5 * dy.magnitude
        ymax = ymin + dy.magnitude
        zmin = -0.5 * dz.magnitude
        zmax = zmin + dz.magnitude
        # Limit selection further by using distance from center
        radial_distance = xyz[indices_close_to_plane] - 0.5 * dataset["amr"]["dx"][
            indices_close_to_plane] * np.sqrt(dataset.meta["ndim"])
        radial_selection = np.ravel(
            np.where(
                np.abs(radial_distance.norm.values) <=
                max(dx.magnitude, dy.magnitude, dz.magnitude) * 0.6 * np.sqrt(2.0)))
        indices_close_to_plane = indices_close_to_plane[radial_selection]

    # Project coordinates onto the plane by taking dot product with axes vectors
    coords = xyz[indices_close_to_plane]
    datax = np.inner(coords, dir_vecs[1])
    datay = np.inner(coords, dir_vecs[2])
    dataz = np.inner(coords, dir_vecs[0])
    datadx = dataset["amr"]["dx"][indices_close_to_plane] * 0.5

    if xmin is None:
        xmin = (datax - datadx).min().values
        xmax = (datax + datadx).max().values
        ymin = (datay - datadx).min().values
        ymax = (datay + datadx).max().values
        zmin = (dataz - datadx).min().values
        zmax = (dataz + datadx).max().values

    scalar_layer = []
    to_binning = []  # contains the variables in cells close to the plane
    for ind in range(len(to_process)):
        if to_render[ind]["mode"] in ["vec", "stream", "lic"]:
            if to_process[ind].ndim < 3:
                uv = to_process[ind].array[indices_close_to_plane]
            else:
                uv = np.inner(
                    to_process[ind].array.take(indices_close_to_plane, axis=0),
                    dir_vecs[1:])
            w = None
            if "color" in to_render[ind]["params"]:
                if isinstance(to_render[ind]["params"]["color"], Array):
                    w = to_render[ind]["params"]["color"].norm.values
                elif isinstance(to_render[ind]["params"]["color"], np.ndarray):
                    w = to_render[ind]["params"]["color"]
            if w is None:
                w = np.linalg.norm(uv, axis=1)
            else:
                w = w.take(indices_close_to_plane, axis=0)
            to_binning.append(apply_mask(uv[:, 0]))
            to_binning.append(apply_mask(uv[:, 1]))
            to_binning.append(w)
            scalar_layer.append(False)
        else:
            to_binning.append(
                apply_mask(to_process[ind].norm.values[indices_close_to_plane]))
            scalar_layer.append(True)

    # Create a grid of pixel centers
    if isinstance(resolution, int):
        resolution = {'x': resolution, 'y': resolution, 'z': resolution}
    xspacing = (xmax - xmin) / (resolution['x'] + 1)
    yspacing = (ymax - ymin) / (resolution['y'] + 1)
    zspacing = (zmax - zmin) / (resolution['z'] + 1)

    xcenters = np.linspace(xmin + 0.5 * xspacing, xmax - 0.5 * xspacing,
                           resolution['x'])
    ycenters = np.linspace(ymin + 0.5 * yspacing, ymax - 0.5 * yspacing,
                           resolution['y'])
    zcenters = np.linspace(zmin + 0.5 * zspacing, zmax - 0.5 * zspacing,
                           resolution['z'])
    xgrid, ygrid, zgrid = np.meshgrid(xcenters, ycenters, zcenters, indexing='xy')
    pixel_positions = xgrid.reshape(xgrid.shape + (1, )) * dir_vecs[1] + ygrid.reshape(
        ygrid.shape + (1, )) * dir_vecs[2] + zgrid.reshape(zgrid.shape +
                                                           (1, )) * dir_vecs[0]

    # Evaluate the values of the data layers at the grid positions
    binned = evaluate_on_grid_3d(
        cell_positions_in_new_basis=np.array(
            [apply_mask(datax.array),
             apply_mask(datay.array),
             apply_mask(dataz.array)]).T,
        cell_positions_in_original_basis=coords.array,
        cell_values=np.array(to_binning),
        cell_sizes=datadx.array,
        grid_lower_edge_in_new_basis=np.array([xmin, ymin, zmin]),
        grid_spacing_in_new_basis=np.array([xspacing, yspacing, zspacing]),
        grid_positions_in_original_basis=pixel_positions,
        ndim=dataset.meta["ndim"]).sum(axis=0)

    # binned = binned.sum(axis=0)

    # Mask NaN values
    mask = np.isnan(binned[-1, ...])
    mask_vec = np.broadcast_to(mask.reshape(*mask.shape, 1), mask.shape + (3, ))

    # Now we fill the arrays to be sent to the renderer, also constructing vectors
    counter = 0
    for ind in range(len(to_render)):
        if scalar_layer[ind]:
            to_render[ind]["data"] = ma.masked_where(mask,
                                                     binned[counter, ...],
                                                     copy=False)
            counter += 1
        else:
            to_render[ind]["data"] = ma.masked_where(mask_vec,
                                                     np.array([
                                                         binned[counter, ...].T,
                                                         binned[counter + 1, ...].T,
                                                         binned[counter + 2, ...].T
                                                     ]).T,
                                                     copy=False)
            counter += 3

    to_return = {"x": xcenters, "y": ycenters, "layers": to_render}
    if plot:
        # Render the map
        figure = render(x=xcenters, y=ycenters, data=to_render, ax=ax)
        figure["ax"].set_xlabel(dataset["amr"]["xyz"].x.label)
        figure["ax"].set_ylabel(dataset["amr"]["xyz"].y.label)
        if ax is None:
            figure["ax"].set_aspect("equal")

        # Add scatter layer
        if len(to_scatter) > 0:
            _add_scatter(to_scatter=to_scatter,
                         origin=origin,
                         dir_vecs=dir_vecs,
                         dx=dx,
                         dy=dy,
                         ax=figure["ax"])

        figure["ax"].set_xlim(xmin, xmax)
        figure["ax"].set_ylim(ymin, ymax)

        to_return.update({"fig": figure["fig"], "ax": figure["ax"]})

    return Plot(**to_return)
