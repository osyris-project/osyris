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
from ..core.tools import to_bin_centers, apply_mask
from scipy.stats import binned_statistic_2d


def _add_scatter(datax, datay, to_scatter, origin, datadx, dir_vecs, dx, dy, ax):
    xyz = to_scatter[0]["data"] - origin
    viewport = max(dx.magnitude, dy.magnitude)
    radius = None
    if "s" in to_scatter[0]["params"]:
        size = to_scatter[0]["params"]["s"]
        if isinstance(size, Array) or isinstance(size, Quantity):
            radius = size.to(datax.unit.units)
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


def plane(*layers,
          direction="z",
          dx=None,
          dy=None,
          fname=None,
          title=None,
          plot=True,
          mode=None,
          norm=None,
          vmin=None,
          vmax=None,
          operation="mean",
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
                                             operation=operation,
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
    if dx is not None and not isinstance(dx, Quantity):
        dx *= dataset["amr"]["xyz"].unit
    if dy is not None and not isinstance(dy, Quantity):
        dy *= dataset["amr"]["xyz"].unit

    dir_vecs, origin = get_slice_direction(direction=direction,
                                           dataset=dataset,
                                           dx=dx,
                                           dy=dy,
                                           origin=origin)

    # Distance to the plane
    xyz = dataset["amr"]["xyz"] - origin
    diagonal = dataset["amr"]["dx"] * np.sqrt(dataset.meta["ndim"]) * 0.5
    dist1 = np.sum(xyz * dir_vecs[0], axis=1)
    # Create an array of indices to allow further narrowing of the selection below
    global_selection = np.arange(len(dataset["amr"]["dx"]))

    # Select cells in contact with plane
    select = np.ravel(np.where(np.abs(dist1) <= diagonal))
    global_selection = global_selection[select]

    if len(select) == 0:
        raise RuntimeError("No cells were selected to construct the plane. "
                           "The resulting figure would be empty.")

    # Project coordinates onto the plane by taking dot product with axes vectors
    coords = xyz[select]
    datax = np.inner(coords, dir_vecs[1])
    datay = np.inner(coords, dir_vecs[2])
    datadx = 0.5 * dataset["amr"]["dx"][select]

    # Get limits
    limits = {
        'xmin': np.amin(datax - datadx).values,
        'xmax': np.amax(datax + datadx).values,
        'ymin': np.amin(datay - datadx).values,
        'ymax': np.amax(datay + datadx).values
    }

    # Define slice extent
    if dx is None:
        xmin = limits['xmin']
        xmax = limits['xmax']
        ymin = limits['ymin']
        ymax = limits['ymax']
    else:
        xmin = -0.5 * dx.magnitude
        xmax = xmin + dx.magnitude
        ymin = -0.5 * dy.magnitude
        ymax = ymin + dy.magnitude
        # Limit selection further by using distance from center
        dist2 = coords - datadx * np.sqrt(dataset.meta["ndim"])
        select2 = np.ravel(
            np.where(
                np.abs(dist2.norm.values) <= max(dx.magnitude, dy.magnitude) * 0.6 *
                np.sqrt(2.0)))
        coords = coords[select2]
        datax = datax[select2]
        datay = datay[select2]
        datadx = datadx[select2]
        global_selection = global_selection[select2]

    scalar_layer = []
    to_binning = []
    for ind in range(len(to_process)):
        if to_render[ind]["mode"] in ["vec", "stream"]:
            if to_process[ind].ndim < 3:
                uv = to_process[ind].array[global_selection]
            else:
                uv = np.inner(to_process[ind].array.take(global_selection, axis=0),
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
                w = w.take(global_selection, axis=0)
            to_binning.append(apply_mask(uv[:, 0]))
            to_binning.append(apply_mask(uv[:, 1]))
            to_binning.append(apply_mask(w))
            scalar_layer.append(False)
        else:
            to_binning.append(apply_mask(to_process[ind].norm.values[global_selection]))
            scalar_layer.append(True)

    # Buffer for counts
    to_binning.append(np.ones_like(to_binning[0]))

    # Construct some bin edges
    if isinstance(resolution, int):
        resolution = {'x': resolution, 'y': resolution}
    xedges = np.linspace(xmin, xmax, resolution['x'] + 1)
    yedges = np.linspace(ymin, ymax, resolution['y'] + 1)

    # In the contour plots, x and y are the centers of the cells, instead of the edges
    xcenters = to_bin_centers(xedges)
    ycenters = to_bin_centers(yedges)

    # First histogram the cell centers into the grid bins
    binned, _, _, _ = binned_statistic_2d(x=apply_mask(datay.array),
                                          y=apply_mask(datax.array),
                                          values=to_binning,
                                          statistic="mean",
                                          bins=[yedges, xedges])

    # Next, find all the empty pixels, find the cell is lies in a apply the value
    condition = np.isnan(binned[-1])
    xgrid, ygrid = np.meshgrid(xcenters, ycenters, indexing='xy')
    # Make array of cell indices (in original array of cells) with image shape
    indices = np.zeros_like(binned[-1], dtype=int)
    # We also need a mask for pixels that find no cells (outside of the data range)
    mask = np.zeros_like(binned[-1], dtype=bool)

    pixel_positions = xgrid.reshape(xgrid.shape + (1, )) * dir_vecs[1] + ygrid.reshape(
        ygrid.shape + (1, )) * dir_vecs[2]
    # We only need to search in the cells above a certain size
    large_cells = np.ravel(
        np.where(datadx >= 0.25 * (min(xedges[1] - xedges[0], yedges[1] - yedges[0]))))
    coords = coords[large_cells]
    large_cells_dx = datadx.array[large_cells]
    global_indices = np.arange(len(datadx))[large_cells]

    # To keep memory usage down to a minimum, we process the image one column at a time
    for i in range(indices.shape[-1]):
        # We know we are looking only at a column of cells, so we make a line from the
        # two end points (x1, x2), compute distance to the line to filter out the cells
        x1 = pixel_positions[0, i, :]
        x2 = pixel_positions[-1, i, :]
        x0_minus_x1 = coords.array - x1
        x0_minus_x2 = coords.array - x2
        x2_minus_x1 = x2 - x1
        distance_to_line = np.cross(x0_minus_x1, x0_minus_x2)
        if distance_to_line.ndim > 1:
            distance_to_line = np.linalg.norm(distance_to_line, axis=-1)
        else:
            distance_to_line = np.abs(distance_to_line)
        distance_to_line /= np.linalg.norm(x2_minus_x1)
        column = np.ravel(np.where(distance_to_line <= np.sqrt(3.0) * large_cells_dx))

        if len(column) > 0:
            distance_to_cell = []
            for n, c in enumerate("xyz"[:xyz.ndim]):
                distance_to_cell.append(pixel_positions[..., i, n:n + 1] -
                                        getattr(coords, c).array[column])
            # Find the cell where the x, y, z distance is smaller than dx/2
            inds = np.logical_and.reduce(
                [np.abs(d) <= large_cells_dx[column] for d in distance_to_cell])
            index_found = inds.max(axis=-1)
            index_value = global_indices[column][inds.argmax(axis=-1)]
            indices[:, i][index_found] = index_value[index_found]
            mask[:, i][np.logical_and(~index_found, condition[:, i])] = True
        else:
            mask[:, i][condition[:, i]] = True

    mask_vec = np.broadcast_to(mask.reshape(*mask.shape, 1), mask.shape + (3, ))

    # Now we fill the arrays to be sent to the renderer, also constructing vectors
    counter = 0
    for ind in range(len(to_render)):
        binned[counter][condition] = to_binning[counter][indices][condition]
        if scalar_layer[ind]:
            to_render[ind]["data"] = ma.masked_where(mask, binned[counter], copy=False)
            counter += 1
        else:
            for j in range(counter + 1, counter + 3):
                binned[j][condition] = to_binning[j][indices][condition]
            to_render[ind]["data"] = ma.masked_where(mask_vec,
                                                     np.array([
                                                         binned[counter].T,
                                                         binned[counter + 1].T,
                                                         binned[counter + 2].T
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
            _add_scatter(datax=datax,
                         datay=datay,
                         to_scatter=to_scatter,
                         origin=origin,
                         datadx=datadx,
                         dir_vecs=dir_vecs,
                         dx=dx,
                         dy=dy,
                         ax=figure["ax"])

        figure["ax"].set_xlim(xmin, xmax)
        figure["ax"].set_ylim(ymin, ymax)

        to_return.update({"fig": figure["fig"], "ax": figure["ax"]})

    return Plot(**to_return)
