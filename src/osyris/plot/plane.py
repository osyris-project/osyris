# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np
import numpy.ma as ma
from pint.quantity import Quantity
from .slice import get_slice_direction
from .render import render
from .parser import parse_layer
from ..core import Plot, Array
from ..core.tools import to_bin_centers, apply_mask
from scipy.stats import binned_statistic_2d
from scipy.ndimage import distance_transform_edt


def plane(*layers,
          direction="z",
          dx=None,
          dy=None,
          fname=None,
          title=None,
          sinks=True,
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
    for layer in layers:
        data, settings, params = parse_layer(layer,
                                             mode=mode,
                                             norm=norm,
                                             vmin=vmin,
                                             vmax=vmax,
                                             operation=operation,
                                             **kwargs)

        to_process.append(data)
        to_render.append({
            "mode": settings["mode"],
            "params": params,
            "unit": data.unit.units,
            "name": data.name
        })
        operations.append(settings["operation"])

    dataset = to_process[0].parent

    # Set window size
    if dy is None:
        dy = dx
    if dx is not None and not isinstance(dx, Quantity):
        dx *= dataset["xyz"].unit
    if dy is not None and not isinstance(dy, Quantity):
        dy *= dataset["xyz"].unit

    dir_vecs, origin = get_slice_direction(direction=direction,
                                           dataset=dataset,
                                           dx=dx,
                                           dy=dy,
                                           origin=origin)

    # Distance to the plane
    xyz = dataset["xyz"] - origin
    diagonal = dataset["dx"] * np.sqrt(dataset.meta["ndim"]) * 0.5
    dist1 = np.sum(xyz * dir_vecs[0], axis=1)

    # Select cells in contact with plane
    cells_in_plane = np.abs(dist1) <= diagonal
    # Project coordinates onto the plane by taking dot product with axes
    # vectors
    select = np.ravel(np.where(cells_in_plane))

    if len(select) == 0:
        raise RuntimeError("No cells were selected to construct the plane. "
                           "The resulting figure would be empty.")
    coords = xyz[select]
    datax = np.inner(coords, dir_vecs[1])
    datay = np.inner(coords, dir_vecs[2])
    datadx = 0.5 * dataset["dx"][select]

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
        dist2 = xyz - diagonal
        select = np.ravel(
            np.where(
                np.logical_and(
                    cells_in_plane,
                    np.abs(dist2.norm) <=
                    max(dx.magnitude, dy.magnitude) * 0.5 * np.sqrt(2.0))))
        coords = xyz[select]
        datax = np.inner(coords, dir_vecs[1])
        datay = np.inner(coords, dir_vecs[2])
        datadx = 0.5 * dataset["dx"][select]

    scalar_layer = []
    to_binning = []
    for ind in range(len(to_process)):
        if to_render[ind]["mode"] in ["vec", "stream"]:
            if to_process[ind].ndim < 3:
                uv = to_process[ind].array[select]
            else:
                uv = np.inner(to_process[ind].array.take(select, axis=0), dir_vecs[1:])
            w = None
            if "color" in to_render[ind]["params"]:
                if isinstance(to_render[ind]["params"]["color"], Array):
                    w = to_render[ind]["params"]["color"].norm
                elif isinstance(to_render[ind]["params"]["color"], np.ndarray):
                    w = to_render[ind]["params"]["color"]
            if w is None:
                w = np.linalg.norm(uv, axis=1)
            else:
                w = w.take(select, axis=0)
            to_binning.append(apply_mask(uv[:, 0]))
            to_binning.append(apply_mask(uv[:, 1]))
            to_binning.append(apply_mask(w))
            scalar_layer.append(False)
        else:
            to_binning.append(apply_mask(to_process[ind].array[select]))
            scalar_layer.append(True)

    # Buffer for counts
    to_binning.append(np.ones_like(to_binning[0]))

    # Construct some bin edges
    if isinstance(resolution, int):
        resolution = {'x': resolution, 'y': resolution}
    xedges = np.linspace(xmin, xmax, resolution['x'] + 1)
    yedges = np.linspace(ymin, ymax, resolution['y'] + 1)

    # In the contour plots, x and y are the centers of the cells, instead of
    # the edges.
    xcenters = to_bin_centers(xedges)
    ycenters = to_bin_centers(yedges)

    binned, _, _, _ = binned_statistic_2d(x=apply_mask(datay.array),
                                          y=apply_mask(datax.array),
                                          values=to_binning,
                                          statistic="mean",
                                          bins=[yedges, xedges])

    # Use Scipy's distance transform to fill blanks with nearest neighbours
    condition = np.isnan(binned[-1])
    print(len(binned[-1]), np.sum(condition))
    nempty = np.sum(condition)

    xgrid, ygrid = np.meshgrid(xcenters, ycenters, indexing='xy')
    # ygrid, xgrid = np.meshgrid(ycenters, xcenters, indexing='xy')
    # xgrid = Array(values=xgrid, unit=xyz.unit)
    # ygrid = Array(values=ygrid, unit=xyz.unit)
    print('xgrid.shape', xgrid.shape)

    indices = np.zeros_like(binned[-1], dtype=int)
    mask = np.zeros_like(binned[-1], dtype=bool)

    pos = (xgrid.reshape(xgrid.shape + (1, )) * dir_vecs[1] +
           ygrid.reshape(ygrid.shape + (1, )) * dir_vecs[2]) + origin.array
    print('dir_vecs', dir_vecs[1], dir_vecs[2])
    print('pos', pos.shape)
    print(pos)
    print('binned', binned[-1].shape)
    fil = np.ravel(np.where(datadx >= 0.5 * (xedges[1] - xedges[0])))
    print(len(fil))

    xcoords = coords.x.array[fil]
    ycoords = coords.y.array[fil]
    filtered_dx = datadx.array[fil]
    global_indices = np.arange(len(datadx))[fil]

    times = [0, 0, 0, 0]
    import time
    print('indices.shape', indices.shape)

    for i in range(indices.shape[-1]):
        # Make a line from the two end points and compute distance to the line to filter out all points
        # first and last points in row
        x1 = pos[0, i, :]
        x2 = pos[-1, i, :]
        x0_minus_x1 = coords.array[fil] - x1
        x0_minus_x2 = coords.array[fil] - x2
        x2_minus_x1 = x2 - x1
        # print('x1', x1, x1.shape, 'x2', x2, x2.shape)
        # print(coords.array[fil].shape)
        # print('x0_minus_x1.shape', x0_minus_x1.shape)
        # print('x0_minus_x2.shape', x0_minus_x2.shape)
        distance_to_line = np.cross(x0_minus_x1, x0_minus_x2)
        if distance_to_line.ndim > 1:
            distance_to_line = np.linalg.norm(distance_to_line, axis=-1)
        else:
            distance_to_line = np.abs(distance_to_line)
        distance_to_line /= np.linalg.norm(x2_minus_x1)
        # print('cross', np.cross(x0_minus_x1, x0_minus_x2).shape)
        # print('outer', np.outer(x0_minus_x1, x0_minus_x2).shape)
        # print(distance_to_line.shape)
        # print('min/max', distance_to_line.min(), distance_to_line.max())
        row = np.ravel(np.where(distance_to_line < np.sqrt(3.0) * filtered_dx))
        # row = np.ravel(np.where(np.abs(ycoords - ycenters[i]) < filtered_dx))
        # print(i, len(row))
        # print(xcoords[row].min(), xcoords[row].max())
        start = time.time()
        distx = pos[..., i, 0:1] - xcoords[row]  #.reshape((len(coords.x.array), 1))
        disty = pos[..., i, 1:2] - ycoords[row]  #.reshape((len(coords.y.array), 1))
        end = time.time()
        times[0] += end - start
        start = end
        # print(distx.shape)

        # ind = np.where(
        #     np.logical_and(
        #         np.abs(distx) <= datadx.array[fil],
        #         np.abs(disty) <= datadx.array[fil]))

        ind = np.logical_and(
            np.abs(distx) <= filtered_dx[row],
            np.abs(disty) <= filtered_dx[row])
        # print('ind.shape', ind.shape)
        end = time.time()
        times[1] += end - start
        start = end
        # print(ind.shape)
        # print(ind.max(axis=-1))
        index_found = ind.max(axis=-1)
        index_value = global_indices[row][ind.argmax(axis=-1)]
        end = time.time()
        times[2] += end - start
        start = end

        cond = index_found == True
        # print('cond', cond)
        # print('index_found', index_found)
        # print('index_value', index_value)
        indices[:, i][cond] = index_value[cond]
        # mask[:, i][np.logical_and(~cond, condition[:, i])] = True
        end = time.time()
        times[3] += end - start
        start = end
        # break

    print('times', times)
    # print(indices.shape)
    # if len(ind) > 0:
    #     indices[index] = ind[0]
    # else:
    #     # print('masking', index)
    #     mask[index] = True

    # for index, x in np.ndenumerate(condition):
    #     if x:
    #         # print(xgrid[index] * dir_vecs[1])
    #         # print(ygrid[index] * dir_vecs[2])
    #         pos = (xgrid[index] * dir_vecs[1] + ygrid[index] * dir_vecs[2]) + origin
    #         # print(pos)
    #         # print(pos[0])
    #         # print(pos.array)
    #         distx = pos.x - coords[..., 0]
    #         disty = pos.y - coords[..., 1]
    #         # print(distx.shape, disty.shape, datadx.shape)
    #         # distz = pos.z - coords[..., 2]
    #         # print(distx)
    #         # ind = np.where(
    #         #     np.logical_and(
    #         #         np.abs(distx) <= datadx,
    #         #         np.logical_and(np.abs(disty) <= datadx,
    #         #                        np.abs(distz) <= datadx)))
    #         ind = np.ravel(
    #             np.where(
    #                 np.logical_and(np.abs(distx) <= datadx,
    #                                np.abs(disty) <= datadx)))
    #         # print(ind)
    #         if len(ind) > 0:
    #             indices[index] = ind[0]
    #         else:
    #             # print('masking', index)
    #             mask[index] = True

    # for i in range(len(binned) - 1):
    #     binned[i][condition] = to_binning[i][indices][condition]
    #     binned[i] = np.ma.masked_where(mask, binned[i], copy=False)
    #     print(i, binned[i])
    #     print(np.ma.masked_where(mask, binned[i], copy=False))

    # transform = tuple(
    #     distance_transform_edt(condition, return_distances=False, return_indices=True))

    # # Define a mask to mask aread outside of the domain range, which have
    # # been filled by the previous transform step
    # xx = np.broadcast_to(xcenters, [resolution, resolution])
    # yy = np.broadcast_to(ycenters, [resolution, resolution]).T
    # mask = np.logical_or.reduce((xx < limits['xmin'], xx > limits['xmax'],
    #                              yy < limits['ymin'], yy > limits['ymax']))
    mask_vec = np.broadcast_to(mask.reshape(*mask.shape, 1), mask.shape + (3, ))
    # print('=====')
    # print(binned[0])
    # print(mask)

    counter = 0
    for ind in range(len(to_render)):
        binned[ind][condition] = to_binning[ind][indices][condition]
        if scalar_layer[ind]:
            to_render[ind]["data"] = ma.masked_where(mask, binned[counter], copy=False)
            # to_render[ind]["data"] = binned[counter]
            counter += 1
        else:
            to_render[ind]["data"] = ma.masked_where(mask_vec,
                                                     np.array([
                                                         binned[counter].T,
                                                         binned[counter + 1].T,
                                                         binned[counter + 2].T
                                                     ]).T,
                                                     copy=False)
            # to_render[ind]["data"] = np.array(
            #     [binned[counter].T, binned[counter + 1].T, binned[counter + 2].T]).T

            counter += 3

    to_return = {"x": xcenters, "y": ycenters, "layers": to_render}
    if plot:
        # Render the map
        figure = render(x=xcenters, y=ycenters, data=to_render, ax=ax)
        figure["ax"].set_xlabel(dataset["xyz"].x.label)
        figure["ax"].set_ylabel(dataset["xyz"].y.label)
        # if ax is None:
        #     figure["ax"].set_aspect("equal")
        to_return.update({"fig": figure["fig"], "ax": figure["ax"]})

    return Plot(**to_return)
