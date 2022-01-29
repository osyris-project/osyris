# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np
from numba import njit, prange


@njit(parallel=True)
def evaluate_on_grid(cell_positions_in_new_basis, cell_positions_in_original_basis,
                     cell_values, cell_sizes, grid_lower_edge_in_new_basis,
                     grid_spacing_in_new_basis, grid_positions_in_original_basis, ndim):

    nz, ny, nx = grid_positions_in_original_basis.shape[:3]
    diagonal = np.sqrt(ndim)
    out = np.full(shape=(cell_values.shape[0], nz, ny, nx),
                  fill_value=np.nan,
                  dtype=np.float64)

    ncells = cell_positions_in_new_basis.shape[0]
    for n in prange(ncells):

        half_size = cell_sizes[n] * diagonal
        ix1 = max(
            int(((cell_positions_in_new_basis[n, 0] - half_size) -
                 grid_lower_edge_in_new_basis[0]) / grid_spacing_in_new_basis[0]), 0)
        ix2 = min(
            int(((cell_positions_in_new_basis[n, 0] + half_size) -
                 grid_lower_edge_in_new_basis[0]) / grid_spacing_in_new_basis[0]) + 1,
            nx)
        iy1 = max(
            int(((cell_positions_in_new_basis[n, 1] - half_size) -
                 grid_lower_edge_in_new_basis[1]) / grid_spacing_in_new_basis[1]), 0)
        iy2 = min(
            int(((cell_positions_in_new_basis[n, 1] + half_size) -
                 grid_lower_edge_in_new_basis[1]) / grid_spacing_in_new_basis[1]) + 1,
            ny)
        iz1 = max(
            int(((cell_positions_in_new_basis[n, 2] - half_size) -
                 grid_lower_edge_in_new_basis[2]) / grid_spacing_in_new_basis[2]), 0)
        iz2 = min(
            int(((cell_positions_in_new_basis[n, 2] + half_size) -
                 grid_lower_edge_in_new_basis[2]) / grid_spacing_in_new_basis[2]) + 1,
            nz)

        for k in range(iz1, iz2):
            for j in range(iy1, iy2):
                for i in range(ix1, ix2):
                    dist = np.abs(grid_positions_in_original_basis[k, j, i, :] -
                                  cell_positions_in_original_basis[n, :])
                    if np.all(dist <= cell_sizes[n]):
                        out[:, k, j, i] = cell_values[:, n]

    return out


@njit(parallel=True)
def hist2d(x, y, values, xmin, xmax, nx, ymin, ymax, ny):

    out = np.zeros(shape=(values.shape[0], ny, nx), dtype=np.float64)
    counts = np.zeros(shape=(ny, nx), dtype=np.int64)
    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny

    for i in prange(len(x)):
        indx = int((x[i] - xmin) / dx)
        indy = int((y[i] - ymin) / dy)
        if (indx >= 0) and (indx < nx) and (indy >= 0) and (indy < ny):
            out[:, indy, indx] += values[:, i]
            counts[indy, indx] += 1

    return out, counts
