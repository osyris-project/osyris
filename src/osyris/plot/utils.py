# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np
import numba


@numba.jit(nopython=True)
def evaluate_on_grid(cell_positions_in_new_basis, cell_positions_in_original_basis,
                     cell_values, cell_sizes, grid_lower_edge_in_new_basis,
                     grid_spacing_in_new_basis, grid_positions_in_original_basis, ndim):

    ny, nx = grid_positions_in_original_basis.shape[:2]
    diagonal = np.sqrt(ndim)
    out = np.full(shape=(cell_values.shape[0], ny, nx),
                  fill_value=np.nan,
                  dtype=np.float64)

    ncells = cell_positions_in_new_basis.shape[0]
    for n in range(ncells):

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

        for j in range(iy1, iy2):
            for i in range(ix1, ix2):
                dist = np.abs(grid_positions_in_original_basis[j, i, :] -
                              cell_positions_in_original_basis[n, :])
                if np.all(dist <= cell_sizes[n]):
                    out[:, j, i] = cell_values[:, n]

    return out
