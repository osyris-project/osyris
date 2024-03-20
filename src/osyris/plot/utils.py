# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)

import numpy as np
from numba import njit, prange


@njit(parallel=True)
def evaluate_on_grid(
    cell_positions_in_new_basis_x,
    cell_positions_in_new_basis_y,
    cell_positions_in_new_basis_z,
    cell_positions_in_original_basis_x,
    cell_positions_in_original_basis_y,
    cell_positions_in_original_basis_z,
    cell_values,
    cell_sizes,
    grid_lower_edge_in_new_basis_x,
    grid_lower_edge_in_new_basis_y,
    grid_lower_edge_in_new_basis_z,
    grid_spacing_in_new_basis_x,
    grid_spacing_in_new_basis_y,
    grid_spacing_in_new_basis_z,
    grid_positions_in_original_basis,
    ndim,
):
    nz, ny, nx = grid_positions_in_original_basis.shape[:3]
    diagonal = np.sqrt(ndim)
    out = np.full(
        shape=(cell_values.shape[0], nz, ny, nx), fill_value=np.nan, dtype=np.float64
    )

    ncells = len(cell_positions_in_new_basis_x)
    for n in prange(ncells):
        half_size = cell_sizes[n] * diagonal
        ix1 = max(
            int(
                (
                    (cell_positions_in_new_basis_x[n] - half_size)
                    - grid_lower_edge_in_new_basis_x
                )
                / grid_spacing_in_new_basis_x
            ),
            0,
        )
        ix2 = min(
            int(
                (
                    (cell_positions_in_new_basis_x[n] + half_size)
                    - grid_lower_edge_in_new_basis_x
                )
                / grid_spacing_in_new_basis_x
            )
            + 1,
            nx,
        )
        iy1 = max(
            int(
                (
                    (cell_positions_in_new_basis_y[n] - half_size)
                    - grid_lower_edge_in_new_basis_y
                )
                / grid_spacing_in_new_basis_y
            ),
            0,
        )
        iy2 = min(
            int(
                (
                    (cell_positions_in_new_basis_y[n] + half_size)
                    - grid_lower_edge_in_new_basis_y
                )
                / grid_spacing_in_new_basis_y
            )
            + 1,
            ny,
        )
        iz1 = max(
            int(
                (
                    (cell_positions_in_new_basis_z[n] - half_size)
                    - grid_lower_edge_in_new_basis_z
                )
                / grid_spacing_in_new_basis_z
            ),
            0,
        )
        iz2 = min(
            int(
                (
                    (cell_positions_in_new_basis_z[n] + half_size)
                    - grid_lower_edge_in_new_basis_z
                )
                / grid_spacing_in_new_basis_z
            )
            + 1,
            nz,
        )

        for k in range(iz1, iz2):
            for j in range(iy1, iy2):
                for i in range(ix1, ix2):
                    ok_x = (
                        np.abs(
                            grid_positions_in_original_basis[k, j, i, 0]
                            - cell_positions_in_original_basis_x[n]
                        )
                        <= cell_sizes[n]
                    )
                    ok_y = True
                    if cell_positions_in_original_basis_y is not None:
                        ok_y = (
                            np.abs(
                                grid_positions_in_original_basis[k, j, i, 1]
                                - cell_positions_in_original_basis_y[n]
                            )
                            <= cell_sizes[n]
                        )
                    ok_z = True
                    if cell_positions_in_original_basis_z is not None:
                        ok_z = (
                            np.abs(
                                grid_positions_in_original_basis[k, j, i, 2]
                                - cell_positions_in_original_basis_z[n]
                            )
                            <= cell_sizes[n]
                        )

                    if ok_x and ok_y and ok_z:
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
