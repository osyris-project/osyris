import numpy as np
import numba


@numba.jit(nopython=True)
def histogram(x, y, values, sizes, xbins, ybins, ndim):
    nx = len(xbins) - 1
    ny = len(ybins) - 1
    out1 = np.zeros((values.shape[0], ny, nx), dtype=np.float64)
    counts1 = np.zeros((ny, nx), dtype=np.int64)
    out2 = np.zeros((values.shape[0], ny, nx), dtype=np.float64)
    counts2 = np.zeros((ny, nx), dtype=np.int64)
    # sizes = 0.5 * sizes
    dx = xbins[1] - xbins[0]
    dy = ybins[1] - ybins[0]
    xmin = xbins[0]
    xmax = xbins[-1]
    ymin = ybins[0]
    ymax = ybins[-1]

    diagonal = np.sqrt(ndim)

    for i in range(len(x)):
        half_size = 0.5 * sizes[i]
        # ix1 = max(int(((x[i] - half_size) - xbins[0]) / dx), 0)
        # ix2 = min(int(((x[i] + half_size) - xbins[0]) / dx), nx)
        # iy1 = max(int(((y[i] - half_size) - ybins[0]) / dy), 0)
        # iy2 = min(int(((y[i] + half_size) - ybins[0]) / dy), ny)
        ix1 = max(int(((x[i] - half_size) - xbins[0]) / dx), 0)
        ix2 = min(max(int(((x[i] + half_size) - xbins[0]) / dx), ix1 + 1), nx)
        iy1 = max(int(((y[i] - half_size) - ybins[0]) / dy), 0)
        iy2 = min(max(int(((y[i] + half_size) - ybins[0]) / dy), iy1 + 1), ny)
        if (ix2 > 0 and iy2 > 0) or (ix1 < nx and iy1 < ny):
            # out[:, iy1:iy2 + 1, ix1:ix2 + 1] += values[:, i]
            # counts[iy1:iy2 + 1, ix1:ix2 + 1] += 1.0
            out1[:, iy1:iy2, ix1:ix2] += values[:, i]
            counts1[iy1:iy2, ix1:ix2] += 1

        half_size = 0.5 * sizes[i] * diagonal
        ix1 = max(int(((x[i] - half_size) - xbins[0]) / dx), 0)
        ix2 = min(max(int(((x[i] + half_size) - xbins[0]) / dx), ix1 + 1), nx)
        iy1 = max(int(((y[i] - half_size) - ybins[0]) / dy), 0)
        iy2 = min(max(int(((y[i] + half_size) - ybins[0]) / dy), iy1 + 1), ny)
        if (ix2 > 0 and iy2 > 0) or (ix1 < nx and iy1 < ny):
            # out[:, iy1:iy2 + 1, ix1:ix2 + 1] += values[:, i]
            # counts[iy1:iy2 + 1, ix1:ix2 + 1] += 1.0
            out2[:, iy1:iy2, ix1:ix2] += values[:, i]
            counts2[iy1:iy2, ix1:ix2] += 1

    for j in range(ny):
        for i in range(nx):
            out1[:, j, i] /= counts1[j, i]
            # # out2[j, i] /= counts2[j, i]
            if counts1[j, i] == 0 and counts2[j, i] > 0:
                out1[:, j, i] = out2[:, j, i] / counts2[j, i]

    # out1 /= counts1
    # out2 /= counts2

    # indices = np.isnan(out1)

    # out1[indices] = out2[indices]

    return out1


@numba.jit(nopython=True)
def interpolate(sampling_positions, cell_positions, cell_sizes, cell_values):
    out = np.zeros((cell_values.shape[0], sampling_positions.shape[0]),
                   dtype=np.float64)
    for i in range(sampling_positions.shape[0]):
        for j in range(cell_positions.shape[0]):
            dist = np.abs(sampling_positions[i, :] - cell_positions[j, :])
            if np.all(dist <= cell_sizes[j]):
                out[:, i] = cell_values[:, j]
                break
    return out


@numba.jit(nopython=True)
def hist_and_fill(x, y, values, sizes, xbins, ybins, ndim, sampling_positions,
                  cell_positions):
    nx = len(xbins) - 1
    ny = len(ybins) - 1

    nmax = 8

    diagonal = np.sqrt(ndim)

    out = np.full(shape=(values.shape[0], ny, nx), fill_value=np.nan, dtype=np.float64)
    # counts = np.zeros((ny, nx), dtype=np.float64)
    buckets = np.zeros(shape=(ny, nx, nmax), dtype=np.int64)
    bucket_counts = np.zeros((ny, nx), dtype=np.int64)

    # sizes = 0.5 * sizes
    dx = xbins[1] - xbins[0]
    dy = ybins[1] - ybins[0]
    for n in range(len(x)):
        # indx = int((x[n] - xbins[0]) / dx)
        # indy = int((y[n] - ybins[0]) / dy)
        # if (indx >= 0 and indx < nx) and (indy >= 0 and indy < ny):
        #     out[:, indy, indx] += values[:, n]
        #     counts[indy, indx] += 1.0

        half_size = sizes[n] * diagonal

        # ix1 = max(int(((x[n] - half_size) - xbins[0]) / dx), 0)
        # ix2 = min(max(int(((x[n] + half_size) - xbins[0]) / dx), ix1 + 1), nx)
        # iy1 = max(int(((y[n] - half_size) - ybins[0]) / dy), 0)
        # iy2 = min(max(int(((y[n] + half_size) - ybins[0]) / dy), iy1 + 1), ny)

        ix1 = max(int(((x[n] - half_size) - xbins[0]) / dx), 0)
        ix2 = min(int(((x[n] + half_size) - xbins[0]) / dx) + 1, nx)
        iy1 = max(int(((y[n] - half_size) - ybins[0]) / dy), 0)
        iy2 = min(int(((y[n] + half_size) - ybins[0]) / dy) + 1, ny)

        #     # out[:, iy1:iy2 + 1, ix1:ix2 + 1] += values[:, i]
        #     # counts[iy1:iy2 + 1, ix1:ix2 + 1] += 1.0
        #     out2[:, iy1:iy2, ix1:ix2] += values[:, i]
        #     counts2[iy1:iy2, ix1:ix2] += 1

        # ix1 = int(((x[i] - half_size) - xbins[0]) / dx)
        # # ix2 = max(int(((x[i] + half_size) - xbins[0]) / dx), ix1 + 1)
        # ix2 = int(((x[i] + half_size) - xbins[0]) / dx) + 1
        # iy1 = int(((y[i] - half_size) - ybins[0]) / dy)
        # # iy2 = max(int(((y[i] + half_size) - ybins[0]) / dy), iy1 + 1)
        # iy2 = int(((y[i] + half_size) - ybins[0]) / dy) + 1
        # # out[:, iy1:iy2 + 1, ix1:ix2 + 1] += values[:, i]
        # # counts[iy1:iy2 + 1, ix1:ix2 + 1] += 1.0
        # if (ix2 > 0 and iy2 > 0) or (ix1 < nx and iy1 < ny):
        for j in range(iy1, iy2):
            for i in range(ix1, ix2):
                dist = np.abs(sampling_positions[j, i, :] - cell_positions[n, :])
                if np.all(dist <= sizes[n]):
                    out[:, j, i] = values[:, n]
                    # break

    #                 if bucket_counts[j, i] < nmax:
    #                     buckets[j, i, bucket_counts[j, i]] = n
    #                     bucket_counts[j, i] += 1
    #                 if

    #     # out[:, iy1:iy2, ix1:ix2] += values[:, i]
    #     # counts[iy1:iy2, ix1:ix2] += 1.0

    # # out /= counts

    # for j in range(ny):
    #     for i in range(nx):
    #         # if counts[j, i] == 0:
    #         for k in range(bucket_counts[j, i]):
    #             n = buckets[j, i, k]
    #             dist = np.abs(sampling_positions[j, i, :] - cell_positions[n, :])
    #             if np.all(dist <= sizes[n]):
    #                 out[:, j, i] = values[:, n]
    #                 break
    #         # else:
    #         #     out[:, j, i] /= counts[j, i]

    return out
