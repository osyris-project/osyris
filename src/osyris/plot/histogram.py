import numpy as np
import numba


@numba.jit(nopython=True)
def histogram(x, y, values, sizes, xbins, ybins):
    out = np.zeros((values.shape[0], len(ybins) - 1, len(xbins) - 1), dtype=np.float64)
    counts = np.zeros((len(ybins) - 1, len(xbins) - 1), dtype=np.float64)
    # sizes = 0.5 * sizes
    dx = xbins[1] - xbins[0]
    dy = ybins[1] - ybins[0]
    for i in range(len(x)):
        half_size = 0.5 * sizes[i]
        ix1 = int(((x[i] - half_size) - xbins[0]) / dx)
        ix2 = max(int(((x[i] + half_size) - xbins[0]) / dx), ix1 + 1)
        iy1 = int(((y[i] - half_size) - ybins[0]) / dy)
        iy2 = max(int(((y[i] + half_size) - ybins[0]) / dy), iy1 + 1)
        # out[:, iy1:iy2 + 1, ix1:ix2 + 1] += values[:, i]
        # counts[iy1:iy2 + 1, ix1:ix2 + 1] += 1.0
        out[:, iy1:iy2, ix1:ix2] += values[:, i]
        counts[iy1:iy2, ix1:ix2] += 1.0
    out /= counts
    return out
