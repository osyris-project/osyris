# SPDX-License-Identifier: BSD-3-Clause

from typing import Union

import numpy as np
import numpy.ma as ma
from numba import njit, prange
from pint import Quantity

from ..core import Array, Layer, Plot, Vector, VectorBasis
from ..core.tools import apply_mask
from .direction import get_direction
from .parser import get_norm, parse_layer
from .render import render
from .scatter import scatter


def _evaluate_on_grid(
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
    nx,
    ny,
    nz,
    ndim,
    ux,
    uy,
    uz,
    vx,
    vy,
    vz,
    nx_vec,
    ny_vec,
    nz_vec,
):
    inv_dx = 1.0 / grid_spacing_in_new_basis_x
    inv_dy = 1.0 / grid_spacing_in_new_basis_y
    inv_dz = 1.0 / grid_spacing_in_new_basis_z

    out = np.full(
        shape=(cell_values.shape[0], nz, ny, nx), fill_value=np.nan, dtype=np.float64
    )

    ncells = len(cell_positions_in_new_basis_x)
    diagonal = np.sqrt(ndim)

    has_y = cell_positions_in_original_basis_y is not None
    has_z = cell_positions_in_original_basis_z is not None

    for n in prange(ncells):
        half_size = cell_sizes[n] * diagonal
        current_val = cell_values[:, n]
        current_size = cell_sizes[n]

        pos_orig_x = cell_positions_in_original_basis_x[n]
        pos_orig_y = cell_positions_in_original_basis_y[n] if has_y else 0.0
        pos_orig_z = cell_positions_in_original_basis_z[n] if has_z else 0.0

        rel_x = cell_positions_in_new_basis_x[n] - grid_lower_edge_in_new_basis_x
        rel_y = cell_positions_in_new_basis_y[n] - grid_lower_edge_in_new_basis_y
        rel_z = cell_positions_in_new_basis_z[n] - grid_lower_edge_in_new_basis_z

        ix1 = int((rel_x - half_size) * inv_dx)
        ix2 = int((rel_x + half_size) * inv_dx) + 1
        iy1 = int((rel_y - half_size) * inv_dy)
        iy2 = int((rel_y + half_size) * inv_dy) + 1
        iz1 = int((rel_z - half_size) * inv_dz)
        iz2 = int((rel_z + half_size) * inv_dz) + 1

        ix1 = max(ix1, 0)
        ix2 = min(ix2, nx)
        iy1 = max(iy1, 0)
        iy2 = min(iy2, ny)
        iz1 = max(iz1, 0)
        iz2 = min(iz2, nz)

        for k in range(iz1, iz2):
            # get z coords
            z_map = (
                grid_lower_edge_in_new_basis_z + (k + 0.5) * grid_spacing_in_new_basis_z
            )

            pz_x = z_map * nx_vec
            pz_y = z_map * ny_vec
            pz_z = z_map * nz_vec

            for j in range(iy1, iy2):
                # get y coords
                y_map = (
                    grid_lower_edge_in_new_basis_y
                    + (j + 0.5) * grid_spacing_in_new_basis_y
                )

                py_x = y_map * vx
                py_y = y_map * vy
                py_z = y_map * vz

                pyz_x = py_x + pz_x
                pyz_y = py_y + pz_y
                pyz_z = py_z + pz_z

                for i in range(ix1, ix2):
                    # get x coords
                    x_map = (
                        grid_lower_edge_in_new_basis_x
                        + (i + 0.5) * grid_spacing_in_new_basis_x
                    )

                    grid_x = x_map * ux + pyz_x

                    dist_x = grid_x - pos_orig_x
                    if np.abs(dist_x) > current_size:
                        continue

                    if has_y:
                        grid_y = x_map * uy + pyz_y
                        dist_y = grid_y - pos_orig_y
                        if np.abs(dist_y) > current_size:
                            continue

                    if has_z:
                        grid_z = x_map * uz + pyz_z
                        dist_z = grid_z - pos_orig_z
                        if np.abs(dist_z) > current_size:
                            continue

                    out[:, k, j, i] = current_val

    return out


# compile two versions of the above function
_evaluate_on_grid_fastmath = njit(parallel=True, fastmath=True)(_evaluate_on_grid)
_evaluate_on_grid_precise = njit(parallel=True, fastmath=False)(_evaluate_on_grid)


def _add_scatter(to_scatter, origin, basis, dx, dy, ax, map_unit):
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
    dist_to_plane = xyz.dot(basis.n)
    global_selection = np.arange(len(to_scatter[0]["data"]))
    select = (np.abs(dist_to_plane) <= radius).values
    global_selection = global_selection[select]
    if len(select) > 0:
        # Project coordinates onto the plane by taking dot product with axes vectors
        coords = xyz[select]
        datax = coords.dot(basis.u)
        datay = coords.dot(basis.v)

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
        datax.name = basis.u.name
        datay.name = basis.v.name
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
    fastmath: bool = False,
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

    to_process = []
    to_render = []
    to_scatter = []
    for layer in layers:
        if not isinstance(layer, Layer):
            raise TypeError(f"Expected Layer object, got {type(layer)} instead. ")
        layer = parse_layer(
            layer,
            mode=mode,
            operation=operation,
            norm=norm,
            vmin=vmin,
            vmax=vmax,
            **kwargs,
        )
        layer.kwargs.update(
            norm=get_norm(norm=layer.norm, vmin=layer.vmin, vmax=layer.vmax)
        )
        if layer.mode == "scatter":
            to_scatter.append({"data": layer.data, "params": layer.kwargs})
        else:
            to_process.append(layer.data)
            to_render.append(
                {
                    "mode": layer.mode,
                    "params": layer.kwargs,
                    "unit": layer.data.unit,
                    "name": layer.data.name,
                }
            )

    position = layers[0]["position"]
    cell_size = layers[0]["dx"]
    ndim = position.nvec

    thick = dz is not None

    spatial_unit = position.unit
    map_unit = spatial_unit

    # Set window size
    if dx is not None:
        map_unit = dx.units
        dx = dx.to(spatial_unit)
    dy = dx if dy is None else dy.to(spatial_unit)
    dz = dx if dz is None else dz.to(spatial_unit)

    if origin is None:
        origin = Vector(*([0] * ndim), unit=spatial_unit)

    if ndim < 3:
        basis = VectorBasis(
            n=Vector(0, 0, name="z"), u=Vector(1, 0, name="x"), v=Vector(0, 1, name="y")
        )
    else:
        basis = get_direction(
            direction=direction,
            data=layers[0],
            dx=dx,
            dy=dy,
            origin=origin,
        )

    # Distance to the plane
    diagonal = np.sqrt(ndim)
    xyz = position - origin
    selection_distance = 0.5 * diagonal * (dz if thick else cell_size)

    normal = basis.n
    vec_u = basis.u
    vec_v = basis.v

    dist_to_plane = xyz.dot(normal)
    # Create an array of indices to allow further narrowing of the selection below
    global_indices = np.arange(len(cell_size))
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
            - 0.5 * cell_size[indices_close_to_plane] * diagonal
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
    datadx = cell_size[indices_close_to_plane] * 0.5

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

    nx_pix = int(resolution["x"])
    ny_pix = int(resolution["y"])

    if thick:
        if "z" not in resolution:
            resolution["z"] = round((zmax - zmin) / (0.5 * (xspacing + yspacing)))
        zspacing = (zmax - zmin) / resolution["z"]
        nz_pix = int(resolution["z"])
    else:
        zmin = -0.5
        zspacing = 1.0
        nz_pix = 1

    # flatten vectors for Numba
    u_vals = np.array(
        [vec_u.x.values, vec_u.y.values, vec_u.z.values if vec_u.z is not None else 0.0]
    )
    v_vals = np.array(
        [vec_v.x.values, vec_v.y.values, vec_v.z.values if vec_v.z is not None else 0.0]
    )
    n_vals = np.array(
        [
            normal.x.values,
            normal.y.values,
            normal.z.values if normal.z is not None else 0.0,
        ]
    )

    cell_values_arr = np.array(to_binning)

    grid_eval = _evaluate_on_grid_fast if fastmath else _evaluate_on_grid_precise

    binned = grid_eval(
        cell_positions_in_new_basis_x=apply_mask(datax.values),
        cell_positions_in_new_basis_y=apply_mask(datay.values),
        cell_positions_in_new_basis_z=apply_mask(dataz.values),
        cell_positions_in_original_basis_x=coords.x.values,
        cell_positions_in_original_basis_y=(
            coords.y.values if coords.y is not None else None
        ),
        cell_positions_in_original_basis_z=(
            coords.z.values if coords.z is not None else None
        ),
        cell_values=cell_values_arr,
        cell_sizes=datadx.values,
        grid_lower_edge_in_new_basis_x=xmin,
        grid_lower_edge_in_new_basis_y=ymin,
        grid_lower_edge_in_new_basis_z=zmin,
        grid_spacing_in_new_basis_x=xspacing,
        grid_spacing_in_new_basis_y=yspacing,
        grid_spacing_in_new_basis_z=zspacing,
        nx=nx_pix,
        ny=ny_pix,
        nz=nz_pix,
        ndim=ndim,
        # Basis vectors
        ux=u_vals[0],
        uy=u_vals[1],
        uz=u_vals[2],
        vx=v_vals[0],
        vy=v_vals[1],
        vz=v_vals[2],
        nx_vec=n_vals[0],
        ny_vec=n_vals[1],
        nz_vec=n_vals[2],
    )

    xcenters = np.linspace(xmin + 0.5 * xspacing, xmax - 0.5 * xspacing, nx_pix)
    ycenters = np.linspace(ymin + 0.5 * yspacing, ymax - 0.5 * yspacing, ny_pix)
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
        figure["ax"].set_xlabel(Array(values=0, unit=map_unit, name=basis.u.name).label)
        figure["ax"].set_ylabel(Array(values=0, unit=map_unit, name=basis.v.name).label)
        if ax is None:
            figure["ax"].set_aspect("equal")

        # Add scatter layer
        if len(to_scatter) > 0:
            _add_scatter(
                to_scatter=to_scatter,
                origin=origin,
                basis=basis,
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
        figure["ax"].set_title(title)

        to_return.update({"fig": figure["fig"], "ax": figure["ax"]})

    return Plot(**to_return)
