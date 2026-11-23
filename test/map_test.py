# SPDX-License-Identifier: BSD-3-Clause

import matplotlib.pyplot as plt
import numpy as np
import pytest

import osyris


@pytest.fixture(scope="module")
def sf_data(starformation_data_path):
    mesh = osyris.RamsesDataset(8, path=starformation_data_path).load()["mesh"]
    ind = np.argmax(mesh["density"])
    center = mesh["position"][ind]
    return mesh, center


def test_map(sf_data):
    data, center = sf_data
    m = osyris.map(
        data.layer("density"),
        norm="log",
        dx=2000 * osyris.units("au"),
        origin=center,
        direction="z",
    )
    layer = m.layers[0]
    assert layer["name"] == "density"
    assert layer["unit"] == osyris.units("g/cm^3")


@pytest.mark.parametrize("mode", ["image", "contour", "contourf"])
def test_map_modes(sf_data, mode):
    data, center = sf_data
    m = osyris.map(
        data.layer("density"),
        norm="log",
        dx=2000 * osyris.units("au"),
        origin=center,
        direction="z",
        mode=mode,
    )
    layer = m.layers[0]
    assert layer["mode"] == mode


def test_map_cmap(sf_data):
    data, center = sf_data
    osyris.map(
        data.layer("density"),
        norm="log",
        dx=2000 * osyris.units("au"),
        origin=center,
        direction="z",
        cmap="magma",
    )


def test_map_vmin_vmax(sf_data):
    data, center = sf_data
    osyris.map(
        data.layer("density"),
        norm="log",
        dx=2000 * osyris.units("au"),
        origin=center,
        direction="z",
        vmin=1.0e-16,
        vmax=1.0e-13,
    )


def test_map_resolution(sf_data):
    data, center = sf_data
    osyris.map(
        data.layer("density"),
        norm="log",
        dx=2000 * osyris.units("au"),
        origin=center,
        direction="z",
        resolution=64,
    )


def test_map_vectors(sf_data):
    data, center = sf_data
    osyris.map(
        data.layer("velocity"),
        mode="vec",
        dx=2000 * osyris.units("au"),
        origin=center,
        direction="z",
        color="k",
    )


def test_map_vectors_with_color_variable(sf_data):
    data, center = sf_data
    osyris.map(
        data.layer("velocity"),
        mode="vec",
        dx=2000 * osyris.units("au"),
        origin=center,
        direction="z",
        color=data["velocity"],
    )


def test_map_streamlines(sf_data):
    data, center = sf_data
    osyris.map(
        data.layer("velocity"),
        mode="stream",
        dx=2000 * osyris.units("au"),
        origin=center,
        direction="z",
        color="g",
    )


def test_map_streamlines_with_color_variable(sf_data):
    data, center = sf_data
    osyris.map(
        data.layer("velocity"),
        mode="stream",
        dx=2000 * osyris.units("au"),
        origin=center,
        direction="z",
        color=data["velocity"],
    )


@pytest.mark.skip(reason="not implemented yet")
def test_map_lic(sf_data):
    pass


@pytest.mark.skip(reason="not implemented yet")
def test_map_lic_with_color_variable(sf_data):
    pass


def test_map_density_and_vectors(sf_data):
    data, center = sf_data
    osyris.map(
        data.layer("density", norm="log"),  # layer 1
        data.layer("velocity", mode="vec"),  # layer 2
        dx=2000 * osyris.units("au"),
        origin=center,
        direction="z",
    )


def test_map_temperature_and_amr_levels(sf_data):
    data, center = sf_data
    osyris.map(
        data.layer(
            "temperature",
            norm="log",
            mode="contourf",
            levels=np.logspace(0.9, 2, 11),
            cmap="hot",
        ),  # layer 1
        data.layer(
            "level",
            mode="contour",
            colors="w",
            levels=[6, 7, 8],
            fmt="%i",
        ),  # layer 2
        dx=2000 * osyris.units("au"),
        origin=center,
        direction="z",
    )


@pytest.mark.parametrize("direction", ["x", "y", "z"])
def test_map_direction(sf_data, direction):
    data, center = sf_data
    osyris.map(
        data.layer("density", norm="log"),
        data.layer("velocity", mode="vec"),
        dx=1000 * osyris.units("au"),
        origin=center,
        direction=direction,
    )


def test_map_direction_from_normal_unit_vector(sf_data):
    data, center = sf_data
    v = osyris.Vector(-1, 1, 3)
    osyris.map(
        data.layer("density", norm="log"),
        data.layer("velocity", mode="vec"),
        dx=500 * osyris.units("au"),
        origin=center,
        direction=v / v.norm,
        cmap="magma",
    )


def test_map_direction_from_normal_vector_unnormalized(sf_data):
    data, center = sf_data
    osyris.map(
        data.layer("density", norm="log"),
        data.layer("velocity", mode="vec"),
        dx=500 * osyris.units("au"),
        origin=center,
        direction=osyris.Vector(-1, 1, 3),
        cmap="magma",
    )


def test_map_top_direction(sf_data):
    data, center = sf_data
    osyris.map(
        data.layer("density", norm="log"),
        data.layer("velocity", mode="vec"),
        dx=500 * osyris.units("au"),
        origin=center,
        direction="top",
    )


def test_map_side_direction(sf_data):
    data, center = sf_data
    osyris.map(
        data.layer("density", norm="log"),
        data.layer("velocity", mode="vec"),
        dx=500 * osyris.units("au"),
        origin=center,
        direction="side",
    )


def test_map_off_center(sf_data):
    data, center = sf_data
    osyris.map(
        data.layer("density", norm="log"),
        data.layer("velocity", mode="vec"),
        dx=500 * osyris.units("au"),
        origin=center + osyris.Vector(0, 0, 20, unit="au"),
        direction="z",
    )


def test_map_disk_cells_only(sf_data):
    data, center = sf_data
    selection = (data["density"] >= (5.0e-14 * osyris.units("g/cm**3"))) & (
        data["density"] <= (5.0e-12 * osyris.units("g/cm**3"))
    )
    osyris.map(
        data[selection].layer("density"),
        dx=500 * osyris.units("au"),
        norm="log",
        origin=center,
    )


def test_map_subplots(sf_data):
    mesh, center = sf_data
    _, ax = plt.subplots(nrows=1, ncols=3, figsize=(12, 4))
    # Define region to plot
    dx = 2000.0 * osyris.units("au")

    osyris.map(
        mesh.layer("density", norm="log"),
        mesh.layer("velocity", mode="vec"),
        dx=dx,
        origin=center,
        direction="z",
        ax=ax[0],
    )

    osyris.map(
        mesh.layer(
            "temperature",
            norm="log",
            mode="contourf",
            levels=np.logspace(0.9, 2, 11),
            cmap="hot",
        ),
        mesh.layer(
            "level",
            mode="contour",
            colors="w",
            levels=[6, 7, 8],
            fmt="%i",
        ),
        dx=dx,
        origin=center,
        direction="z",
        ax=ax[1],
    )

    osyris.map(
        mesh.layer("density", norm="log", cmap="magma"),
        mesh.layer("velocity", mode="stream", color="w"),
        dx=0.25 * dx,
        origin=center,
        direction="side",
        ax=ax[2],
    )
