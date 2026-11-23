# SPDX-License-Identifier: BSD-3-Clause

import matplotlib.pyplot as plt
import numpy as np
import pytest

import osyris


@pytest.fixture(scope="module")
def sf_data(starformation_data_path):
    return osyris.RamsesDataset(8, path=starformation_data_path).load(cpu_list=[9])[
        "mesh"
    ]


def test_hist2d(sf_data):
    h = osyris.hist2d(sf_data["density"], sf_data["B_field"])
    assert len(h.ax.collections) > 0


def test_hist2d_lognorm(sf_data):
    h = osyris.hist2d(sf_data["density"], sf_data["B_field"], norm="log")
    assert len(h.ax.collections) > 0


def test_hist2d_logx(sf_data):
    h = osyris.hist2d(sf_data["density"], sf_data["B_field"], logx=True)
    assert len(h.ax.collections) > 0


def test_hist2d_logy(sf_data):
    h = osyris.hist2d(sf_data["density"], sf_data["B_field"], logy=True)
    assert len(h.ax.collections) > 0


def test_hist2d_loglog(sf_data):
    h = osyris.hist2d(sf_data["density"], sf_data["B_field"], loglog=True)
    assert len(h.ax.collections) > 0


def test_hist2d_cmap(sf_data):
    h = osyris.hist2d(sf_data["density"], sf_data["B_field"], cmap="magma")
    assert len(h.ax.collections) > 0


def test_hist2d_bins_int(sf_data):
    h = osyris.hist2d(sf_data["density"], sf_data["B_field"], bins=10)
    assert len(h.ax.collections) > 0


def test_hist2d_bins_list(sf_data):
    h = osyris.hist2d(sf_data["density"], sf_data["B_field"], bins=[10, 20])
    assert len(h.ax.collections) > 0


def test_hist2d_bins_array(sf_data):
    h = osyris.hist2d(
        sf_data["density"],
        sf_data["B_field"],
        norm="log",
        loglog=True,
        bins=(np.logspace(-19, -16, 301), np.logspace(-5, -3, 301)),
    )
    assert len(h.ax.collections) > 0


def test_hist2d_layer_as_color(sf_data):
    h = osyris.hist2d(
        sf_data["density"],
        sf_data["B_field"],
        sf_data.layer("mass", norm="log"),
        loglog=True,
    )
    assert "mass" in h.fig.get_axes()[1].get_ylabel()


@pytest.mark.parametrize("operation", ["mean", "sum", "max"])
def test_hist2d_other_ops(sf_data, operation):
    h = osyris.hist2d(
        sf_data["density"],
        sf_data["B_field"],
        sf_data.layer("level", operation=operation),
        loglog=True,
    )
    assert h.layers[0]["data"].max() < 10


def test_hist2d_multiple_layers(sf_data):
    h = osyris.hist2d(
        sf_data["density"],
        sf_data["B_field"],
        sf_data.layer("mass", norm="log"),  # layer 1
        sf_data.layer(
            "level",
            operation="mean",
            fmt="%i",
            mode="contour",
            colors="k",
            levels=[5, 6, 7, 8, 9],
        ),  # layer 2
        loglog=True,
    )
    assert len(h.ax.collections) == 2


def test_hist2d_subplots(sf_data):
    _, ax = plt.subplots(1, 2, figsize=(12, 4))
    osyris.hist2d(
        sf_data["density"], sf_data["B_field"], norm="log", loglog=True, ax=ax[0]
    )
    osyris.hist2d(
        sf_data["density"],
        sf_data["velocity"],
        sf_data["mass"],
        norm="log",
        loglog=True,
        cmap="magma",
        ax=ax[1],
    )
