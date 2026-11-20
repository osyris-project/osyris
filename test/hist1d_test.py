# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
import pytest

import osyris


@pytest.fixture(scope="module")
def sf_data(starformation_data_path):
    return osyris.RamsesDataset(8, path=starformation_data_path).load(cpu_list=[9])[
        "mesh"
    ]


def test_hist1d(sf_data):
    h = osyris.hist1d(sf_data["velocity"].norm)
    assert len(h.ax.patches) > 0


def test_hist1d_logx(sf_data):
    h = osyris.hist1d(sf_data["velocity"].norm, logx=True)
    assert h.ax.get_xscale() == "log"
    assert h.ax.get_yscale() == "linear"


def test_hist1d_logy(sf_data):
    h = osyris.hist1d(sf_data["velocity"].norm, logy=True)
    assert h.ax.get_xscale() == "linear"
    assert h.ax.get_yscale() == "log"


def test_hist1d_loglog(sf_data):
    h = osyris.hist1d(sf_data["velocity"].norm, loglog=True)
    assert h.ax.get_xscale() == "log"
    assert h.ax.get_yscale() == "log"


def test_hist1d_custom_bins(sf_data):
    h = osyris.hist1d(sf_data["density"], logx=True, bins=np.logspace(-18.0, -13.0, 10))
    assert len(h.ax.patches) == 9


def test_hist1d_color(sf_data):
    h = osyris.hist1d(sf_data["velocity"].norm, color="red")
    for patch in h.ax.patches:
        assert patch.get_facecolor()[:3] == (1.0, 0.0, 0.0)


def test_hist1d_multiple_layers_from_layers(sf_data):
    bins = np.linspace(-0.15, 0.15, 40)
    osyris.hist1d(
        sf_data.layer("B_field", alpha=0.5).x,
        sf_data.layer("B_field", alpha=0.5).y,
        sf_data.layer("B_field", alpha=0.5).z,
        logy=True,
        bins=bins,
    )


def test_hist1d_multiple_layers_from_array(sf_data):
    bins = np.linspace(-0.15, 0.15, 40)
    osyris.hist1d(
        sf_data["B_field"].x,
        sf_data["B_field"].y,
        sf_data["B_field"].z,
        logy=True,
        bins=bins,
        alpha=0.5,
    )


def test_hist1d_multiple_layers_different_bins(sf_data):
    osyris.hist1d(
        sf_data.layer("B_field", alpha=0.5, bins=40).x,
        sf_data.layer("B_field", alpha=0.5, bins=10).y,
        logy=True,
    )
