# SPDX-License-Identifier: BSD-3-Clause

import numpy as np

import osyris


def test_hist1d(starformation_data_path):
    data = osyris.RamsesDataset(8, path=starformation_data_path).load()["mesh"]
    h = osyris.hist1d(data["velocity"].norm)
    assert len(h.ax.patches) > 0


def test_hist1d_logx(starformation_data_path):
    data = osyris.RamsesDataset(8, path=starformation_data_path).load()["mesh"]
    h = osyris.hist1d(data["velocity"].norm, logx=True)
    assert h.ax.get_xscale() == "log"
    assert h.ax.get_yscale() == "linear"


def test_hist1d_logy(starformation_data_path):
    data = osyris.RamsesDataset(8, path=starformation_data_path).load()["mesh"]
    h = osyris.hist1d(data["velocity"].norm, logy=True)
    assert h.ax.get_xscale() == "linear"
    assert h.ax.get_yscale() == "log"


def test_hist1d_loglog(starformation_data_path):
    data = osyris.RamsesDataset(8, path=starformation_data_path).load()["mesh"]
    h = osyris.hist1d(data["velocity"].norm, loglog=True)
    assert h.ax.get_xscale() == "log"
    assert h.ax.get_yscale() == "log"


def test_hist1d_custom_bins(starformation_data_path):
    data = osyris.RamsesDataset(8, path=starformation_data_path).load()["mesh"]
    h = osyris.hist1d(data["density"], logx=True, bins=np.logspace(-18.0, -13.0, 10))
    assert len(h.ax.patches) == 9


def test_hist1d_color(starformation_data_path):
    data = osyris.RamsesDataset(8, path=starformation_data_path).load()["mesh"]
    h = osyris.hist1d(data["velocity"].norm, color="red")
    for patch in h.ax.patches:
        assert patch.get_facecolor()[:3] == (1.0, 0.0, 0.0)


def test_hist1d_multiple_layers_from_layers(starformation_data_path):
    data = osyris.RamsesDataset(8, path=starformation_data_path).load()["mesh"]
    bins = np.linspace(-0.15, 0.15, 40)
    osyris.hist1d(
        data.layer("B_field", alpha=0.5).x,
        data.layer("B_field", alpha=0.5).y,
        data.layer("B_field", alpha=0.5).z,
        logy=True,
        bins=bins,
    )


def test_hist1d_multiple_layers_from_array(starformation_data_path):
    data = osyris.RamsesDataset(8, path=starformation_data_path).load()["mesh"]
    bins = np.linspace(-0.15, 0.15, 40)
    osyris.hist1d(
        data["B_field"].x,
        data["B_field"].y,
        data["B_field"].z,
        logy=True,
        bins=bins,
        alpha=0.5,
    )


def test_hist1d_multiple_layers_different_bins(starformation_data_path):
    data = osyris.RamsesDataset(8, path=starformation_data_path).load()["mesh"]
    osyris.hist1d(
        data.layer("B_field", alpha=0.5, bins=40).x,
        data.layer("B_field", alpha=0.5, bins=10).y,
        logy=True,
    )
