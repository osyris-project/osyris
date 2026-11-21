# SPDX-License-Identifier: BSD-3-Clause

import matplotlib.pyplot as plt
import numpy as np
import pytest

import osyris


@pytest.fixture(scope="module")
def data1(sod_data_path):
    return osyris.RamsesDataset(1, path=sod_data_path).load()["mesh"]


@pytest.fixture(scope="module")
def data2(sod_data_path):
    return osyris.RamsesDataset(2, path=sod_data_path).load()["mesh"]


def test_plot_line_profile(data2):
    f = osyris.plot(data2["position_x"], data2["density"])
    layer = f.layers[0]
    assert np.array_equal(layer["x"].values, data2["position_x"].values)
    assert np.array_equal(layer["y"].values, data2["density"].values)
    (line,) = f.ax.lines
    inds = np.argsort(data2["position_x"].values)
    assert np.array_equal(line.get_xdata(), data2["position_x"].values[inds])
    assert np.array_equal(line.get_ydata(), data2["density"].values[inds])
    assert f.ax.get_xlabel() == "position_x [cm]"
    assert f.ax.get_ylabel() == "density [g / cm ** 3]"


def test_plot_line_profile_no_x(data2):
    f = osyris.plot(data2["density"])
    layer = f.layers[0]
    arange = np.arange(len(data2["density"].values))
    assert np.array_equal(layer["x"].values, arange)
    assert np.array_equal(layer["y"].values, data2["density"].values)
    (line,) = f.ax.lines
    assert np.array_equal(line.get_xdata(), arange)
    assert np.array_equal(line.get_ydata(), data2["density"].values)
    assert f.ax.get_ylabel() == "density [g / cm ** 3]"


def test_plot_linestyle(data2):
    f = osyris.plot(data2["position_x"], data2["density"], ls="dashed")
    (line,) = f.ax.lines
    assert line.get_linestyle() == "--"


def test_plot_color(data2):
    f = osyris.plot(data2["position_x"], data2["density"], color="red")
    (line,) = f.ax.lines
    assert line.get_color() == "red"


def test_plot_marker(data2):
    f = osyris.plot(data2["position_x"], data2["density"], marker="o", ls="None")
    (line,) = f.ax.lines
    assert line.get_marker() == "o"
    assert line.get_linestyle() == "None"


def test_plot_logx(data2):
    f = osyris.plot(data2["position_x"], data2["density"], logx=True)
    assert f.ax.get_xscale() == "log"
    assert f.ax.get_yscale() == "linear"


def test_plot_logy(data2):
    f = osyris.plot(data2["position_x"], data2["density"], logy=True)
    assert f.ax.get_xscale() == "linear"
    assert f.ax.get_yscale() == "log"


def test_plot_log_log(data2):
    f = osyris.plot(data2["position_x"], data2["density"], logx=True, logy=True)
    assert f.ax.get_xscale() == "log"
    assert f.ax.get_yscale() == "log"


def test_plot_two_lines(data2):
    f = osyris.plot(data2["position_x"], data2["density"], data2["density"] * 2)
    layer1, layer2 = f.layers
    assert np.array_equal(layer1["x"].values, data2["position_x"].values)
    assert np.array_equal(layer2["x"].values, data2["position_x"].values)
    assert np.array_equal(layer1["y"].values, data2["density"].values)
    assert np.array_equal(layer2["y"].values, data2["density"].values * 2)

    lines = f.ax.lines
    assert len(lines) == 2
    inds = np.argsort(data2["position_x"].values)
    assert np.array_equal(lines[0].get_xdata(), data2["position_x"].values[inds])
    assert np.array_equal(lines[0].get_ydata(), data2["density"].values[inds])
    assert np.array_equal(lines[1].get_xdata(), data2["position_x"].values[inds])
    assert np.array_equal(lines[1].get_ydata(), (data2["density"].values * 2)[inds])


def test_plot_two_lines_with_incompatible_units_raises(data2):
    with pytest.raises(ValueError):
        osyris.plot(data2["position_x"], data2["density"], data2["pressure"])


def test_plot_two_lines_with_different_x(data1, data2):
    f = osyris.plot(
        {"x": data1["position_x"], "y": data1["density"]},
        {"x": data2["position_x"], "y": data2["density"]},
    )
    layer1, layer2 = f.layers
    assert np.array_equal(layer1["x"].values, data1["position_x"].values)
    assert np.array_equal(layer2["x"].values, data2["position_x"].values)
    assert np.array_equal(layer1["y"].values, data1["density"].values)
    assert np.array_equal(layer2["y"].values, data2["density"].values)
    assert len(f.ax.lines) == 2


def test_plot_two_lines_with_different_x_common_linestyles(data1, data2):
    f = osyris.plot(
        {"x": data1["position_x"], "y": data1["density"]},
        {"x": data2["position_x"], "y": data2["density"]},
        marker="o",
        ls="None",
    )
    for line in f.ax.lines:
        assert line.get_linestyle() == "None"
        assert line.get_marker() == "o"


def test_plot_two_lines_with_different_x_different_linestyles(data1, data2):
    f = osyris.plot(
        {
            "x": data1["position_x"],
            "y": data1["density"],
            "ls": "dashed",
            "color": "blue",
        },
        {
            "x": data2["position_x"],
            "y": data2["density"],
            "ls": "dotted",
            "color": "red",
        },
    )
    line1, line2 = f.ax.lines
    assert line1.get_linestyle() == "--"
    assert line1.get_color() == "blue"
    assert line2.get_linestyle() == ":"
    assert line2.get_color() == "red"


def test_plot_using_subplots(data1):
    _, ax = plt.subplots(2, 2, figsize=(12, 9))
    osyris.plot(
        data1["position_x"], data1["density"], marker="o", ls="None", ax=ax[0, 0]
    )
    osyris.plot(data1["position_x"], data1["velocity_x"], color="green", ax=ax[0, 1])
    osyris.plot(data1["position_x"], data1["pressure"], color="red", ax=ax[1, 0])
    osyris.plot(
        data1["position_x"], data1["level"], color="black", ls="dotted", ax=ax[1, 1]
    )
