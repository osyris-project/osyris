# SPDX-License-Identifier: BSD-3-Clause


import matplotlib.pyplot as plt
import numpy as np
import pytest

import osyris


def test_plot_line_profile(sod_data_path):
    data = osyris.RamsesDataset(2, path=sod_data_path).load()["mesh"]
    f = osyris.plot(data["position_x"], data["density"])
    layer = f.layers[0]
    assert np.array_equal(layer["x"].values, data["position_x"].values)
    assert np.array_equal(layer["y"].values, data["density"].values)
    (line,) = f.ax.lines
    inds = np.argsort(data["position_x"].values)
    assert np.array_equal(line.get_xdata(), data["position_x"].values[inds])
    assert np.array_equal(line.get_ydata(), data["density"].values[inds])
    assert f.ax.get_xlabel() == "position_x [cm]"
    assert f.ax.get_ylabel() == "density [g / cm ** 3]"


def test_plot_line_profile_no_x(sod_data_path):
    data = osyris.RamsesDataset(2, path=sod_data_path).load()["mesh"]
    f = osyris.plot(data["density"])
    layer = f.layers[0]
    arange = np.arange(len(data["density"].values))
    assert np.array_equal(layer["x"].values, arange)
    assert np.array_equal(layer["y"].values, data["density"].values)
    (line,) = f.ax.lines
    assert np.array_equal(line.get_xdata(), arange)
    assert np.array_equal(line.get_ydata(), data["density"].values)
    assert f.ax.get_ylabel() == "density [g / cm ** 3]"


def test_plot_linestyle(sod_data_path):
    data = osyris.RamsesDataset(2, path=sod_data_path).load()["mesh"]
    f = osyris.plot(data["position_x"], data["density"], ls="dashed")
    (line,) = f.ax.lines
    assert line.get_linestyle() == "--"


def test_plot_color(sod_data_path):
    data = osyris.RamsesDataset(2, path=sod_data_path).load()["mesh"]
    f = osyris.plot(data["position_x"], data["density"], color="red")
    (line,) = f.ax.lines
    assert line.get_color() == "red"


def test_plot_marker(sod_data_path):
    data = osyris.RamsesDataset(2, path=sod_data_path).load()["mesh"]
    f = osyris.plot(data["position_x"], data["density"], marker="o", ls="None")
    (line,) = f.ax.lines
    assert line.get_marker() == "o"
    assert line.get_linestyle() == "None"


def test_plot_logx(sod_data_path):
    data = osyris.RamsesDataset(2, path=sod_data_path).load()["mesh"]
    f = osyris.plot(data["position_x"], data["density"], logx=True)
    assert f.ax.get_xscale() == "log"
    assert f.ax.get_yscale() == "linear"


def test_plot_logy(sod_data_path):
    data = osyris.RamsesDataset(2, path=sod_data_path).load()["mesh"]
    f = osyris.plot(data["position_x"], data["density"], logy=True)
    assert f.ax.get_xscale() == "linear"
    assert f.ax.get_yscale() == "log"


def test_plot_log_log(sod_data_path):
    data = osyris.RamsesDataset(2, path=sod_data_path).load()["mesh"]
    f = osyris.plot(data["position_x"], data["density"], logx=True, logy=True)
    assert f.ax.get_xscale() == "log"
    assert f.ax.get_yscale() == "log"


def test_plot_two_lines(sod_data_path):
    data = osyris.RamsesDataset(2, path=sod_data_path).load()["mesh"]
    f = osyris.plot(data["position_x"], data["density"], data["density"] * 2)
    layer1, layer2 = f.layers
    assert np.array_equal(layer1["x"].values, data["position_x"].values)
    assert np.array_equal(layer2["x"].values, data["position_x"].values)
    assert np.array_equal(layer1["y"].values, data["density"].values)
    assert np.array_equal(layer2["y"].values, data["density"].values * 2)

    lines = f.ax.lines
    assert len(lines) == 2
    inds = np.argsort(data["position_x"].values)
    assert np.array_equal(lines[0].get_xdata(), data["position_x"].values[inds])
    assert np.array_equal(lines[0].get_ydata(), data["density"].values[inds])
    assert np.array_equal(lines[1].get_xdata(), data["position_x"].values[inds])
    assert np.array_equal(lines[1].get_ydata(), (data["density"].values * 2)[inds])


def test_plot_two_lines_with_incompatible_units_raises(sod_data_path):
    data = osyris.RamsesDataset(2, path=sod_data_path).load()["mesh"]
    with pytest.raises(ValueError):
        osyris.plot(
            data["position_x"],
            data["density"],
            data["pressure"],
        )


def test_plot_two_lines_with_different_x(sod_data_path):
    data1 = osyris.RamsesDataset(1, path=sod_data_path).load()["mesh"]
    data2 = osyris.RamsesDataset(2, path=sod_data_path).load()["mesh"]
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


def test_plot_two_lines_with_different_x_common_linestyles(sod_data_path):
    data1 = osyris.RamsesDataset(1, path=sod_data_path).load()["mesh"]
    data2 = osyris.RamsesDataset(2, path=sod_data_path).load()["mesh"]
    f = osyris.plot(
        {"x": data1["position_x"], "y": data1["density"]},
        {"x": data2["position_x"], "y": data2["density"]},
        marker="o",
        ls="None",
    )
    for line in f.ax.lines:
        assert line.get_linestyle() == "None"
        assert line.get_marker() == "o"


def test_plot_two_lines_with_different_x_different_linestyles(sod_data_path):
    data1 = osyris.RamsesDataset(1, path=sod_data_path).load()["mesh"]
    data2 = osyris.RamsesDataset(2, path=sod_data_path).load()["mesh"]
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


def test_plot_using_subplots(sod_data_path):
    data = osyris.RamsesDataset(1, path=sod_data_path).load()["mesh"]
    _, ax = plt.subplots(2, 2, figsize=(12, 9))
    osyris.plot(data["position_x"], data["density"], marker="o", ls="None", ax=ax[0, 0])
    osyris.plot(data["position_x"], data["velocity_x"], color="green", ax=ax[0, 1])
    osyris.plot(data["position_x"], data["pressure"], color="red", ax=ax[1, 0])
    osyris.plot(
        data["position_x"], data["level"], color="black", ls="dotted", ax=ax[1, 1]
    )
