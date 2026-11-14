# SPDX-License-Identifier: BSD-3-Clause

import numpy as np

import osyris


def test_load_ramses_data(starformation_data_path):
    data = osyris.RamsesDataset(8, path=starformation_data_path).load()
    assert "mesh" in data.keys()


def test_select_cells_based_on_value(starformation_data_path):
    data = osyris.RamsesDataset(8, path=starformation_data_path).load(
        select={"mesh": {"density": lambda d: d > 1.0e-15 * osyris.units("g/cm**3")}}
    )
    assert (data["mesh"]["density"] > 1.0e-15 * osyris.units("g/cm**3")).values.all()


def test_multiple_selection_criteria(starformation_data_path):
    data = osyris.RamsesDataset(8, path=starformation_data_path).load(
        select={
            "mesh": {
                "density": lambda d: d > 1.0e-16 * osyris.units("g/cm**3"),
                "position_x": lambda x: x > 1500.0 * osyris.units("au"),
            },
        }
    )
    assert (data["mesh"]["density"] > 1.0e-16 * osyris.units("g/cm**3")).values.all()
    assert (data["mesh"]["position"].x > 1500.0 * osyris.units("au")).values.all()


def test_select_only_certain_groups(starformation_data_path):
    data = osyris.RamsesDataset(8, path=starformation_data_path).load(["mesh", "sink"])
    assert len(data) == 2
    assert "mesh" in data.keys()
    assert "sink" in data.keys()


def test_select_only_certain_variables_inside_a_group(starformation_data_path):
    data = osyris.RamsesDataset(8, path=starformation_data_path).load(
        select={"mesh": ["density", "position_x", "velocity_y"]}
    )
    assert "mesh" in data.keys()
    assert "density" in data["mesh"].keys()
    assert "position_x" in data["mesh"].keys()
    assert "velocity_y" in data["mesh"].keys()


def test_select_with_maximum_amr_level(starformation_data_path):
    data = osyris.RamsesDataset(8, path=starformation_data_path).load(
        select={"mesh": {"level": lambda lev: lev < 7}}
    )
    assert (data["mesh"]["level"].values < 7).all()


def test_select_range_of_amr_levels(starformation_data_path):
    data = osyris.RamsesDataset(8, path=starformation_data_path).load(
        select={"mesh": {"level": lambda lev: np.logical_and(lev > 5, lev < 9)}}
    )
    assert (data["mesh"]["level"].values > 5).all()
    assert (data["mesh"]["level"].values < 9).all()


def test_select_cpu_list(starformation_data_path):
    data = osyris.RamsesDataset(8, path=starformation_data_path).load(
        cpu_list=[1, 2, 10, 4, 5]
    )
    assert np.unique(data["mesh"]["cpu"].values).tolist() == [1, 2, 4, 5, 10]
