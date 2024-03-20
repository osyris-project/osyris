# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)
from common import arrayequal
from osyris import Array, Datagroup, Dataset
from copy import copy, deepcopy


def test_dataset_creation():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="s")
    dg1 = Datagroup({"a": a, "b": b})
    dg2 = Datagroup({"a": a[:-1], "b": b[:-1]})
    ds = Dataset()
    ds["one"] = dg1
    ds["two"] = dg2
    assert "one" in ds
    assert "two" in ds


def test_dataset_clear():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="s")
    dg1 = Datagroup({"a": a, "b": b})
    dg2 = Datagroup({"a": a[:-1], "b": b[:-1]})
    ds = Dataset()
    ds["one"] = dg1
    ds["two"] = dg2
    ds.clear()
    assert len(ds) == 0


def test_dataset_get():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="s")
    dg1 = Datagroup({"a": a, "b": b})
    dg2 = Datagroup({"a": a[:-1], "b": b[:-1]})
    ds = Dataset()
    ds["one"] = dg1
    ds["two"] = dg2
    assert ds.get("one", 42) == dg1
    assert ds.get("three", 42) == 42


def test_dataset_pop():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="s")
    dg1 = Datagroup({"a": a, "b": b})
    dg2 = Datagroup({"a": a[:-1], "b": b[:-1]})
    ds = Dataset()
    ds["one"] = dg1
    ds["two"] = dg2
    c = ds.pop("one")
    assert "one" not in ds
    assert c == dg1


def test_dataset_update():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="s")
    c = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="K")
    dg1 = Datagroup({"a": a, "b": b})
    dg2 = Datagroup({"a": a[:-1], "b": b[:-1]})
    dg3 = Datagroup({"a": a, "c": c})
    ds = Dataset()
    ds["one"] = dg1
    ds["two"] = dg2
    ds.update({"one": Datagroup({"a": a * 2.0}), "three": dg3})
    assert "three" in ds
    assert "b" not in ds["one"]
    assert "c" in ds["three"]


def test_copy():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="s")
    dg1 = Datagroup({"a": a, "b": b})
    dg2 = Datagroup({"a": a[:-1], "b": b[:-1]})
    ds1 = Dataset()
    ds1["one"] = dg1
    ds1["two"] = dg2
    ds1.meta = {"metadata1": [1, 2, 3], "metadata2": "4,5,6"}
    ds2 = ds1.copy()
    del ds1["two"]
    assert "two" in ds2
    del dg1["b"]
    assert "b" not in ds2["one"]
    dg1["a"] *= 10.0
    assert arrayequal(
        ds2["one"]["a"], Array(values=[10.0, 20.0, 30.0, 40.0, 50.0], unit="m")
    )
    del ds1.meta["metadata2"]
    assert "metadata2" in ds2.meta
    ds1.meta["metadata1"][0] = 0
    assert ds2.meta["metadata1"][0] == 0


def test_copy_overload():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="s")
    dg1 = Datagroup({"a": a, "b": b})
    dg2 = Datagroup({"a": a[:-1], "b": b[:-1]})
    ds1 = Dataset()
    ds1["one"] = dg1
    ds1["two"] = dg2
    ds1.meta = {"metadata1": [1, 2, 3], "metadata2": "4,5,6"}
    ds2 = copy(ds1)
    del ds1["two"]
    assert "two" in ds2
    del dg1["b"]
    assert "b" not in ds2["one"]
    dg1["a"] *= 10.0
    assert arrayequal(
        ds2["one"]["a"], Array(values=[10.0, 20.0, 30.0, 40.0, 50.0], unit="m")
    )
    del ds1.meta["metadata2"]
    assert "metadata2" in ds2.meta
    ds1.meta["metadata1"][0] = 0
    assert ds2.meta["metadata1"][0] == 0


def test_deepcopy():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="s")
    dg1 = Datagroup({"a": a, "b": b})
    dg2 = Datagroup({"a": a[:-1], "b": b[:-1]})
    ds1 = Dataset()
    ds1["one"] = dg1
    ds1["two"] = dg2
    ds1.meta = {"metadata1": [1, 2, 3], "metadata2": "4,5,6"}
    ds2 = deepcopy(ds1)
    del ds1["two"]
    assert "two" in ds2
    del dg1["b"]
    assert "b" in ds2["one"]
    dg1["a"] *= 10.0
    assert arrayequal(
        ds2["one"]["a"], Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    )
    del ds1.meta["metadata2"]
    assert "metadata2" in ds2.meta
    ds1.meta["metadata1"][0] = 0
    assert ds2.meta["metadata1"][0] == 1
