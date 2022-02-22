# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)
from osyris import Array, Datagroup, Dataset


def test_dataset_creation():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[6., 7., 8., 9., 10.], unit='s')
    dg1 = Datagroup({'a': a, 'b': b})
    dg2 = Datagroup({'a': a[:-1], 'b': b[:-1]})
    ds = Dataset()
    ds['one'] = dg1
    ds['two'] = dg2
    assert 'one' in ds
    assert 'two' in ds


def test_dataset_clear():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[6., 7., 8., 9., 10.], unit='s')
    dg1 = Datagroup({'a': a, 'b': b})
    dg2 = Datagroup({'a': a[:-1], 'b': b[:-1]})
    ds = Dataset()
    ds['one'] = dg1
    ds['two'] = dg2
    ds.clear()
    assert len(ds) == 0


def test_dataset_get():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[6., 7., 8., 9., 10.], unit='s')
    dg1 = Datagroup({'a': a, 'b': b})
    dg2 = Datagroup({'a': a[:-1], 'b': b[:-1]})
    ds = Dataset()
    ds['one'] = dg1
    ds['two'] = dg2
    assert ds.get('one', 42) == dg1
    assert ds.get('three', 42) == 42


def test_dataset_pop():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[6., 7., 8., 9., 10.], unit='s')
    dg1 = Datagroup({'a': a, 'b': b})
    dg2 = Datagroup({'a': a[:-1], 'b': b[:-1]})
    ds = Dataset()
    ds['one'] = dg1
    ds['two'] = dg2
    c = ds.pop('one')
    assert 'one' not in ds
    assert c == dg1


def test_dataset_update():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[6., 7., 8., 9., 10.], unit='s')
    c = Array(values=[11., 12., 13., 14., 15.], unit='K')
    dg1 = Datagroup({'a': a, 'b': b})
    dg2 = Datagroup({'a': a[:-1], 'b': b[:-1]})
    dg3 = Datagroup({'a': a, 'c': c})
    ds = Dataset()
    ds['one'] = dg1
    ds['two'] = dg2
    ds.update({'one': Datagroup({'a': a * 2.0}), 'three': dg3})
    assert 'three' in ds
    assert 'b' not in ds['one']
    assert 'c' in ds['three']
