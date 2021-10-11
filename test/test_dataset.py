# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)
import osyris


def test_dataset_creation():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = osyris.Array(values=[6., 7., 8., 9., 10.], unit='s')
    dg1 = osyris.Datagroup()
    dg1['a'] = a
    dg1['b'] = b
    dg2 = osyris.Datagroup()
    dg2['a'] = a[:-1]
    dg2['b'] = b[:-1]
    ds = osyris.Dataset()
    ds['one'] = dg1
    ds['two'] = dg2
    assert 'one' in ds
    assert 'two' in ds
