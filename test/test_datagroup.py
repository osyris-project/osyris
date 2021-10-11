# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)
import osyris
import pytest


def test_datagroup_creation():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = osyris.Array(values=[6., 7., 8., 9., 10.], unit='s')
    dg = osyris.Datagroup()
    dg['a'] = a
    dg['b'] = b
    assert 'a' in dg.keys()
    assert 'b' in dg.keys()
    assert dg['a'].name == 'a'
    assert dg['b'].name == 'b'
    assert len(dg) == 2


def test_datagroup_creation_from_dict():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = osyris.Array(values=[6., 7., 8., 9., 10.], unit='s')
    dg = osyris.Datagroup({'a': a, 'b': b})
    assert 'a' in dg.keys()
    assert 'b' in dg.keys()
    assert dg['a'].name == 'a'
    assert dg['b'].name == 'b'
    assert len(dg) == 2


def test_datagroup_insert_vector():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = osyris.Array(values=[[6., 7., 8.], [9., 10., 11.], [12., 13., 14.],
                             [15., 16., 17.], [18., 19., 20.]],
                     unit='s')
    dg = osyris.Datagroup()
    dg['a'] = a
    dg['b'] = b
    assert 'a' in dg.keys()
    assert 'b' in dg.keys()
    assert dg['a'].name == 'a'
    assert dg['b'].name == 'b'
    assert len(dg) == 2


def test_datagroup_bad_length_insertion():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = osyris.Array(values=[6., 7., 8., 9.], unit='s')
    dg = osyris.Datagroup()
    dg['a'] = a
    with pytest.raises(RuntimeError):
        dg['b'] = b


def test_datagroup_delitem():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = osyris.Array(values=[6., 7., 8., 9., 10.], unit='s')
    dg = osyris.Datagroup()
    dg['a'] = a
    dg['b'] = b
    del dg['a']
    assert 'a' not in dg.keys()
    assert 'b' in dg.keys()
    assert len(dg) == 1


def test_datagroup_slice():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = osyris.Array(values=[6., 7., 8., 9., 10.], unit='s')
    dg = osyris.Datagroup()
    dg['a'] = a
    dg['b'] = b
    expected = osyris.Datagroup({
        'a': osyris.Array(values=[2.], unit='m'),
        'b': osyris.Array(values=[7.], unit='s')
    })
    sliced = dg[1]
    assert sliced['a'] == expected['a']
    assert sliced['b'] == expected['b']


def test_datagroup_slice_range():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = osyris.Array(values=[6., 7., 8., 9., 10.], unit='s')
    dg = osyris.Datagroup()
    dg['a'] = a
    dg['b'] = b
    expected = osyris.Datagroup({
        'a': osyris.Array(values=[2., 3., 4.], unit='m'),
        'b': osyris.Array(values=[7., 8., 9.], unit='s')
    })
    sliced = dg[1:4]
    assert all(sliced['a'] == expected['a'])
    assert all(sliced['b'] == expected['b'])
