# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)
from osyris import Array, Datagroup
import pytest


def test_datagroup_creation():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[6., 7., 8., 9., 10.], unit='s')
    dg = Datagroup()
    dg['a'] = a
    dg['b'] = b
    assert 'a' in dg.keys()
    assert 'b' in dg.keys()
    assert dg['a'].name == 'a'
    assert dg['b'].name == 'b'
    assert len(dg) == 2


def test_datagroup_creation_from_dict():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[6., 7., 8., 9., 10.], unit='s')
    dg = Datagroup({'a': a, 'b': b})
    assert 'a' in dg.keys()
    assert 'b' in dg.keys()
    assert dg['a'].name == 'a'
    assert dg['b'].name == 'b'
    assert len(dg) == 2


def test_datagroup_insert_vector():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[[6., 7., 8.], [9., 10., 11.], [12., 13., 14.], [15., 16., 17.],
                      [18., 19., 20.]],
              unit='s')
    dg = Datagroup()
    dg['a'] = a
    dg['b'] = b
    assert 'a' in dg.keys()
    assert 'b' in dg.keys()
    assert dg['a'].name == 'a'
    assert dg['b'].name == 'b'
    assert len(dg) == 2


def test_datagroup_bad_length_insertion():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[6., 7., 8., 9.], unit='s')
    dg = Datagroup()
    dg['a'] = a
    with pytest.raises(RuntimeError):
        dg['b'] = b


def test_datagroup_delitem():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[6., 7., 8., 9., 10.], unit='s')
    dg = Datagroup()
    dg['a'] = a
    dg['b'] = b
    del dg['a']
    assert 'a' not in dg.keys()
    assert 'b' in dg.keys()
    assert len(dg) == 1


def test_datagroup_slice():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[6., 7., 8., 9., 10.], unit='s')
    dg = Datagroup()
    dg['a'] = a
    dg['b'] = b
    expected = Datagroup({
        'a': Array(values=[2.], unit='m'),
        'b': Array(values=[7.], unit='s')
    })
    sliced = dg[1]
    assert sliced['a'] == expected['a']
    assert sliced['b'] == expected['b']


def test_datagroup_slice_range():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[6., 7., 8., 9., 10.], unit='s')
    dg = Datagroup()
    dg['a'] = a
    dg['b'] = b
    expected = Datagroup({
        'a': Array(values=[2., 3., 4.], unit='m'),
        'b': Array(values=[7., 8., 9.], unit='s')
    })
    sliced = dg[1:4]
    assert all(sliced['a'] == expected['a'])
    assert all(sliced['b'] == expected['b'])


def test_datagroup_sortby_key():
    a = Array(values=[2., 3., 1., 5., 4.], unit='m')
    b = Array(values=[6., 7., 8., 9., 10.], unit='s')
    dg = Datagroup({'a': a, 'b': b})
    dg.sortby('a')
    expected = Datagroup({
        'a': Array(values=[1., 2., 3., 4., 5.], unit='m'),
        'b': Array(values=[8., 6., 7., 10., 9.], unit='s')
    })
    assert all(dg['a'] == expected['a'])
    assert all(dg['b'] == expected['b'])


def test_datagroup_sortby_indices():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[6., 7., 8., 9., 10.], unit='s')
    dg = Datagroup({'a': a, 'b': b})
    dg.sortby([4, 3, 1, 2, 0])
    expected = Datagroup({
        'a': Array(values=[5., 4., 2., 3., 1.], unit='m'),
        'b': Array(values=[10., 9., 7., 8., 6.], unit='s')
    })
    assert all(dg['a'] == expected['a'])
    assert all(dg['b'] == expected['b'])
