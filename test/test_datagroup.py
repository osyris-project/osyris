# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)
from common import arrayequal
from osyris import Array, Datagroup
from copy import copy, deepcopy
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


def test_datagroup_bad_length_insertion():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[6., 7., 8., 9.], unit='s')
    dg = Datagroup()
    dg['a'] = a
    with pytest.raises(ValueError):
        dg['b'] = b


def test_datagroup_delitem():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[6., 7., 8., 9., 10.], unit='s')
    dg = Datagroup({'a': a, 'b': b})
    del dg['a']
    assert 'a' not in dg.keys()
    assert 'b' in dg.keys()
    assert len(dg) == 1


def test_datagroup_equal_operator():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[6., 7., 8., 9., 10.], unit='s')
    assert Datagroup({'a': a, 'b': b}) == Datagroup({'a': a, 'b': b})


def test_datagroup_slice():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[6., 7., 8., 9., 10.], unit='s')
    dg = Datagroup({'a': a, 'b': b})
    expected = Datagroup({
        'a': Array(values=[2.], unit='m'),
        'b': Array(values=[7.], unit='s')
    })
    assert dg[1] == expected


def test_datagroup_slice_range():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[6., 7., 8., 9., 10.], unit='s')
    dg = Datagroup({'a': a, 'b': b})
    expected = Datagroup({
        'a': Array(values=[2., 3., 4.], unit='m'),
        'b': Array(values=[7., 8., 9.], unit='s')
    })
    assert dg[1:4] == expected


def test_datagroup_sortby_key():
    a = Array(values=[2., 3., 1., 5., 4.], unit='m')
    b = Array(values=[6., 7., 8., 9., 10.], unit='s')
    dg = Datagroup({'a': a, 'b': b})
    dg.sortby('a')
    expected = Datagroup({
        'a': Array(values=[1., 2., 3., 4., 5.], unit='m'),
        'b': Array(values=[8., 6., 7., 10., 9.], unit='s')
    })
    assert arrayequal(dg['a'], expected['a'])
    assert arrayequal(dg['b'], expected['b'])


def test_datagroup_sortby_indices():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[6., 7., 8., 9., 10.], unit='s')
    dg = Datagroup({'a': a, 'b': b})
    dg.sortby([4, 3, 1, 2, 0])
    expected = Datagroup({
        'a': Array(values=[5., 4., 2., 3., 1.], unit='m'),
        'b': Array(values=[10., 9., 7., 8., 6.], unit='s')
    })
    assert arrayequal(dg['a'], expected['a'])
    assert arrayequal(dg['b'], expected['b'])


def test_datagroup_clear():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[6., 7., 8., 9., 10.], unit='s')
    dg = Datagroup({'a': a, 'b': b})
    dg.clear()
    assert len(dg) == 0


def test_datagroup_get():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[6., 7., 8., 9., 10.], unit='s')
    c = Array(values=[11., 12., 13., 14., 15.], unit='K')
    dg = Datagroup({'a': a, 'b': b})
    assert arrayequal(dg.get('a', 42), a)
    assert dg.get('c', 42) == 42
    assert arrayequal(dg.get('c', c), c)


def test_datagroup_pop():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[6., 7., 8., 9., 10.], unit='s')
    dg = Datagroup({'a': a, 'b': b})
    c = dg.pop('a')
    assert 'a' not in dg
    assert arrayequal(c, a)


def test_datagroup_update():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[6., 7., 8., 9., 10.], unit='s')
    c = Array(values=[11., 12., 13., 14., 15.], unit='K')
    dg = Datagroup({'a': a, 'b': b})
    dg.update({'b': 2.0 * b, 'c': c})
    assert arrayequal(dg['b'], Array(values=[12., 14., 16., 18., 20.], unit='s'))
    assert arrayequal(dg['c'], c)


def test_copy():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[6., 7., 8., 9., 10.], unit='s')
    dg1 = Datagroup({'a': a, 'b': b})
    dg2 = dg1.copy()
    del dg1['b']
    assert 'b' in dg2
    dg1['a'] *= 10.
    assert arrayequal(dg2['a'], Array(values=[10., 20., 30., 40., 50.], unit='m'))


def test_copy_overload():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[6., 7., 8., 9., 10.], unit='s')
    dg1 = Datagroup({'a': a, 'b': b})
    dg2 = copy(dg1)
    del dg1['b']
    assert 'b' in dg2
    dg1['a'] *= 10.
    assert arrayequal(dg2['a'], Array(values=[10., 20., 30., 40., 50.], unit='m'))


def test_deepcopy():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[6., 7., 8., 9., 10.], unit='s')
    dg1 = Datagroup({'a': a, 'b': b})
    dg2 = deepcopy(dg1)
    del dg1['b']
    assert 'b' in dg2
    dg1['a'] *= 10.
    assert arrayequal(dg2['a'], Array(values=[1., 2., 3., 4., 5.], unit='m'))


def test_datagroup_shape():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    dg = Datagroup()
    dg['a'] = a
    assert dg.shape == (5, )
    a = Array(values=[[1., 2., 3.], [4., 5., 6.]], unit='m')
    dg = Datagroup()
    dg['a'] = a
    assert dg.shape == (2, 3)
