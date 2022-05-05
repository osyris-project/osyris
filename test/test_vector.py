# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)
from common import (vectorclose, vectortrue, vectorequal, arraytrue, arrayequal,
                    arrayclose)
from osyris import Array, Vector, units
from copy import copy, deepcopy
import numpy as np
from pint.errors import DimensionalityError
import pytest


def test_constructor():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    assert len(v) == len(x)
    assert v.shape == x.shape


def test_constructor_bad_length():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14.], unit='m')
    with pytest.raises(ValueError):
        _ = Vector(x, y, z)


def test_constructor_bad_unit():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='cm')
    with pytest.raises(ValueError):
        _ = Vector(x, y, z)


def test_components():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y)
    assert arrayequal(v.x, x)
    assert arrayequal(v.y, y)
    assert v.z is None
    v.z = z
    assert arrayequal(v.z, z)
    assert v.x.name == "_x"
    assert v.y.name == "_y"


def test_equal():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v1 = Vector(x, y, z)
    v2 = Vector(x, y, z)
    assert vectortrue(v1 == v2)


def test_not_equal():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v1 = Vector(x, y, z)
    z = Array(values=[11., 12., 13., 19., 15.], unit='m')
    v2 = Vector(x, y, z)
    result = v1 != v2
    assert all(result.x.values == [False, False, False, False, False])
    assert all(result.z.values == [False, False, False, True, False])


def test_addition():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    expected = Vector(x=x + x, y=y + y, z=z + z)
    assert vectorequal(v + v, expected)


def test_addition_array():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    expected = Vector(x=x + x, y=y + x, z=z + x)
    assert vectorequal(v + x, expected)


def test_addition_quantity():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    q = 3.5 * units('m')
    expected = Vector(x=x + q, y=y + q, z=z + q)
    assert vectorequal(v + q, expected)


def test_addition_inplace():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    expected = Vector(x=x + x, y=y + y, z=z + z)
    v += v
    assert vectorequal(v, expected)


def test_addition_array_inplace():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    a = Array(values=[1.1, 2.2, 3.3, 4.4, 5.5], unit='m')
    expected = Vector(x=x + a, y=y + a, z=z + a)
    v += a
    assert vectorequal(v, expected)


def test_addition_quantity_inplace():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    q = 3.5 * units('m')
    expected = Vector(x=x + q, y=y + q, z=z + q)
    v += q
    assert vectorequal(v, expected)


def test_subtraction():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    expected = Vector(x=x - x, y=y - y, z=z - z)
    assert vectorequal(v - v, expected)


def test_subtraction_array():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    expected = Vector(x=x - x, y=y - x, z=z - x)
    assert vectorequal(v - x, expected)


def test_subtraction_quantity():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    q = 3.5 * units('m')
    expected = Vector(x=x - q, y=y - q, z=z - q)
    assert vectorequal(v - q, expected)


def test_subtraction_inplace():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    expected = Vector(x=x - x, y=y - y, z=z - z)
    v -= v
    assert vectorequal(v, expected)


def test_subtraction_array_inplace():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    a = Array(values=[1.1, 2.2, 3.3, 4.4, 5.5], unit='m')
    expected = Vector(x=x - a, y=y - a, z=z - a)
    v -= a
    assert vectorequal(v, expected)


def test_subtraction_quantity_inplace():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    q = 3.5 * units('m')
    expected = Vector(x=x - q, y=y - q, z=z - q)
    v -= q
    assert vectorequal(v, expected)


def test_multiplication():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    expected = Vector(x=x * x, y=y * y, z=z * z)
    assert vectorequal(v * v, expected)


def test_multiplication_array():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    expected = Vector(x=x * x, y=y * x, z=z * x)
    assert vectorequal(v * x, expected)
    assert vectorequal(x * v, expected)


def test_multiplication_float():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    f = 3.5
    expected = Vector(x=x * f, y=y * f, z=z * f)
    assert vectorequal(v * f, expected)
    assert vectorequal(f * v, expected)


def test_multiplication_ndarray():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    a = np.arange(5.)
    expected = Vector(x=x * a, y=y * a, z=z * a)
    assert vectorequal(v * a, expected)


def test_multiplication_quantity():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    q = 3.5 * units('m')
    expected = Vector(x=x * q, y=y * q, z=z * q)
    assert vectorequal(v * q, expected)


def test_multiplication_inplace():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    expected = Vector(x=x * x, y=y * y, z=z * z)
    v *= v
    assert vectorequal(v, expected)


def test_multiplication_array_inplace():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    a = Array(values=[1.1, 2.2, 3.3, 4.4, 5.5], unit='m')
    expected = Vector(x=x * a, y=y * a, z=z * a)
    v *= a
    assert vectorequal(v, expected)


def test_multiplication_float_inplace():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    f = 3.5
    expected = Vector(x=x * f, y=y * f, z=z * f)
    v *= f
    assert vectorequal(v, expected)


def test_multiplication_ndarray_inplace():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    a = np.arange(5.)
    expected = Vector(x=x * a, y=y * a, z=z * a)
    v *= a
    assert vectorequal(v, expected)


def test_multiplication_quantity_inplace():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    q = 3.5 * units('m')
    expected = Vector(x=x * q, y=y * q, z=z * q)
    v *= q
    assert vectorequal(v, expected)


def test_division():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v1 = Vector(x, y, z)
    w = Array(values=[16., 17., 18., 19., 20.], unit='m')
    v2 = Vector(w, y, w)
    expected = Vector(x=x / w, y=y / y, z=z / w)
    assert vectorequal(v1 / v2, expected)


def test_division_array():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    expected = Vector(x=x / x, y=y / x, z=z / x)
    assert vectorequal(v / x, expected)


def test_division_float():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    f = 3.5
    expected = Vector(x=x / f, y=y / f, z=z / f)
    assert vectorequal(v / f, expected)


def test_division_ndarray():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    a = np.arange(3., 8.)
    expected = Vector(x=x / a, y=y / a, z=z / a)
    assert vectorequal(v / a, expected)


def test_division_quantity():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    q = 3.5 * units('m')
    expected = Vector(x=x / q, y=y / q, z=z / q)
    assert vectorequal(v / q, expected)


def test_division_inplace():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v1 = Vector(x, y, z)
    w = Array(values=[16., 17., 18., 19., 20.], unit='m')
    v2 = Vector(w, y, w)
    expected = Vector(x=x / w, y=y / y, z=z / w)
    v1 /= v2
    assert vectorequal(v1, expected)


def test_division_array_inplace():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    a = Array(values=[1.1, 2.2, 3.3, 4.4, 5.5], unit='m')
    expected = Vector(x=x / a, y=y / a, z=z / a)
    v /= a
    assert vectorequal(v, expected)


def test_division_float_inplace():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    f = 3.5
    expected = Vector(x=x / f, y=y / f, z=z / f)
    v /= f
    assert vectorequal(v, expected)


def test_division_ndarray_inplace():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    a = np.arange(3., 8.)
    expected = Vector(x=x / a, y=y / a, z=z / a)
    v /= a
    assert vectorequal(v, expected)


def test_division_quantity_inplace():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    q = 3.5 * units('m')
    expected = Vector(x=x / q, y=y / q, z=z / q)
    v /= q
    assert vectorequal(v, expected)


def test_norm():
    x = Array(values=[1., 2., 3., 4., 5.], unit='s')
    y = Array(values=[6., 7., 8., 9., 10.], unit='s')
    z = Array(values=[11., 12., 13., 14., 15.], unit='s')
    v = Vector(x, y)
    assert arrayclose(v.norm, Array(values=np.sqrt([37., 53., 73., 97., 125.]),
                                    unit='s'))
    v.z = z
    assert arrayclose(v.norm,
                      Array(values=np.sqrt([158., 197., 242., 293., 350.]), unit='s'))


def test_power():
    x = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v = Vector(x, y, z)
    expected = Vector(x=x**3, y=y**3, z=z**3)
    assert vectorequal(v**3, expected)


def test_less_than():
    x1 = Array(values=[1., 2., 3., 4., 5.], unit='m')
    y1 = Array(values=[6., 7., 8., 9., 10.], unit='m')
    z1 = Array(values=[11., 12., 13., 14., 15.], unit='m')
    v1 = Vector(x1, y1, z1)
    x2 = Array(values=[-1., 2., 13., 4., 2.], unit='m')
    y2 = Array(values=[6., -7., 8., 3.3, 10.1], unit='m')
    z2 = Array(values=[11., 22., 7., 7., 15.], unit='m')
    v2 = Vector(x2, y2, z2)
    exp_x = [False, False, True, False, False]
    exp_y = [False, False, False, False, True]
    exp_z = [False, True, False, False, False]
    result = v1 < v2
    assert all(result.x.values == exp_x)
    assert all(result.y.values == exp_y)
    assert all(result.z.values == exp_z)


def test_less_than_conversion():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[600., 700., 100., 400., 1000.], unit='cm')
    expected = [True, True, False, False, True]
    assert all((a < b).values == expected)


def test_less_than_bad_units():
    a = Array(values=[1., 2., 3., 4., 5.], unit='s')
    b = Array(values=[6., 7., 1., 4., 10.], unit='m')
    with pytest.raises(DimensionalityError):
        _ = a < b


def test_less_equal():
    a = Array(values=[1., 2., 3., 4., 5.], unit='s')
    b = Array(values=[6., 7., 1., 4., 10.], unit='s')
    expected = [True, True, False, True, True]
    assert all((a <= b).values == expected)


def test_less_equal_bad_units():
    a = Array(values=[1., 2., 3., 4., 5.], unit='s')
    b = Array(values=[6., 7., 1., 4., 10.], unit='m')
    with pytest.raises(DimensionalityError):
        _ = a <= b


def test_greater_than():
    a = Array(values=[1., 2., 3., 4., 5.], unit='s')
    b = Array(values=[6., 7., 1., 4., 10.], unit='s')
    expected = [True, True, False, False, True]
    assert all((b > a).values == expected)


def test_greater_than_bad_units():
    a = Array(values=[1., 2., 3., 4., 5.], unit='s')
    b = Array(values=[6., 7., 1., 4., 10.], unit='K')
    with pytest.raises(DimensionalityError):
        _ = b > a


def test_greater_equal():
    a = Array(values=[1., 2., 3., 4., 5.], unit='s')
    b = Array(values=[6., 7., 1., 4., 10.], unit='s')
    expected = [True, True, False, True, True]
    assert all((b >= a).values == expected)


def test_greater_equal_bad_units():
    a = Array(values=[1., 2., 3., 4., 5.], unit='s')
    b = Array(values=[6., 7., 1., 4., 10.], unit='K')
    with pytest.raises(DimensionalityError):
        _ = b >= a


def test_to():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = Array(values=[1.0e-3, 2.0e-3, 3.0e-3, 4.0e-3, 5.0e-3], unit='km')
    assert allclose(a.to('km'), b)
    assert a.unit.units == units('m')


def test_to_bad_units():
    a = Array(values=[1., 2., 3., 4., 5.], unit='m')
    with pytest.raises(DimensionalityError):
        _ = a.to('s')


def test_min():
    a = Array(values=[1., -2., 3., 0.4, 0.5, 0.6], unit='m')
    assert (a.min() == Array(values=-2., unit='m')).values
    b = Array(values=np.array([1., -2., 3., 0.4, 0.5, 0.6]).reshape(2, 3), unit='m')
    assert (b.min() == Array(values=-2., unit='m')).values


def test_max():
    a = Array(values=[1., 2., 3., -15., 5., 6.], unit='m')
    assert (a.max() == Array(values=6.0, unit='m')).values
    b = Array(values=np.array([1., 2., 3., -15., 5., 6.]).reshape(2, 3), unit='m')
    assert (b.max() == Array(values=6.0, unit='m')).values


def test_reshape():
    a = Array(values=[1., 2., 3., 4., 5., 6.], unit='m')
    expected = Array(values=[[1., 2., 3.], [4., 5., 6.]], unit='m')
    assert alltrue(np.ravel(a.reshape(2, 3) == expected))


def test_slicing():
    a = Array(values=[11., 12., 13., 14., 15.], unit='m')
    assert a[2] == Array(values=[13.], unit='m')
    assert alltrue(a[:4] == Array(values=[11., 12., 13., 14.], unit='m'))
    assert alltrue(a[2:4] == Array(values=[13., 14.], unit='m'))


def test_slicing_vector():
    a = Array(values=np.arange(12.).reshape(4, 3), unit='m')
    assert alltrue(np.ravel(a[2:3] == Array(values=[[6., 7., 8.]], unit='m')))
    assert a[2:3].shape == (1, 3)
    assert alltrue(
        np.ravel(a[:2] == Array(values=[[0., 1., 2.], [3., 4., 5.]], unit='m')))


def test_copy():
    a = Array(values=[11., 12., 13., 14., 15.], unit='m')
    b = a.copy()
    a *= 10.
    assert alltrue(b == Array(values=[11., 12., 13., 14., 15.], unit='m'))


def test_copy_overload():
    a = Array(values=[11., 12., 13., 14., 15.], unit='m')
    b = copy(a)
    a *= 10.
    assert alltrue(b == Array(values=[11., 12., 13., 14., 15.], unit='m'))


def test_deepcopy():
    a = Array(values=[11., 12., 13., 14., 15.], unit='m')
    b = deepcopy(a)
    a *= 10.
    assert alltrue(b == Array(values=[11., 12., 13., 14., 15.], unit='m'))


def test_numpy_unary():
    values = [1., 2., 3., 4., 5.]
    a = Array(values=values, unit='m')
    expected = np.log10(values)
    result = np.log10(a)
    assert np.allclose(result.values, expected)
    assert result.unit == units('m')


def test_numpy_sqrt():
    values = [1., 2., 3., 4., 5.]
    a = Array(values=values, unit='m*m')
    expected = np.sqrt(values)
    result = np.sqrt(a)
    assert np.allclose(result.values, expected)
    assert result.unit == units('m')


def test_numpy_binary():
    a_buf = [1., 2., 3., 4., 5.]
    b_buf = [6., 7., 8., 9., 10.]
    a = Array(values=a_buf, unit='m')
    b = Array(values=b_buf, unit='m')
    expected = np.dot(a_buf, b_buf)
    result = np.dot(a, b)
    assert result.values == expected
    assert result.unit == units('m')
