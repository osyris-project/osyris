# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)
import osyris
import numpy as np
import pint
import pytest


def test_1d_array():
    a = np.arange(100.)
    array = osyris.Array(values=a, unit='m')
    assert array.unit == osyris.units('m')
    assert len(array) == len(a)
    assert array.shape == a.shape
    assert np.allclose(array.values, a)


def test_2d_array():
    a = np.random.random([100, 3])
    array = osyris.Array(values=a, unit='m')
    assert array.unit == osyris.units('m')
    assert np.allclose(array.x.values, a[:, 0])
    assert np.allclose(array.y.values, a[:, 1])
    assert np.allclose(array.z.values, a[:, 2])


def test_equal():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    c = osyris.Array(values=[100., 200., 300., 400., 500.], unit='cm')
    assert all(a == b)
    assert all(a == c)


def test_not_equal():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = osyris.Array(values=[1., 2., 3., 4., 5.], unit='cm')
    c = osyris.Array(values=[100., 200., 300., 400., 500.], unit='m')
    assert all(a != b)
    assert all(a != c)


def test_addition():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = osyris.Array(values=[6., 7., 8., 9., 10.], unit='m')
    expected = osyris.Array(values=[7., 9., 11., 13., 15.], unit='m')
    assert all(a + b == expected)


def test_addition_bad_units():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = osyris.Array(values=[6., 7., 8., 9., 10.], unit='s')
    with pytest.raises(pint.errors.DimensionalityError):
        _ = a + b
    with pytest.raises(TypeError):
        _ = a + 3.0


def test_addition_quantity():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = 3.5 * osyris.units('m')
    expected = osyris.Array(values=[4.5, 5.5, 6.5, 7.5, 8.5], unit='m')
    assert all(a + b == expected)


def test_addition_inplace():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = osyris.Array(values=[6., 7., 8., 9., 10.], unit='m')
    expected = osyris.Array(values=[7., 9., 11., 13., 15.], unit='m')
    a += b
    assert all(a == expected)


def test_addition_quantity_inplace():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = 3.5 * osyris.units('m')
    expected = osyris.Array(values=[4.5, 5.5, 6.5, 7.5, 8.5], unit='m')
    a += b
    assert all(a == expected)


def test_subtraction():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = osyris.Array(values=[6., 7., 8., 9., 10.], unit='m')
    expected = osyris.Array(values=[5., 5., 5., 5., 5.], unit='m')
    assert all(b - a == expected)


def test_subtraction_bad_units():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = osyris.Array(values=[6., 7., 8., 9., 10.], unit='s')
    with pytest.raises(pint.errors.DimensionalityError):
        _ = a - b
    with pytest.raises(TypeError):
        _ = a - 3.0


def test_subtraction_quantity():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = 3.5 * osyris.units('m')
    expected = osyris.Array(values=[-2.5, -1.5, -0.5, 0.5, 1.5], unit='m')
    assert all(a - b == expected)


def test_subtraction_inplace():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = osyris.Array(values=[6., 7., 8., 9., 10.], unit='m')
    expected = osyris.Array(values=[5., 5., 5., 5., 5.], unit='m')
    b -= a
    assert all(b == expected)


def test_subtraction_quantity_inplace():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = 3.5 * osyris.units('m')
    expected = osyris.Array(values=[-2.5, -1.5, -0.5, 0.5, 1.5], unit='m')
    a -= b
    assert all(a == expected)


def test_multiplication():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = osyris.Array(values=[6., 7., 8., 9., 10.], unit='m')
    expected = osyris.Array(values=[6., 14., 24., 36., 50.], unit='m*m')
    assert all(a * b == expected)


def test_multiplication_float():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = 3.0
    expected = osyris.Array(values=[3., 6., 9., 12., 15.], unit='m')
    assert all(a * b == expected)


def test_multiplication_quantity():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = 3.5 * osyris.units('s')
    expected = osyris.Array(values=[3.5, 7.0, 10.5, 14.0, 17.5], unit='m*s')
    assert all(a * b == expected)


def test_multiplication_inplace():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = osyris.Array(values=[6., 7., 8., 9., 10.], unit='m')
    expected = osyris.Array(values=[6., 14., 24., 36., 50.], unit='m*m')
    a *= b
    assert all(a == expected)


def test_multiplication_float_inplace():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = 3.0
    expected = osyris.Array(values=[3., 6., 9., 12., 15.], unit='m')
    a *= b
    assert all(a == expected)


def test_multiplication_quantity_inplace():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = 3.5 * osyris.units('s')
    expected = osyris.Array(values=[3.5, 7.0, 10.5, 14.0, 17.5], unit='m*s')
    a *= b
    assert all(a == expected)


def test_division():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='s')
    b = osyris.Array(values=[6., 7., 8., 9., 10.], unit='m')
    expected = osyris.Array(values=[600., 350., 800. / 3., 225.0, 200.], unit='cm/s')
    assert all(b / a == expected)


def test_division_float():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='s')
    b = 3.0
    expected = osyris.Array(values=[1. / 3., 2. / 3., 1., 4. / 3., 5. / 3.], unit='s')
    assert all(a / b == expected)


def test_division_quantity():
    a = osyris.Array(values=[0., 2., 4., 6., 200.], unit='s')
    b = 2.0 * osyris.units('s')
    expected = osyris.Array(values=[0., 1., 2., 3., 100.], unit='dimensionless')
    assert all(a / b == expected)


def test_division_inplace():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='s')
    b = osyris.Array(values=[6., 7., 8., 9., 10.], unit='m')
    expected = osyris.Array(values=[600., 350., 800. / 3., 225.0, 200.], unit='cm/s')
    b /= a
    assert all(b == expected)


def test_division_float_inplace():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='s')
    b = 3.0
    expected = osyris.Array(values=[1. / 3., 2. / 3., 1., 4. / 3., 5. / 3.], unit='s')
    a /= b
    assert all(a == expected)


def test_division_quantity_inplace():
    a = osyris.Array(values=[0., 2., 4., 6., 200.], unit='s')
    b = 2.0 * osyris.units('s')
    expected = osyris.Array(values=[0., 1., 2., 3., 100.], unit='dimensionless')
    a /= b
    assert all(a == expected)


def test_norm():
    a2d = osyris.Array(values=np.array([[1., 2.], [3., 4.], [5., 6.], [7., 8.]]),
                       unit='s')
    a3d = osyris.Array(values=np.array([[1., 2., 3.], [4., 5., 6.]]), unit='g')
    assert all(a2d.norm == osyris.Array(values=np.sqrt([5., 25., 61., 113.]), unit='s'))
    assert all(a3d.norm == osyris.Array(values=np.sqrt([14., 77.]), unit='g'))


def test_broadcast():
    a1d = osyris.Array(values=np.array([1., 2., 3., 4., 5.]), unit='s')
    a3d = osyris.Array(values=np.array([[1., 2., 3.], [4., 5., 6.], [7., 8., 9.],
                                        [10., 11., 12.], [13., 14., 15.]]),
                       unit='g')
    expected = osyris.Array(values=np.array([[1., 2., 3.], [8., 10., 12.],
                                             [21., 24., 27.], [40., 44., 48.],
                                             [65., 70., 75.]]),
                            unit='g*s')
    assert all(np.ravel(a1d * a3d == expected))


def test_power():
    a = osyris.Array(values=[1., 2., 4., 6., 200.], unit='s')
    expected = osyris.Array(values=[1., 8., 64., 216., 8.0e6], unit='s**3')
    assert all(a**3 == expected)


def test_less_than():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='s')
    b = osyris.Array(values=[6., 7., 1., 4., 10.], unit='s')
    expected = [True, True, False, False, True]
    np.array_equal(a < b, expected)


def test_less_than_bad_units():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='s')
    b = osyris.Array(values=[6., 7., 1., 4., 10.], unit='m')
    with pytest.raises(pint.errors.DimensionalityError):
        _ = a < b


def test_less_equal():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='s')
    b = osyris.Array(values=[6., 7., 1., 4., 10.], unit='s')
    expected = [True, True, False, True, True]
    np.array_equal(a <= b, expected)


def test_less_equal_bad_units():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='s')
    b = osyris.Array(values=[6., 7., 1., 4., 10.], unit='m')
    with pytest.raises(pint.errors.DimensionalityError):
        _ = a <= b


def test_greater_than():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='s')
    b = osyris.Array(values=[6., 7., 1., 4., 10.], unit='s')
    expected = [True, True, False, False, True]
    np.array_equal(b > a, expected)


def test_greater_than_bad_units():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='s')
    b = osyris.Array(values=[6., 7., 1., 4., 10.], unit='K')
    with pytest.raises(pint.errors.DimensionalityError):
        _ = b > a


def test_greater_equal():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='s')
    b = osyris.Array(values=[6., 7., 1., 4., 10.], unit='s')
    expected = [True, True, False, True, True]
    np.array_equal(b <= a, expected)


def test_greater_equal_bad_units():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='s')
    b = osyris.Array(values=[6., 7., 1., 4., 10.], unit='K')
    with pytest.raises(pint.errors.DimensionalityError):
        _ = b >= a


def test_to():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    b = osyris.Array(values=[1.0e-3, 2.0e-3, 3.0e-3, 4.0e-3, 5.0e-3], unit='km')
    assert all(a.to('km') == b)


def test_to_bad_units():
    a = osyris.Array(values=[1., 2., 3., 4., 5.], unit='m')
    with pytest.raises(pint.errors.DimensionalityError):
        _ = a.to('s')


def test_min():
    a = osyris.Array(values=[1., 2., 0.3, 4., 5.], unit='m')
    b = osyris.Array(values=[0.3], unit='m')
    assert all(a.min() == b)


def test_max():
    a = osyris.Array(values=[1., 2., 3., 45., 5.], unit='m')
    b = osyris.Array(values=[45.0], unit='m')
    assert all(a.max() == b)


def test_reshape():
    a = osyris.Array(values=[1., 2., 3., 4., 5., 6.], unit='m')
    expected = osyris.Array(values=[[1., 2., 3.], [4., 5., 6.]], unit='m')
    assert all(np.ravel(a.reshape(2, 3) == expected))
