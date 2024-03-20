# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)
from common import arrayclose, arraytrue, arrayequal
from osyris import Array, units
from copy import copy, deepcopy
import numpy as np
from pint.errors import DimensionalityError
import pytest


def test_constructor_ndarray():
    a = np.arange(100.0)
    array = Array(values=a, unit="m")
    assert array.unit == units("m")
    assert len(array) == len(a)
    assert array.shape == a.shape
    assert np.array_equal(array.values, a)


def test_constructor_list():
    alist = [1.0, 2.0, 3.0, 4.0, 5.0]
    array = Array(values=alist, unit="s")
    assert array.unit == units("s")
    assert np.array_equal(array.values, alist)


def test_constructor_int():
    num = 15
    array = Array(values=num, unit="m")
    assert array.unit == units("m")
    assert np.array_equal(array.values, np.array(num))


def test_constructor_float():
    num = 154.77
    array = Array(values=num, unit="m")
    assert array.unit == units("m")
    assert np.array_equal(array.values, np.array(num))


def test_constructor_quantity():
    q = 6.7 * units("K")
    array = Array(values=q)
    assert array.unit == units("K")
    assert np.array_equal(array.values, np.array(q.magnitude))


def test_bad_constructor_quantity_with_unit():
    q = 6.7 * units("K")
    with pytest.raises(ValueError):
        _ = Array(values=q, unit="s")


def test_constructor_masked_array():
    a = np.arange(5.0)
    b = np.ma.masked_where(a > 2, a)
    array = Array(values=b, unit="m")
    assert array.unit == units("m")
    assert len(array) == len(b)
    assert array.shape == b.shape
    assert np.array_equal(array.values, b)
    assert np.array_equal(array.values.mask, [False, False, False, True, True])


def test_constructor_2d():
    a = np.arange(12.0).reshape(4, 3)
    array = Array(values=a, unit="m")
    assert array.unit == units("m")
    assert len(array) == len(a)
    assert array.shape == a.shape
    assert np.array_equal(array.values, a)


def test_addition():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    expected = Array(values=[7.0, 9.0, 11.0, 13.0, 15.0], unit="m")
    assert arrayclose(a + b, expected)


def test_addition_conversion():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="cm")
    expected = Array(values=[1.06, 2.07, 3.08, 4.09, 5.1], unit="m")
    assert arrayclose(a + b, expected)


def test_addition_bad_units():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="s")
    with pytest.raises(DimensionalityError):
        _ = a + b
    with pytest.raises(TypeError):
        _ = a + 3.0


def test_addition_quantity():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = 3.5 * units("m")
    expected = Array(values=[4.5, 5.5, 6.5, 7.5, 8.5], unit="m")
    assert arrayclose(a + b, expected)


def test_addition_inplace():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    expected = Array(values=[7.0, 9.0, 11.0, 13.0, 15.0], unit="m")
    a += b
    assert arrayclose(a, expected)


def test_addition_quantity_inplace():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = 3.5 * units("m")
    expected = Array(values=[4.5, 5.5, 6.5, 7.5, 8.5], unit="m")
    a += b
    assert arrayclose(a, expected)


def test_subtraction():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    expected = Array(values=[5.0, 5.0, 5.0, 5.0, 5.0], unit="m")
    assert arrayclose(b - a, expected)


def test_subtraction_bad_units():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="s")
    with pytest.raises(DimensionalityError):
        _ = a - b
    with pytest.raises(TypeError):
        _ = a - 3.0


def test_subtraction_quantity():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = 3.5 * units("m")
    expected = Array(values=[-2.5, -1.5, -0.5, 0.5, 1.5], unit="m")
    assert arrayclose(a - b, expected)


def test_subtraction_inplace():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    expected = Array(values=[5.0, 5.0, 5.0, 5.0, 5.0], unit="m")
    b -= a
    assert arrayclose(b, expected)


def test_subtraction_quantity_inplace():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = 3.5 * units("m")
    expected = Array(values=[-2.5, -1.5, -0.5, 0.5, 1.5], unit="m")
    a -= b
    assert arrayclose(a, expected)


def test_multiplication():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    expected = Array(values=[6.0, 14.0, 24.0, 36.0, 50.0], unit="m*m")
    assert arrayclose(a * b, expected)


def test_multiplication_conversion():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="cm")
    expected = Array(values=[0.06, 0.14, 0.24, 0.36, 0.5], unit="m*m")
    assert arrayclose(a * b, expected)


def test_multiplication_float():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = 3.0
    expected = Array(values=[3.0, 6.0, 9.0, 12.0, 15.0], unit="m")
    assert arrayclose(a * b, expected)
    assert arrayclose(b * a, expected)


def test_multiplication_ndarray():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = np.arange(5.0)
    expected = Array(values=[0.0, 2.0, 6.0, 12.0, 20.0], unit="m")
    assert arrayclose(a * b, expected)
    assert arrayclose(b * a, expected)


def test_multiplication_quantity():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = 3.5 * units("s")
    expected = Array(values=[3.5, 7.0, 10.5, 14.0, 17.5], unit="m*s")
    assert arrayclose(a * b, expected)


def test_multiplication_inplace():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    expected = Array(values=[6.0, 14.0, 24.0, 36.0, 50.0], unit="m*m")
    a *= b
    assert arrayclose(a, expected)


def test_multiplication_float_inplace():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = 3.0
    expected = Array(values=[3.0, 6.0, 9.0, 12.0, 15.0], unit="m")
    a *= b
    assert arrayclose(a, expected)


def test_multiplication_ndarray_inplace():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = np.arange(5.0)
    expected = Array(values=[0.0, 2.0, 6.0, 12.0, 20.0], unit="m")
    a *= b
    assert arrayclose(a, expected)


def test_multiplication_quantity_inplace():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = 3.5 * units("s")
    expected = Array(values=[3.5, 7.0, 10.5, 14.0, 17.5], unit="m*s")
    a *= b
    assert arrayclose(a, expected)


def test_division():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="s")
    b = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    expected = Array(values=[6.0, 3.5, 8.0 / 3.0, 2.25, 2.0], unit="m/s")
    assert arrayclose(b / a, expected)


def test_division_float():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="s")
    b = 3.0
    expected = Array(values=[1.0 / 3.0, 2.0 / 3.0, 1.0, 4.0 / 3.0, 5.0 / 3.0], unit="s")
    assert arrayclose(a / b, expected)
    expected = Array(values=[3.0, 3.0 / 2.0, 1.0, 3.0 / 4.0, 3.0 / 5.0], unit="1/s")
    assert arrayclose(b / a, expected)


def test_division_ndarray():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="s")
    b = np.arange(5.0, 10.0)
    expected = Array(
        values=[1.0 / 5.0, 2.0 / 6.0, 3.0 / 7.0, 4.0 / 8.0, 5.0 / 9.0], unit="s"
    )
    assert arrayclose(a / b, expected)
    # expected = Array(values=[3., 3. / 2., 1., 3. / 4., 3. / 5.], unit='1/s')
    # assert arrayclose(b / a, expected)


def test_division_quantity():
    a = Array(values=[0.0, 2.0, 4.0, 6.0, 200.0], unit="s")
    b = 2.0 * units("s")
    expected = Array(values=[0.0, 1.0, 2.0, 3.0, 100.0], unit="dimensionless")
    assert arrayclose(a / b, expected)


def test_division_inplace():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="s")
    b = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    expected = Array(values=[6.0, 3.5, 8.0 / 3.0, 2.25, 2.0], unit="m/s")
    b /= a
    assert arrayclose(b, expected)


def test_division_float_inplace():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="s")
    b = 3.0
    expected = Array(values=[1.0 / 3.0, 2.0 / 3.0, 1.0, 4.0 / 3.0, 5.0 / 3.0], unit="s")
    a /= b
    assert arrayclose(a, expected)


def test_division_ndarray_inplace():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="s")
    b = np.arange(5.0, 10.0)
    expected = Array(
        values=[1.0 / 5.0, 2.0 / 6.0, 3.0 / 7.0, 4.0 / 8.0, 5.0 / 9.0], unit="s"
    )
    a /= b
    assert arrayclose(a, expected)
    # expected = Array(values=[3., 3. / 2., 1., 3. / 4., 3. / 5.], unit='1/s')
    # assert arrayclose(b / a, expected)


def test_division_quantity_inplace():
    a = Array(values=[0.0, 2.0, 4.0, 6.0, 200.0], unit="s")
    b = 2.0 * units("s")
    expected = Array(values=[0.0, 1.0, 2.0, 3.0, 100.0], unit="dimensionless")
    a /= b
    assert arrayclose(a, expected)


def test_power():
    a = Array(values=[1.0, 2.0, 4.0, 6.0, 200.0], unit="s")
    expected = Array(values=[1.0, 8.0, 64.0, 216.0, 8.0e6], unit="s**3")
    assert arrayclose(a**3, expected)


def test_negative():
    a = Array(values=[1.0, 2.0, 4.0, 6.0, 200.0], unit="s")
    expected = Array(values=[-1.0, -2.0, -4.0, -6.0, -200.0], unit="s")
    assert arrayequal(-a, expected)


def test_equal():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = Array(values=[11.0, 2.0, 3.0, 4.1, 5.0], unit="m")
    expected = [False, True, True, False, True]
    assert all((a == b).values == expected)


def test_equal_conversion():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = Array(values=[1100.0, 200.0, 300.0, 410.0, 500.0], unit="cm")
    expected = [False, True, True, False, True]
    assert all((a == b).values == expected)


def test_equal_bad_units():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = Array(values=[11.0, 2.0, 3.0, 4.1, 5.0], unit="s")
    with pytest.raises(DimensionalityError):
        _ = a == b


def test_equal_ndarray():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0])
    b = np.array([11.0, 2.0, 3.0, 4.1, 5.0])
    expected = [False, True, True, False, True]
    assert all((a == b).values == expected)
    a.unit = "m"
    with pytest.raises(DimensionalityError):
        _ = a == b


def test_equal_float():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0])
    b = 3.0
    expected = [False, False, True, False, False]
    assert all((a == b).values == expected)
    a.unit = "m"
    with pytest.raises(DimensionalityError):
        _ = a == b


def test_equal_quantity():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = 3.0 * units("m")
    expected = [False, False, True, False, False]
    assert all((a == b).values == expected)
    b = 3.0 * units("s")
    with pytest.raises(DimensionalityError):
        _ = a == b


def test_not_equal():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = Array(values=[11.0, 2.0, 3.0, 4.1, 5.0], unit="m")
    expected = [True, False, False, True, False]
    assert all((a != b).values == expected)


def test_not_equal_conversion():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = Array(values=[1100.0, 200.0, 300.0, 410.0, 500.0], unit="cm")
    expected = [True, False, False, True, False]
    assert all((a != b).values == expected)


def test_not_equal_bad_units():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = Array(values=[11.0, 2.0, 3.0, 4.1, 5.0], unit="s")
    with pytest.raises(DimensionalityError):
        _ = a != b


def test_not_equal_ndarray():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0])
    b = np.array([11.0, 2.0, 3.0, 4.1, 5.0])
    expected = [True, False, False, True, False]
    assert all((a != b).values == expected)
    a.unit = "m"
    with pytest.raises(DimensionalityError):
        _ = a != b


def test_not_equal_float():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0])
    b = 3.0
    expected = [True, True, False, True, True]
    assert all((a != b).values == expected)
    a.unit = "m"
    with pytest.raises(DimensionalityError):
        _ = a != b


def test_not_equal_quantity():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = 3.0 * units("m")
    expected = [True, True, False, True, True]
    assert all((a != b).values == expected)
    b = 3.0 * units("s")
    with pytest.raises(DimensionalityError):
        _ = a != b


def test_less_than():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="s")
    b = Array(values=[6.0, 7.0, 1.0, 4.0, 10.0], unit="s")
    expected = [True, True, False, False, True]
    assert all((a < b).values == expected)


def test_less_than_conversion():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = Array(values=[600.0, 700.0, 100.0, 400.0, 1000.0], unit="cm")
    expected = [True, True, False, False, True]
    assert all((a < b).values == expected)


def test_less_than_bad_units():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="s")
    b = Array(values=[6.0, 7.0, 1.0, 4.0, 10.0], unit="m")
    with pytest.raises(DimensionalityError):
        _ = a < b


def test_less_than_ndarray():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0])
    b = np.array([6.0, 7.0, 1.0, 4.0, 10.0])
    expected = [True, True, False, False, True]
    assert all((a < b).values == expected)
    a.unit = "m"
    with pytest.raises(DimensionalityError):
        _ = a < b


def test_less_than_float():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0])
    b = 3.0
    expected = [True, True, False, False, False]
    assert all((a < b).values == expected)
    a.unit = "m"
    with pytest.raises(DimensionalityError):
        _ = a < b


def test_less_than_quantity():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = 3.0 * units("m")
    expected = [True, True, False, False, False]
    assert all((a < b).values == expected)
    b = 3.0 * units("s")
    with pytest.raises(DimensionalityError):
        _ = a < b


def test_less_equal():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="s")
    b = Array(values=[6.0, 7.0, 1.0, 4.0, 10.0], unit="s")
    expected = [True, True, False, True, True]
    assert all((a <= b).values == expected)


def test_less_equal_bad_units():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="s")
    b = Array(values=[6.0, 7.0, 1.0, 4.0, 10.0], unit="m")
    with pytest.raises(DimensionalityError):
        _ = a <= b


def test_less_equal_ndarray():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0])
    b = np.array([6.0, 7.0, 1.0, 4.0, 10.0])
    expected = [True, True, False, True, True]
    assert all((a <= b).values == expected)
    a.unit = "m"
    with pytest.raises(DimensionalityError):
        _ = a < b


def test_less_equal_float():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0])
    b = 3.0
    expected = [True, True, True, False, False]
    assert all((a <= b).values == expected)
    a.unit = "m"
    with pytest.raises(DimensionalityError):
        _ = a < b


def test_less_equal_quantity():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = 3.0 * units("m")
    expected = [True, True, True, False, False]
    assert all((a <= b).values == expected)
    b = 3.0 * units("s")
    with pytest.raises(DimensionalityError):
        _ = a < b


def test_greater_than():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="s")
    b = Array(values=[6.0, 7.0, 1.0, 4.0, 10.0], unit="s")
    expected = [True, True, False, False, True]
    assert all((b > a).values == expected)


def test_greater_than_bad_units():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="s")
    b = Array(values=[6.0, 7.0, 1.0, 4.0, 10.0], unit="K")
    with pytest.raises(DimensionalityError):
        _ = b > a


def test_greater_than_ndarray():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0])
    b = np.array([6.0, 7.0, 1.0, 4.0, 10.0])
    expected = [False, False, True, False, False]
    assert all((a > b).values == expected)
    a.unit = "m"
    with pytest.raises(DimensionalityError):
        _ = a > b


def test_greater_than_float():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0])
    b = 3.0
    expected = [False, False, False, True, True]
    assert all((a > b).values == expected)
    a.unit = "m"
    with pytest.raises(DimensionalityError):
        _ = a > b


def test_greater_than_quantity():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = 3.0 * units("m")
    expected = [False, False, False, True, True]
    assert all((a > b).values == expected)
    b = 3.0 * units("s")
    with pytest.raises(DimensionalityError):
        _ = a > b


def test_greater_equal():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="s")
    b = Array(values=[6.0, 7.0, 1.0, 4.0, 10.0], unit="s")
    expected = [True, True, False, True, True]
    assert all((b >= a).values == expected)


def test_greater_equal_bad_units():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="s")
    b = Array(values=[6.0, 7.0, 1.0, 4.0, 10.0], unit="K")
    with pytest.raises(DimensionalityError):
        _ = b >= a


def test_greater_equal_ndarray():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0])
    b = np.array([6.0, 7.0, 1.0, 4.0, 10.0])
    expected = [False, False, True, True, False]
    assert all((a >= b).values == expected)
    a.unit = "m"
    with pytest.raises(DimensionalityError):
        _ = a >= b


def test_greater_equal_float():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0])
    b = 3.0
    expected = [False, False, True, True, True]
    assert all((a >= b).values == expected)
    a.unit = "m"
    with pytest.raises(DimensionalityError):
        _ = a >= b


def test_greater_equal_quantity():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = 3.0 * units("m")
    expected = [False, False, True, True, True]
    assert all((a >= b).values == expected)
    b = 3.0 * units("s")
    with pytest.raises(DimensionalityError):
        _ = a >= b


def test_logical_and():
    a = Array(values=[True, True, True, False, False, False])
    b = Array(values=[True, False, True, False, True, False])
    expected = [True, False, True, False, False, False]
    assert all((b & a).values == expected)


def test_logical_or():
    a = Array(values=[True, True, True, False, False, False])
    b = Array(values=[True, False, True, False, True, False])
    expected = [True, True, True, False, True, False]
    assert all((b | a).values == expected)


def test_logical_xor():
    a = Array(values=[True, True, True, False, False, False])
    b = Array(values=[True, False, True, False, True, False])
    expected = [False, True, False, False, True, False]
    assert all((b ^ a).values == expected)


def test_logical_invert():
    a = Array(values=[True, True, False, False, True, False])
    expected = [False, False, True, True, False, True]
    assert all((~a).values == expected)


def test_to():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    b = Array(values=[1.0e-3, 2.0e-3, 3.0e-3, 4.0e-3, 5.0e-3], unit="km")
    assert arrayclose(a.to("km"), b)
    assert a.unit == units("m")


def test_to_bad_units():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    with pytest.raises(DimensionalityError):
        _ = a.to("s")


def test_min():
    a = Array(values=[1.0, -2.0, 3.0, 0.4, 0.5, 0.6], unit="m")
    assert (a.min() == Array(values=-2.0, unit="m")).values
    b = Array(values=np.array([1.0, -2.0, 3.0, 0.4, 0.5, 0.6]).reshape(2, 3), unit="m")
    assert (b.min() == Array(values=-2.0, unit="m")).values


def test_max():
    a = Array(values=[1.0, 2.0, 3.0, -15.0, 5.0, 6.0], unit="m")
    assert (a.max() == Array(values=6.0, unit="m")).values
    b = Array(values=np.array([1.0, 2.0, 3.0, -15.0, 5.0, 6.0]).reshape(2, 3), unit="m")
    assert (b.max() == Array(values=6.0, unit="m")).values


def test_reshape():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0, 6.0], unit="m")
    expected = Array(values=[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], unit="m")
    assert arrayequal(a.reshape(2, 3), expected)


def test_slicing():
    a = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    assert a[2] == Array(values=[13.0], unit="m")
    assert arrayequal(a[:4], Array(values=[11.0, 12.0, 13.0, 14.0], unit="m"))
    assert arrayequal(a[2:4], Array(values=[13.0, 14.0], unit="m"))


def test_slicing_2d():
    a = Array(values=np.arange(12.0).reshape(4, 3), unit="m")
    assert arrayequal(a[0, :], Array(values=[0.0, 1.0, 2.0], unit="m"))
    assert arrayequal(a[1], Array(values=[3.0, 4.0, 5.0], unit="m"))
    assert arrayequal(a[2:3], Array(values=[[6.0, 7.0, 8.0]], unit="m"))
    assert a[2:3].shape == (1, 3)
    assert arrayequal(a[:2], Array(values=[[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]], unit="m"))


def test_slicing_with_array():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    select = a > (3.0 * units("m"))
    assert arrayequal(a[select], Array(values=[4.0, 5.0], unit="m"))
    b = Array(values=[0, 2, 4])
    assert arrayequal(a[b], Array(values=[1.0, 3.0, 5.0], unit="m"))


def test_slicing_with_array_bad_dtype():
    a = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    with pytest.raises(TypeError):
        _ = a[Array(values=[0.0, 2.0, 4.0])]


def test_copy():
    a = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    b = a.copy()
    a *= 10.0
    assert arraytrue(b == Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m"))


def test_copy_overload():
    a = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    b = copy(a)
    a *= 10.0
    assert arraytrue(b == Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m"))


def test_deepcopy():
    a = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    b = deepcopy(a)
    a *= 10.0
    assert arraytrue(b == Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m"))


def test_numpy_unary():
    values = [1.0, 2.0, 3.0, 4.0, 5.0]
    a = Array(values=values, unit="m")
    expected = np.log10(values)
    result = np.log10(a)
    assert np.allclose(result.values, expected)
    assert result.unit == units("m")


def test_numpy_sqrt():
    values = [1.0, 2.0, 3.0, 4.0, 5.0]
    a = Array(values=values, unit="m*m")
    expected = np.sqrt(values)
    result = np.sqrt(a)
    assert np.allclose(result.values, expected)
    assert result.unit == units("m")


def test_numpy_binary():
    a_buf = [1.0, 2.0, 3.0, 4.0, 5.0]
    b_buf = [6.0, 7.0, 8.0, 9.0, 10.0]
    a = Array(values=a_buf, unit="m")
    b = Array(values=b_buf, unit="m")
    expected = np.dot(a_buf, b_buf)
    result = np.dot(a, b)
    assert result.values == expected
    assert result.unit == units("m")


def test_numpy_iterable():
    a_buf = [1.0, 2.0, 3.0, 4.0, 5.0]
    b_buf = [6.0, 7.0, 8.0, 9.0, 10.0]
    a = Array(values=a_buf, unit="m")
    b = Array(values=b_buf, unit="m")
    expected = np.concatenate([a_buf, b_buf])
    result = np.concatenate([a, b])
    assert np.array_equal(result.values, expected)
    assert result.unit == units("m")


def test_numpy_multiply_with_ndarray():
    a_buf = [1.0, 2.0, 3.0, 4.0, 5.0]
    a = Array(values=a_buf, unit="m")
    b = np.array([6.0, 7.0, 8.0, 9.0, 10.0])
    expected = np.multiply(a_buf, b)
    result = np.multiply(a, b)
    assert np.array_equal(result.values, expected)
    assert result.unit == units("m")
    result = np.multiply(b, a)
    assert np.array_equal(result.values, expected)
    assert result.unit == units("m")


def test_numpy_multiply_with_quantity():
    a_buf = [1.0, 2.0, 3.0, 4.0, 5.0]
    a = Array(values=a_buf, unit="m")
    b = 3.5 * units("s")
    expected = np.multiply(a_buf, b.magnitude)
    result = np.multiply(a, b)
    assert np.array_equal(result.values, expected)
    assert result.unit == units("m*s")


def test_numpy_multiply_with_float():
    a_buf = [1.0, 2.0, 3.0, 4.0, 5.0]
    a = Array(values=a_buf, unit="m")
    b = 3.5
    expected = np.multiply(a_buf, b)
    result = np.multiply(a, b)
    assert np.array_equal(result.values, expected)
    assert result.unit == units("m")
    result = np.multiply(b, a)
    assert np.array_equal(result.values, expected)
    assert result.unit == units("m")


def test_numpy_divide_with_ndarray():
    a_buf = [1.0, 2.0, 3.0, 4.0, 5.0]
    a = Array(values=a_buf, unit="m")
    b = np.array([6.0, 7.0, 8.0, 9.0, 10.0])
    expected = np.divide(a_buf, b)
    result = np.divide(a, b)
    assert np.array_equal(result.values, expected)
    assert result.unit == units("m")
    expected = np.divide(b, a_buf)
    result = np.divide(b, a)
    assert np.array_equal(result.values, expected)
    assert result.unit == units("1/m")


def test_numpy_divide_with_quantity():
    a_buf = [1.0, 2.0, 3.0, 4.0, 5.0]
    a = Array(values=a_buf, unit="m")
    b = 3.5 * units("s")
    expected = np.divide(a_buf, b.magnitude)
    result = np.divide(a, b)
    assert np.array_equal(result.values, expected)
    assert result.unit == units("m/s")


def test_numpy_divide_with_float():
    a_buf = [1.0, 2.0, 3.0, 4.0, 5.0]
    a = Array(values=a_buf, unit="m")
    b = 3.5
    expected = np.divide(a_buf, b)
    result = np.divide(a, b)
    assert np.array_equal(result.values, expected)
    assert result.unit == units("m")
    expected = np.divide(b, a_buf)
    result = np.divide(b, a)
    assert np.array_equal(result.values, expected)
    assert result.unit == units("1/m")
