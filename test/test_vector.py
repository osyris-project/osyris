# SPDX-License-Identifier: BSD-3-Clause
from copy import copy, deepcopy

import numpy as np
import pytest
from common import arrayclose, arrayequal, vectorclose, vectorequal

from osyris import Array, Vector, units


def test_constructor_from_arrays():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    assert len(v) == len(x)
    assert v.shape == x.shape


def test_constructor_bad_length():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0], unit="m")
    with pytest.raises(ValueError):
        _ = Vector(x, y, z)


def test_constructor_bad_unit():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="cm")
    with pytest.raises(ValueError):
        _ = Vector(x, y, z)


def test_components():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y)
    assert arrayequal(v.x, x)
    assert arrayequal(v.y, y)
    assert v.z is None
    v.z = z
    assert arrayequal(v.z, z)
    assert v.x.name == "_x"
    assert v.y.name == "_y"


def test_bad_constructor_from_arrays_with_unit():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0])
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0])
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0])
    with pytest.raises(ValueError):
        _ = Vector(x, y, z, unit="K")


def test_constructor_from_floats():
    v = Vector(1.0, 2.0, 3.0, unit="s")
    assert len(v) == 0
    exp_x = Array(1.0, unit="s")
    exp_y = Array(2.0, unit="s")
    exp_z = Array(3.0, unit="s")
    assert arrayequal(v.x, exp_x)
    assert arrayequal(v.y, exp_y)
    assert arrayequal(v.z, exp_z)


def test_constructor_from_ndarrays():
    x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    y = np.array([6.0, 7.0, 8.0, 9.0, 10.0])
    z = np.array([11.0, 12.0, 13.0, 14.0, 15.0])
    v = Vector(x, y, z, unit="s")
    assert len(v) == 5
    assert v.shape == (5,)
    exp_x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="s")
    exp_y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="s")
    exp_z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="s")
    assert arrayequal(v.x, exp_x)
    assert arrayequal(v.y, exp_y)
    assert arrayequal(v.z, exp_z)


def test_addition():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    w = Vector(y, x, z)
    expected = Vector(x=x + y, y=y + x, z=z + z)
    assert vectorequal(v + w, expected)
    assert vectorequal(w + v, expected)


def test_addition_array():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    expected = Vector(x=x + x, y=y + x, z=z + x)
    assert vectorequal(v + x, expected)
    assert vectorequal(x + v, expected)


def test_addition_quantity():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    q = 3.5 * units("m")
    expected = Vector(x=x + q, y=y + q, z=z + q)
    assert vectorequal(v + q, expected)
    assert vectorequal(q + v, expected)


def test_addition_inplace():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    expected = Vector(x=x + x, y=y + y, z=z + z)
    v += v
    assert vectorequal(v, expected)


def test_addition_array_inplace():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    a = Array(values=[1.1, 2.2, 3.3, 4.4, 5.5], unit="m")
    expected = Vector(x=x + a, y=y + a, z=z + a)
    v += a
    assert vectorequal(v, expected)


def test_addition_quantity_inplace():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    q = 3.5 * units("m")
    expected = Vector(x=x + q, y=y + q, z=z + q)
    v += q
    assert vectorequal(v, expected)


def test_subtraction():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    w = Vector(y, x, z)
    expected = Vector(x=x - y, y=y - x, z=z - z)
    assert vectorequal(v - w, expected)
    assert vectorequal(w - v, -expected)


def test_subtraction_array():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    expected = Vector(x=x - x, y=y - x, z=z - x)
    assert vectorequal(v - x, expected)
    assert vectorequal(x - v, -expected)


def test_subtraction_quantity():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    q = 3.5 * units("m")
    expected = Vector(x=x - q, y=y - q, z=z - q)
    assert vectorequal(v - q, expected)
    assert vectorequal(q - v, -expected)


def test_subtraction_inplace():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    expected = Vector(x=x - x, y=y - y, z=z - z)
    v -= v
    assert vectorequal(v, expected)


def test_subtraction_array_inplace():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    a = Array(values=[1.1, 2.2, 3.3, 4.4, 5.5], unit="m")
    expected = Vector(x=x - a, y=y - a, z=z - a)
    v -= a
    assert vectorequal(v, expected)


def test_subtraction_quantity_inplace():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    q = 3.5 * units("m")
    expected = Vector(x=x - q, y=y - q, z=z - q)
    v -= q
    assert vectorequal(v, expected)


def test_multiplication():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    w = Vector(z, x, y)
    expected = Vector(x=x * z, y=y * x, z=z * y)
    assert vectorequal(v * w, expected)
    assert vectorequal(w * v, expected)


def test_multiplication_array():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    expected = Vector(x=x * x, y=y * x, z=z * x)
    assert vectorequal(v * x, expected)
    assert vectorequal(x * v, expected)


def test_multiplication_float():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    f = 3.5
    expected = Vector(x=x * f, y=y * f, z=z * f)
    assert vectorequal(v * f, expected)
    assert vectorequal(f * v, expected)


def test_multiplication_ndarray():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    a = np.arange(5.0)
    expected = Vector(x=x * a, y=y * a, z=z * a)
    assert vectorequal(v * a, expected)
    assert vectorequal(a * v, expected)


def test_multiplication_quantity():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    q = 3.5 * units("m")
    expected = Vector(x=x * q, y=y * q, z=z * q)
    assert vectorequal(v * q, expected)
    assert vectorequal(q * v, expected)


def test_multiplication_inplace():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    expected = Vector(x=x * x, y=y * y, z=z * z)
    v *= v
    assert vectorequal(v, expected)


def test_multiplication_array_inplace():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    a = Array(values=[1.1, 2.2, 3.3, 4.4, 5.5], unit="m")
    expected = Vector(x=x * a, y=y * a, z=z * a)
    v *= a
    assert vectorequal(v, expected)


def test_multiplication_float_inplace():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    f = 3.5
    expected = Vector(x=x * f, y=y * f, z=z * f)
    v *= f
    assert vectorequal(v, expected)


def test_multiplication_ndarray_inplace():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    a = np.arange(5.0)
    expected = Vector(x=x * a, y=y * a, z=z * a)
    v *= a
    assert vectorequal(v, expected)


def test_multiplication_quantity_inplace():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    q = 3.5 * units("m")
    expected = Vector(x=x * q, y=y * q, z=z * q)
    v *= q
    assert vectorequal(v, expected)


def test_division():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v1 = Vector(x, y, z)
    w = Array(values=[16.0, 17.0, 18.0, 19.0, 20.0], unit="m")
    v2 = Vector(w, y, w)
    expected = Vector(x=x / w, y=y / y, z=z / w)
    assert vectorequal(v1 / v2, expected)
    expected = Vector(x=w / x, y=y / y, z=w / z)
    assert vectorequal(v2 / v1, expected)


def test_division_array():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    expected = Vector(x=x / x, y=y / x, z=z / x)
    assert vectorequal(v / x, expected)
    expected = Vector(x=x / x, y=x / y, z=x / z)
    assert vectorequal(x / v, expected)


def test_division_float():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    f = 3.5
    expected = Vector(x=x / f, y=y / f, z=z / f)
    assert vectorequal(v / f, expected)
    expected = Vector(x=f / x, y=f / y, z=f / z)
    assert vectorequal(f / v, expected)


def test_division_ndarray():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    a = np.arange(3.0, 8.0)
    expected = Vector(x=x / a, y=y / a, z=z / a)
    assert vectorequal(v / a, expected)
    expected = Vector(x=a / x, y=a / y, z=a / z)
    assert vectorequal(a / v, expected)


def test_division_quantity():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    q = 3.5 * units("m")
    expected = Vector(x=x / q, y=y / q, z=z / q)
    assert vectorequal(v / q, expected)
    expected = Vector(x=q / x, y=q / y, z=q / z)
    assert vectorequal(q / v, expected)


def test_division_inplace():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v1 = Vector(x, y, z)
    w = Array(values=[16.0, 17.0, 18.0, 19.0, 20.0], unit="m")
    v2 = Vector(w, y, w)
    expected = Vector(x=x / w, y=y / y, z=z / w)
    v1 /= v2
    assert vectorequal(v1, expected)


def test_division_array_inplace():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    a = Array(values=[1.1, 2.2, 3.3, 4.4, 5.5], unit="m")
    expected = Vector(x=x / a, y=y / a, z=z / a)
    v /= a
    assert vectorequal(v, expected)


def test_division_float_inplace():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    f = 3.5
    expected = Vector(x=x / f, y=y / f, z=z / f)
    v /= f
    assert vectorequal(v, expected)


def test_division_ndarray_inplace():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    a = np.arange(3.0, 8.0)
    expected = Vector(x=x / a, y=y / a, z=z / a)
    v /= a
    assert vectorequal(v, expected)


def test_division_quantity_inplace():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    q = 3.5 * units("m")
    expected = Vector(x=x / q, y=y / q, z=z / q)
    v /= q
    assert vectorequal(v, expected)


def test_norm():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="s")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="s")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="s")
    v = Vector(x, y)
    assert arrayclose(
        v.norm, Array(values=np.sqrt([37.0, 53.0, 73.0, 97.0, 125.0]), unit="s")
    )
    v.z = z
    assert arrayclose(
        v.norm, Array(values=np.sqrt([158.0, 197.0, 242.0, 293.0, 350.0]), unit="s")
    )


def test_power():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    expected = Vector(x=x**3, y=y**3, z=z**3)
    assert vectorequal(v**3, expected)


def test_less_than():
    x1 = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y1 = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z1 = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v1 = Vector(x1, y1, z1)
    x2 = Array(values=[-1.0, 2.0, 13.0, 4.0, 2.0], unit="m")
    y2 = Array(values=[6.0, -7.0, 8.0, 3.3, 10.1], unit="m")
    z2 = Array(values=[11.0, 22.0, 7.0, 7.0, 15.0], unit="m")
    v2 = Vector(x2, y2, z2)
    exp_x = [False, False, True, False, False]
    exp_y = [False, False, False, False, True]
    exp_z = [False, True, False, False, False]
    result = v1 < v2
    assert all(result.x.values == exp_x)
    assert all(result.y.values == exp_y)
    assert all(result.z.values == exp_z)


def test_less_than_array():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    a = Array(values=[-1.0, 2.0, 17.0, 4.0, 2.0], unit="m")
    exp_x = [False, False, True, False, False]
    exp_y = [False, False, True, False, False]
    exp_z = [False, False, True, False, False]
    result = v < a
    assert all(result.x.values == exp_x)
    assert all(result.y.values == exp_y)
    assert all(result.z.values == exp_z)


def test_less_equal():
    x1 = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y1 = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z1 = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v1 = Vector(x1, y1, z1)
    x2 = Array(values=[-1.0, 2.0, 13.0, 4.0, 2.0], unit="m")
    y2 = Array(values=[6.0, -7.0, 8.0, 3.3, 10.1], unit="m")
    z2 = Array(values=[11.0, 22.0, 7.0, 7.0, 15.0], unit="m")
    v2 = Vector(x2, y2, z2)
    exp_x = [False, True, True, True, False]
    exp_y = [True, False, True, False, True]
    exp_z = [True, True, False, False, True]
    result = v1 <= v2
    assert all(result.x.values == exp_x)
    assert all(result.y.values == exp_y)
    assert all(result.z.values == exp_z)


def test_equal():
    x1 = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y1 = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z1 = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v1 = Vector(x1, y1, z1)
    x2 = Array(values=[-1.0, 2.0, 13.0, 4.0, 2.0], unit="m")
    y2 = Array(values=[6.0, -7.0, 8.0, 3.3, 10.1], unit="m")
    z2 = Array(values=[11.0, 22.0, 7.0, 7.0, 15.0], unit="m")
    v2 = Vector(x2, y2, z2)
    exp_x = [False, True, False, True, False]
    exp_y = [True, False, True, False, False]
    exp_z = [True, False, False, False, True]
    result = v1 == v2
    assert all(result.x.values == exp_x)
    assert all(result.y.values == exp_y)
    assert all(result.z.values == exp_z)


def test_not_equal():
    x1 = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y1 = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z1 = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v1 = Vector(x1, y1, z1)
    x2 = Array(values=[-1.0, 2.0, 13.0, 4.0, 2.0], unit="m")
    y2 = Array(values=[6.0, -7.0, 8.0, 3.3, 10.1], unit="m")
    z2 = Array(values=[11.0, 22.0, 7.0, 7.0, 15.0], unit="m")
    v2 = Vector(x2, y2, z2)
    exp_x = [True, False, True, False, True]
    exp_y = [False, True, False, True, True]
    exp_z = [False, True, True, True, False]
    result = v1 != v2
    assert all(result.x.values == exp_x)
    assert all(result.y.values == exp_y)
    assert all(result.z.values == exp_z)


def test_greater_than():
    x1 = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y1 = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z1 = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v1 = Vector(x1, y1, z1)
    x2 = Array(values=[-1.0, 2.0, 13.0, 4.0, 2.0], unit="m")
    y2 = Array(values=[6.0, -7.0, 8.0, 3.3, 10.1], unit="m")
    z2 = Array(values=[11.0, 22.0, 7.0, 7.0, 15.0], unit="m")
    v2 = Vector(x2, y2, z2)
    exp_x = [True, False, False, False, True]
    exp_y = [False, True, False, True, False]
    exp_z = [False, False, True, True, False]
    result = v1 > v2
    assert all(result.x.values == exp_x)
    assert all(result.y.values == exp_y)
    assert all(result.z.values == exp_z)


def test_greater_equal():
    x1 = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y1 = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z1 = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v1 = Vector(x1, y1, z1)
    x2 = Array(values=[-1.0, 2.0, 13.0, 4.0, 2.0], unit="m")
    y2 = Array(values=[6.0, -7.0, 8.0, 3.3, 10.1], unit="m")
    z2 = Array(values=[11.0, 22.0, 7.0, 7.0, 15.0], unit="m")
    v2 = Vector(x2, y2, z2)
    exp_x = [True, True, False, True, True]
    exp_y = [True, True, True, True, False]
    exp_z = [True, False, True, True, True]
    result = v1 >= v2
    assert all(result.x.values == exp_x)
    assert all(result.y.values == exp_y)
    assert all(result.z.values == exp_z)


def test_logical_and():
    x1 = Array(values=[True, True, True, False, False])
    y1 = Array(values=[True, False, True, False, True])
    z1 = Array(values=[False, False, True, True, True])
    v1 = Vector(x1, y1, z1)
    x2 = Array(values=[False, True, False, True, False])
    y2 = Array(values=[True, True, False, False, False])
    z2 = Array(values=[True, False, False, False, True])
    v2 = Vector(x2, y2, z2)
    exp_x = [False, True, False, False, False]
    exp_y = [True, False, False, False, False]
    exp_z = [False, False, False, False, True]
    result = v1 & v2
    assert all(result.x.values == exp_x)
    assert all(result.y.values == exp_y)
    assert all(result.z.values == exp_z)


def test_logical_or():
    x1 = Array(values=[True, True, True, False, False])
    y1 = Array(values=[True, False, True, False, True])
    z1 = Array(values=[False, False, True, True, True])
    v1 = Vector(x1, y1, z1)
    x2 = Array(values=[False, True, False, True, False])
    y2 = Array(values=[True, True, False, False, False])
    z2 = Array(values=[True, False, False, False, True])
    v2 = Vector(x2, y2, z2)
    exp_x = [True, True, True, True, False]
    exp_y = [True, True, True, False, True]
    exp_z = [True, False, True, True, True]
    result = v1 | v2
    assert all(result.x.values == exp_x)
    assert all(result.y.values == exp_y)
    assert all(result.z.values == exp_z)


def test_logical_xor():
    x1 = Array(values=[True, True, True, False, False])
    y1 = Array(values=[True, False, True, False, True])
    z1 = Array(values=[False, False, True, True, True])
    v1 = Vector(x1, y1, z1)
    x2 = Array(values=[False, True, False, True, False])
    y2 = Array(values=[True, True, False, False, False])
    z2 = Array(values=[True, False, False, False, True])
    v2 = Vector(x2, y2, z2)
    exp_x = [True, False, True, True, False]
    exp_y = [False, True, True, False, True]
    exp_z = [True, False, True, True, False]
    result = v1 ^ v2
    assert all(result.x.values == exp_x)
    assert all(result.y.values == exp_y)
    assert all(result.z.values == exp_z)


def test_logical_invert():
    x = Array(values=[True, True, True, False, False])
    y = Array(values=[True, False, True, False, True])
    z = Array(values=[False, False, True, True, True])
    v = Vector(x, y, z)
    exp_x = [False, False, False, True, True]
    exp_y = [False, True, False, True, False]
    exp_z = [True, True, False, False, False]
    result = ~v
    assert all(result.x.values == exp_x)
    assert all(result.y.values == exp_y)
    assert all(result.z.values == exp_z)


def test_to():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    km = "km"
    expected = Vector(x=x.to(km), y=y.to(km), z=z.to(km))
    assert vectorclose(v.to(km), expected)


def test_min():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="s")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="s")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="s")
    v = Vector(x, y, z)
    expected = Vector(x=x.min(), y=y.min(), z=z.min())
    assert vectorequal(v.min(), expected)


def test_max():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="s")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="s")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="s")
    v = Vector(x, y, z)
    expected = Vector(x=x.max(), y=y.max(), z=z.max())
    assert vectorequal(v.max(), expected)


def test_reshape():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0, 6.0], unit="s")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0, 11.0], unit="s")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0, 16.0], unit="s")
    v = Vector(x, y, z)
    shape = (2, 3)
    expected = Vector(x=x.reshape(*shape), y=y.reshape(*shape), z=z.reshape(*shape))
    assert vectorequal(v.reshape(*shape), expected)


def test_slicing():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    assert vectorequal(
        v[2],
        Vector(
            x=Array(values=3.0, unit="m"),
            y=Array(values=8.0, unit="m"),
            z=Array(values=13.0, unit="m"),
        ),
    )
    assert vectorequal(
        v[:4],
        Vector(
            x=Array(values=[1.0, 2.0, 3.0, 4.0], unit="m"),
            y=Array(values=[6.0, 7.0, 8.0, 9.0], unit="m"),
            z=Array(values=[11.0, 12.0, 13.0, 14.0], unit="m"),
        ),
    )
    assert vectorequal(
        v[2:4],
        Vector(
            x=Array(values=[3.0, 4.0], unit="m"),
            y=Array(values=[8.0, 9.0], unit="m"),
            z=Array(values=[13.0, 14.0], unit="m"),
        ),
    )
    assert vectorequal(v[2:4], Vector(x=x[2:4], y=y[2:4], z=z[2:4]))


def test_slicing_with_array():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    select = v.norm > (15.0 * units("m"))
    exp_x = Array(values=[3.0, 4.0, 5.0], unit="m")
    exp_y = Array(values=[8.0, 9.0, 10.0], unit="m")
    exp_z = Array(values=[13.0, 14.0, 15.0], unit="m")
    expected = Vector(exp_x, exp_y, exp_z)
    assert vectorequal(v[select], expected)


def test_slicing_with_vector_raises():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    select = v > (16.0 * units("m"))
    with pytest.raises(ValueError):
        _ = v[select]


def test_copy():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    a = Vector(x, y, z)
    b = a.copy()
    a *= 10.0
    original = Vector(
        x=Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m"),
        y=Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m"),
        z=Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m"),
    )
    assert vectorequal(b, original)


def test_copy_overload():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    a = Vector(x, y, z)
    b = copy(a)
    a *= 10.0
    original = Vector(
        x=Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m"),
        y=Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m"),
        z=Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m"),
    )
    assert vectorequal(b, original)


def test_vector_values():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    val = v.values
    expected = np.array(
        [
            [1.0, 6.0, 11.0],
            [2.0, 7.0, 12.0],
            [3.0, 8.0, 13.0],
            [4.0, 9.0, 14.0],
            [5.0, 10.0, 15.0],
        ]
    )
    assert val.shape == (5, 3)
    assert np.array_equal(val, expected)
    assert np.array_equal(val[:, 0], x.values)
    assert np.array_equal(val[:, 1], y.values)
    assert np.array_equal(val[:, 2], z.values)


def test_deepcopy():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    a = Vector(x, y, z)
    b = deepcopy(a)
    a *= 10.0
    original = Vector(
        x=Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m"),
        y=Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m"),
        z=Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m"),
    )
    assert vectorequal(b, original)


def test_numpy_unary():
    values_x = [1.0, 2.0, 3.0, 4.0, 5.0]
    values_y = [6.0, 7.0, 8.0, 9.0, 10.0]
    values_z = [11.0, 12.0, 13.0, 14.0, 15.0]
    x = Array(values=values_x, unit="m")
    y = Array(values=values_y, unit="m")
    z = Array(values=values_z, unit="m")
    v = Vector(x, y, z)
    exp_x = np.log10(values_x)
    exp_y = np.log10(values_y)
    exp_z = np.log10(values_z)
    result = np.log10(v)
    assert np.allclose(result.x.values, exp_x)
    assert np.allclose(result.y.values, exp_y)
    assert np.allclose(result.z.values, exp_z)
    assert result.unit == units("m")


def test_numpy_sqrt():
    values_x = [1.0, 2.0, 3.0, 4.0, 5.0]
    values_y = [6.0, 7.0, 8.0, 9.0, 10.0]
    values_z = [11.0, 12.0, 13.0, 14.0, 15.0]
    x = Array(values=values_x, unit="m*m")
    y = Array(values=values_y, unit="m*m")
    z = Array(values=values_z, unit="m*m")
    v = Vector(x, y, z)
    exp_x = np.sqrt(values_x)
    exp_y = np.sqrt(values_y)
    exp_z = np.sqrt(values_z)
    result = np.sqrt(v)
    assert np.allclose(result.x.values, exp_x)
    assert np.allclose(result.y.values, exp_y)
    assert np.allclose(result.z.values, exp_z)
    assert result.unit == units("m")


def test_numpy_binary():
    a_values_x = [1.0, 2.0, 3.0, 4.0, 5.0]
    a_values_y = [6.0, 7.0, 8.0, 9.0, 10.0]
    a_values_z = [11.0, 12.0, 13.0, 14.0, 15.0]
    b_values_x = [-1.0, -2.0, -3.0, -4.0, 5.0]
    b_values_y = [6.1, 7.2, 8.3, 9.4, 10.5]
    b_values_z = [111.0, 112.0, 113.0, -114.0, 115.0]
    a_x = Array(values=a_values_x, unit="m")
    a_y = Array(values=a_values_y, unit="m")
    a_z = Array(values=a_values_z, unit="m")
    b_x = Array(values=b_values_x, unit="m")
    b_y = Array(values=b_values_y, unit="m")
    b_z = Array(values=b_values_z, unit="m")
    a = Vector(a_x, a_y, a_z)
    b = Vector(b_x, b_y, b_z)
    exp_x = np.add(a_values_x, b_values_x)
    exp_y = np.add(a_values_y, b_values_y)
    exp_z = np.add(a_values_z, b_values_z)
    result = np.add(a, b)
    assert np.allclose(result.x.values, exp_x)
    assert np.allclose(result.y.values, exp_y)
    assert np.allclose(result.z.values, exp_z)
    assert result.unit == units("m")


def test_numpy_vstack():
    x = Array(values=[1.0, 2.0, 3.0, 4.0, 5.0], unit="m")
    y = Array(values=[6.0, 7.0, 8.0, 9.0, 10.0], unit="m")
    z = Array(values=[11.0, 12.0, 13.0, 14.0, 15.0], unit="m")
    v = Vector(x, y, z)
    w = Vector(y, x, z)
    result = np.vstack([v, w])
    assert vectorequal(result[0], v)
    assert vectorequal(result[1], w)
