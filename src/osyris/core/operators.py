# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)
import numpy as np
from pint.quantity import Quantity


def _maybe_broadcast(lhs, rhs):
    print(type(lhs), type(rhs))
    print("_maybe_broadcast 1", lhs, rhs)
    if (not rhs.shape) or (lhs.ndim == rhs.ndim):
        print("_maybe_broadcast 2")
        return rhs
    print("_maybe_broadcast 3")
    if hasattr(lhs, "x") and (not hasattr(rhs, "x")):
        print("_maybe_broadcast 4")
        return rhs.reshape(rhs.shape + (1, ))
    else:
        print("_maybe_broadcast 5")
        return NotImplemented


def _maybe_broadcast(lhs, rhs):
    if (lhs.ndim == rhs.ndim) or (len(lhs.shape) == 0) or (len(rhs.shape) == 0):
        return lhs, rhs
    if lhs.ndim > rhs.ndim:
        ind = np.argmax(np.array(lhs.shape) == rhs.shape[0])
        if ind == 0:
            return lhs, rhs.reshape(rhs.shape + tuple([1]))
        else:
            return lhs, rhs.reshape(tuple([1]) + rhs.shape)
    else:
        ind = np.argmax(np.array(rhs.shape) == lhs.shape[0])
        if ind == 0:
            return lhs.reshape(lhs.shape + tuple([1])), rhs
        else:
            return lhs.reshape(tuple([1]) + lhs.shape), rhs


def _operator(op, lhs, rhs, mul_or_div=False, out=None):
    result = None
    if isinstance(rhs, Quantity):
        rhs = lhs.__class__(values=rhs.magnitude, unit=1.0 * rhs.units)
    if isinstance(rhs, (int, float)):
        rhs = lhs.__class__(values=rhs)
    # if not isinstance(rhs, lhs.__class__):
    #     return NotImplemented
    rhs = _maybe_broadcast(lhs, rhs)
    if rhs is NotImplemented:
        print("NotImplemented!")
        return rhs
    ratio = op(lhs._unit, rhs._unit) if mul_or_div else rhs.unit.to(lhs._unit.units)
    result = op(lhs._array, rhs._array, out=out)
    result *= ratio.magnitude
    if out is None:
        return lhs.__class__(values=result, unit=1.0 * ratio.units)
    else:
        if mul_or_div:
            lhs._unit = 1.0 * ratio.units
        return lhs


def add(lhs, rhs):
    return _operator(np.add, lhs, rhs)


def iadd(lhs, rhs):
    return _operator(np.add, lhs, rhs, out=lhs._array)


def sub(lhs, rhs):
    return _operator(np.subtract, lhs, rhs)


def isub(lhs, rhs):
    return _operator(np.subtract, lhs, rhs, out=lhs._array)


def mul(lhs, rhs):
    return _operator(np.multiply, lhs, rhs, mul_or_div=True)


def imul(lhs, rhs):
    return _operator(np.multiply, lhs, rhs, mul_or_div=True, out=lhs._array)


def div(lhs, rhs):
    return _operator(np.divide, lhs, rhs, mul_or_div=True)


def idiv(lhs, rhs):
    return _operator(np.divide, lhs, rhs, mul_or_div=True, out=lhs._array)


def comp(op, lhs, rhs):
    if isinstance(rhs, lhs.__class__):
        scale_r = rhs.unit.to(lhs._unit.units)
        return op(lhs._array, rhs._array * scale_r.magnitude)
    if isinstance(rhs, Quantity):
        return op(lhs._array, rhs.to(lhs._unit.units).magnitude)
    return op(lhs._array, rhs)


def binary_op(op, lhs, rhs, mul_or_div=False, out=None):
    if isinstance(rhs, Quantity):
        rhs = lhs.__class__(values=rhs.magnitude, unit=1.0 * rhs.units)
    if isinstance(rhs, (int, float, np.ndarray)):
        rhs = lhs.__class__(values=rhs)
    if (len(rhs) > 1) and (lhs.ndim != rhs.ndim):
        return NotImplemented

    ratio = op(lhs._unit, rhs._unit) if mul_or_div else rhs.unit.to(lhs._unit.units)
    result = op(lhs._array, rhs._array, out=out)
    result *= ratio.magnitude
    if out is None:
        return lhs.__class__(values=result, unit=1.0 * ratio.units)
    else:
        if mul_or_div:
            lhs._unit = 1.0 * ratio.units
        return lhs


def comparison_op(op, lhs, rhs):
    if isinstance(rhs, lhs.__class__):
        scale_r = rhs.unit.to(lhs._unit.units)
        return op(lhs._array, rhs._array * scale_r.magnitude)
    if isinstance(rhs, Quantity):
        return op(lhs._array, rhs.to(lhs._unit.units).magnitude)
    return op(lhs._array, rhs)
