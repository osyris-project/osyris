# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)
import numpy as np
from pint.quantity import Quantity


def _maybe_broadcast(lhs, rhs):
    if (not rhs.shape) or (lhs.ndim == rhs.ndim):
        return rhs
    if hasattr(lhs, "x") and (not hasattr(rhs, "x")):
        rhs = rhs.reshape(rhs.shape + (1, ))
    else:
        return NotImplemented


def _operator(op, lhs, rhs, to_base_units=False, out=None):
    result = None
    if isinstance(rhs, Quantity):
        rhs = lhs.__class__(values=rhs.magnitude, unit=1.0 * rhs.units)
    if isinstance(rhs, (int, float)):
        rhs = lhs.__class__(values=rhs)
    if not isinstance(rhs, lhs.__class__):
        return NotImplemented
    rhs = _maybe_broadcast(lhs, rhs)
    if rhs is NotImplemented:
        return rhs
    ratio = op(lhs._unit, rhs._unit) if to_base_units else rhs.unit.to(lhs._unit.units)
    result = op(lhs._array, rhs._array, out=out)
    result *= ratio.magnitude
    if out is None:
        return lhs.__class__(values=result, unit=1.0 * ratio.units)
    else:
        if to_base_units:
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
    return _operator(np.multiply, lhs, rhs, to_base_units=True)


def imul(lhs, rhs):
    return _operator(np.multiply, lhs, rhs, to_base_units=True, out=lhs._array)


def div(lhs, rhs):
    return _operator(np.divide, lhs, rhs, to_base_units=True)


def idiv(lhs, rhs):
    return _operator(np.divide, lhs, rhs, to_base_units=True, out=lhs._array)


def comp(lhs, rhs, op):
    func = getattr(np, op)
    if isinstance(rhs, Array):
        scale_r = rhs.unit.to(lhs._unit.units)
        return func(lhs._array, rhs._array * scale_r.magnitude)
    if isinstance(rhs, Quantity):
        return func(lhs._array, rhs.to(lhs._unit.units).magnitude)
    return func(lhs._array, rhs)
