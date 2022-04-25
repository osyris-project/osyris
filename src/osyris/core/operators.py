# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)
import numpy as np
from pint.quantity import Quantity

# def _raise_incompatible_types_error(lhs, rhs, op):
#     raise TypeError("Could not {} types {} and {}.".format(op, lhs, rhs))


def _maybe_broadcast(lhs, rhs):
    if (not rhs.shape) or (lhs.ndim == rhs.ndim):
        return rhs
    # if lhs.ndim != rhs.ndim:
    if hasattr(lhs, "x") and (not hasattr(rhs, "x")):
        rhs = rhs.reshape(rhs.shape + (1, ))
    else:
        return NotImplemented
    # return rhs


def _to_lhs_unit(lhs, rhs):
    return rhs.unit.to(lhs._unit.units)


def _to_base_units(op, lhs, rhs):
    scale_l = lhs._unit.to_base_units()
    scale_r = rhs._unit.to_base_units()
    return op(scale_l, scale_r)


# def _add_subtract(op, lhs, rhs, out=None):
#     result = None
#     if isinstance(rhs, Quantity):
#         rhs = lhs.__class__(values=rhs.magnitude, unit=1.0 * rhs.units)
#         # print("here", rhs)
#     if not isinstance(rhs, lhs.__class__):
#         # print("bad")
#         return NotImplemented
#     rhs = _maybe_broadcast(lhs, rhs)
#     # print("after boradcast", rhs)
#     if rhs is NotImplemented:
#         return rhs
#     # print("two", rhs)
#     ratio = _to_lhs_unit(lhs, rhs)
#     # print("ratio", ratio)
#     # rhs = rhs._array * ratio
#     result = op(lhs._array, rhs._array, out=out)
#     result *= ratio
#     # # return lhs.__class__(values=op(lhs._array, rhs), unit=lhs._unit)
#     # if isinstance(rhs, Quantity):
#     #     result = op(lhs._array, rhs._array, out=out)

#     #     return lhs.__class__(values=op(lhs._array,
#     #                                    rhs.to(lhs._unit.units).magnitude),
#     #                          unit=lhs._unit)
#     # if result is not None:
#     if out is None:
#         return lhs.__class__(values=result, unit=lhs._unit)
#     else:
#         return lhs

#     # _raise_incompatible_types_error(lhs, rhs, op.__name__)

# def _multiply_divide(op, lhs, rhs, out=None):
#     result = None
#     if isinstance(rhs, Quantity):
#         rhs = lhs.__class__(values=rhs.magnitude, unit=1.0 * rhs.units)
#         # print("here", rhs)
#     if not isinstance(rhs, lhs.__class__):
#         # print("bad")
#         return NotImplemented
#     rhs = _maybe_broadcast(lhs, rhs)
#     # print("after boradcast", rhs)
#     if rhs is NotImplemented:
#         return rhs
#     # print("two", rhs)
#     ratio = _to_base_units(op, lhs, rhs)
#     # print("ratio", ratio)
#     # rhs = rhs._array * ratio
#     result = op(lhs._array, rhs._array, out=out)
#     result *= ratio.magnitude
#     # # return lhs.__class__(values=op(lhs._array, rhs), unit=lhs._unit)
#     # if isinstance(rhs, Quantity):
#     #     result = op(lhs._array, rhs._array, out=out)

#     #     return lhs.__class__(values=op(lhs._array,
#     #                                    rhs.to(lhs._unit.units).magnitude),
#     #                          unit=lhs._unit)
#     # if result is not None:
#     if out is None:
#         return lhs.__class__(values=result, unit=1.0 * ratio.units)
#     else:
#         return lhs


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
    ratio = _to_base_units(op, lhs, rhs) if to_base_units else _to_lhs_unit(lhs, rhs)
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
