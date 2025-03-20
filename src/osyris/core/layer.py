# SPDX-License-Identifier: BSD-3-Clause
from collections.abc import Iterable
from typing import Dict, Optional, Union

from .array import Array


class Layer:
    def __init__(
        self,
        data: Array,
        aux: Dict[str, Array] = None,
        *,
        mode: Optional[str] = None,
        operation: Optional[str] = None,
        norm: Optional[str] = None,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        bins: Optional[Union[int, Iterable]] = None,
        weights: Optional[Array] = None,
        **kwargs,
    ):
        self.key = data.name
        self.arrays = {self.key: data}
        if aux:
            self.arrays.update(aux)
        self.mode = mode
        self.operation = operation
        self.norm = norm
        self.vmin = vmin
        self.vmax = vmax
        self.bins = bins
        self.weights = weights
        self.kwargs = kwargs

    def __getitem__(self, key: str) -> Array:
        return self.arrays[key]

    def keys(self) -> Iterable[str]:
        return self.arrays.keys()

    def values(self) -> Iterable[Array]:
        return self.arrays.values()

    def items(self) -> Iterable[Dict[str, Array]]:
        return self.arrays.items()

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return (
            f"Layer({self.data.name}, mode={self.mode}, "
            f"operation={self.operation}, {self.kwargs})"
        )

    def copy(self):
        return self.__class__(
            self.data,
            aux={
                key: value  # .copy()
                for key, value in self.arrays.items()
                if key != self.key
            },
            mode=self.mode,
            operation=self.operation,
            norm=self.norm,
            vmin=self.vmin,
            vmax=self.vmax,
            bins=self.bins,
            weights=self.weights,
            **self.kwargs,
        )

    @property
    def data(self) -> Array:
        return self.arrays[self.key]

    @property
    def x(self):
        out = self.copy()
        out.arrays[self.key] = self.data.x
        return out

    @property
    def y(self):
        out = self.copy()
        out.arrays[self.key] = self.data.y
        return out

    @property
    def z(self):
        out = self.copy()
        out.arrays[self.key] = self.data.z
        return out

    def update(
        self,
        mode: Optional[str] = None,
        operation: Optional[str] = None,
        norm: Optional[str] = None,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        bins: Optional[Union[int, Iterable]] = None,
        weights: Optional[Array] = None,
        **kwargs,
    ):
        if self.mode is None:
            self.mode = mode
        if self.operation is None:
            self.operation = operation
        if self.norm is None:
            self.norm = norm
        if self.vmin is None:
            self.vmin = vmin
        if self.vmax is None:
            self.vmax = vmax
        if self.bins is None:
            self.bins = bins
        if self.weights is None:
            self.weights = weights
        self.kwargs.update(
            {key: value for key, value in kwargs.items() if key not in self.kwargs}
        )
