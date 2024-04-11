# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2024 Osyris contributors (https://github.com/osyris-project/osyris)
from collections.abc import Iterable
from typing import Dict, Optional

from ..core.array import Array


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

    @property
    def data(self) -> Array:
        return self.arrays[self.key]

    def update(
        self,
        mode: Optional[str] = None,
        operation: Optional[str] = None,
        norm: Optional[str] = None,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
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
        # self.norm = get_norm(norm=self.norm, vmin=self.vmin, vmax=self.vmax)
        self.kwargs.update(**kwargs)
