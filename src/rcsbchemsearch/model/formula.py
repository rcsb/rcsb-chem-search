# SPDX-FileCopyrightText: Copyright 2024, Contributors to rcsb-chem
# SPDX-PackageHomePage: https://github.com/rcsb/rcsb-chem
# SPDX-License-Identifier: BSD-3-Clause

"""
Models.
"""
import re
from enum import auto, StrEnum, Flag
from pydantic import BaseModel, ConfigDict
from typing import Self, TypeVar


__all__ = [
    "Formula",
    "ElementReq",
]
_maximum = 1_000_000_000_000_000  # 1 quadrillion
_T = TypeVar("_T")


class ElementReq(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid", strict=True)
    elements: set[str]
    minimum: int = 0
    maximum: int = _maximum

    def __str__(self: Self) -> str:
        return repr(self)

    def __repr__(self: Self) -> str:
        s = ""
        if self.minimum > 0:
            s += f"{self.minimum}<="
        if len(self.elements) == 1:
            s += next(iter(self.elements))
        else:
            s += f"({"|".join(self.elements)})"
        if self.maximum < _maximum:
            s += f"<={self.maximum}"
        return s

    def __contains__(self: Self, element: str) -> bool:
        return element in self.elements

    def matches(self: Self, element: str, count: int) -> bool:
        return element in self.elements and self.minimum <= count <= self.maximum


class Formula(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid", strict=True)
    elements: list[ElementReq]

    # B, C, N, O, P, S, F, Cl, Br, or I,
    def matches(self: Self, formula: str) -> bool:
        re.compile(r"\(+\)")
