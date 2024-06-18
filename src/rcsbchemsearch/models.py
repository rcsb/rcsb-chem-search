# SPDX-FileCopyrightText: Copyright 2024, the RCSB PDB and contributors
# SPDX-PackageHomePage: https://github.com/rcsb/rcsb-chem-search
# SPDX-License-Identifier: BSD-3-Clause
#
"""
Models.
"""
from enum import auto, StrEnum
from dataclasses import dataclass
from typing import Self

__all__ = [
    "Consider",
    "Formula",
    "Metric",
    "AlgorithmParams",
    "SearchType",
]


class SearchType(StrEnum):
    """"""
    fingerprint_ecfp = auto()
    # ...


class Metric(StrEnum):
    """"""
    dice = auto()
    hamming = auto()
    # ...


@dataclass(frozen=True, slots=True)
class AlgorithmParams:
    """"""
    values: dict[str, str]


@dataclass(frozen=True, slots=True, order=True)
class ElementReq:
    element: str
    minimum: int | None
    maximum: int | None


@dataclass(frozen=True, slots=True, order=True)
class Formula:
    """"""
    text: str

    @property
    def allows_more(self: Self) -> bool:
        return False

    @property
    def elements(self: Self) -> list[ElementReq]:
        return []


@dataclass(frozen=True, slots=True, order=True)
class Consider:
    """"""
    optical: bool = True
    cis_trans: bool = True
    tautomer: bool = True
    protonation: bool = True
    hydrogen: bool = True
