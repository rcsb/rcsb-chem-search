# SPDX-FileCopyrightText: Copyright 2024, Contributors to rcsb-chem-search
# SPDX-PackageHomePage: https://github.com/rcsb/rcsb-chem-search
# SPDX-License-Identifier: BSD-3-Clause

"""
Removing info.
"""
from copy import copy
from dataclasses import dataclass, field
from enum import auto, Enum
from functools import partial
from typing import Self, Callable, ClassVar
from rdkit import Chem
from rdkit.Chem import Mol
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.MolStandardize.rdMolStandardize import MetalDisconnector

__all__ = ["Standardizer"]

from errors import ChemError
from models import Item


class Step(Enum):
    STEREO = auto()
    OPTICAL = auto()
    ISOTOPE = auto()
    CHARGE = auto()
    TAUTOMER = auto()
    METALLIC = auto()
    SANITIZE = auto()


class StandardizationError(ChemError):
    def __init__(self: Self, identifier: str, errors: dict[Step, Exception]) -> None:
        self.identifier = identifier
        self.errors = errors

    def __str__(self: Self) -> str:
        return

    def __repr__(self) -> str:
        return f"""
        Failed on {self.identifier} with errors:
        ===============================================================
        {"\n".join([str(e) for e in self.errors])}
        ===============================================================
        """


@dataclass(slots=True)
class Standardizer:
    metal_disconnector: ClassVar[MetalDisconnector] = MetalDisconnector()
    identifier: str
    mol: Mol
    errors: dict[Step, Exception] = field(init=False)

    def __post_init__(self: Self) -> None:
        self.errors = {}
        self.mol = copy(self.mol)

    def __call__(
        self: Self,
        stereo: bool = True,
        optical: bool = True,
        isotope: bool = True,
        charge: bool = False,
        tautomer: bool = False,
        metallic: bool = False,
    ) -> None:
        if not stereo:
            optical = False
            tautomer = False
        self._step(Step.METALLIC, metallic, self.metal_disconnector.Disconnect)
        self._step(Step.ISOTOPE, isotope, partial(rdMolStandardize.IsotopeParent, skipStandardize=True))
        self._step(Step.TAUTOMER, tautomer, partial(rdMolStandardize.TautomerParent, skipStandardize=True))
        self._step(Step.STEREO, stereo, partial(rdMolStandardize.StereoParent, skipStandardize=True))
        self._step(Step.CHARGE, charge, partial(rdMolStandardize.ChargeParent, skipStandardize=True))
        self._step(Step.SANITIZE, True, Chem.SanitizeMol)

    def throw(self: Self) -> None:
        if len(self.errors) > 0:
            raise StandardizationError(self.identifier, self.errors)

    def _step(self: Self, step: Step, predicate: bool, fn: Callable[[Mol], Mol]) -> bool:
        if predicate:
            try:
                self.mol = fn(self.mol)
                return True
            except Exception as e:
                self.errors[step] = e
                return False
