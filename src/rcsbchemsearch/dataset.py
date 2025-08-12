# SPDX-FileCopyrightText: Copyright 2024, Contributors to rcsb-chem
# SPDX-PackageHomePage: https://github.com/rcsb/rcsb-chem
# SPDX-License-Identifier: BSD-3-Clause

"""
Dataset.
"""

import re
from collections import defaultdict
from functools import cached_property
from typing import Self, Callable, Iterator

from pydantic import BaseModel, ConfigDict, computed_field
from rdkit import Chem
from rdkit.Chem import Mol, rdMolDescriptors
from errors import ChemInputError
from models import ExactFlag

_INCHIKEY_PATTERN = re.compile(r"^[A-Z0-9]{14}-[A-Z0-9]{10}-[SN]$")


class Wrapped(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid", strict=True)
    mol: Mol
    inchikey: str
    inchi: str
    smiles: str

    @computed_field
    @cached_property
    def formula(self: Self) -> str:
        return rdMolDescriptors.CalcMolFormula(self.mol)

    @cached_property
    def element_counts(self: Self) -> dict[str, int]:
        counts = defaultdict(lambda: 0)
        for atom in self.mol.GetAtoms():
            counts[atom.GetSymbol()] += 1
        return counts


class Dataset(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid", strict=True)
    by_inchikey: dict[str, Wrapped]

    def test(self: Self) -> bool:
        return True  # TODO

    def convert(self: Self, chem: str) -> Wrapped:
        received = self._convert(chem)
        if isinstance(received, Wrapped):
            return received
        return self._wrap(received)

    def _wrap(self: Self, mol: Mol) -> Wrapped:
        inchi, aux = Chem.MolToInchiAndAuxInfo(mol, treatWarningAsError=False)
        return Wrapped(
            mol=mol,
            inchikey=Chem.InchiToInchiKey(inchi),
            inchi=inchi,
            smiles=Chem.MolToSmiles(mol),
        )

    def _convert(self: Self, chem: str) -> Mol | Wrapped:
        is_inchi = chem.startswith("InChI=")
        is_possibly_inchikey = _INCHIKEY_PATTERN.fullmatch(chem)
        if is_possibly_inchikey:
            try:
                return self.by_inchikey[chem]
            except KeyError:
                pass  # try again
            try:
                return self.by_inchikey[chem[:-1] + "S"]
            except KeyError:
                pass  # we'll handle
        if is_inchi:
            try:
                return Chem.MolFromInchi(chem, treatWarningAsError=False)
            except Exception:
                raise ChemInputError(chem, detail=f"The InChI is invalid")
        try:
            return Chem.MolFromSmiles(chem)
        except Exception:
            if is_possibly_inchikey:
                raise ChemInputError(chem, detail=f"The InChI Key is not in the database")
            raise ChemInputError(chem, detail=f"The SMILES is invalid")


class DatasetLoader:

    def __init__(self: Self, loader: Callable[[], Dataset]) -> None:
        self._loader = loader
        self._percent_progress = 0
        self._dataset: Dataset | None = None

    async def load(self: Self) -> Dataset:
        if self._dataset:
            return self._dataset
        self._dataset = self._loader()
        self._percent_progress = 100
        return self._dataset

    @property
    def is_ready(self: Self) -> bool:
        return bool(self._dataset)


class DatasetFactory:
    async def __call__(self: Self) -> Iterator[Wrapped]:
        pass


class DatasetSource:
    async def __call__(self: Self, flags: ExactFlag) -> DatasetFactory:
        pass
