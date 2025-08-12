# SPDX-FileCopyrightText: Copyright 2024, Contributors to rcsb-chem-search
# SPDX-PackageHomePage: https://github.com/rcsb/rcsb-chem-search
# SPDX-License-Identifier: BSD-3-Clause

"""
Models.
"""
from enum import auto, StrEnum
from pydantic import BaseModel
from typing import Self, ClassVar
from rdkit import DataStructs
from rdkit.Chem import Mol
from rdkit.Chem import AllChem

__all__ = [
    "Item",
    "SearchMetaType",
    "SearchType",
    "Use",
    "Formula",
    "Metric",
    "ElementReq",
    "JsonPrimitive",
    "JsonArray",
    "JsonObject",
    "JsonValue",
    "Fingerprint",
    "RdkitFingerprint",
    "Morgan2Fingerprint",
    "Morgan3Fingerprint",
    "Morgan4Fingerprint",
    "TopologicalTorsionFingerprint",
]


class _ApiEnum(StrEnum):

    @staticmethod
    def _generate_next_value_(name: str, start: int, count: int, last_values: list[str]) -> str:
        return name.upper().replace("_", "-")


class SearchMetaType(StrEnum):
    STRUCTURAL = auto()
    SMARTS = auto()
    FINGERPRINT = auto()


class SearchType(_ApiEnum):
    EQUIVALENT = auto()
    SMARTS = auto()
    SUBSTRUCTURE = auto()
    SUPERSTRUCTURE = auto()
    FINGERPRINT_RDKIT = auto()
    FINGERPRINT_MORGAN2 = auto()
    FINGERPRINT_MORGAN3 = auto()
    FINGERPRINT_MORGAN4 = auto()
    FINGERPRINT_TOPOLOGICAL_TORSION = auto()
    #FINGERPRINT_AVALON = auto()
    #FINGERPRINT_E3FP = auto()
    #FINGERPRINT_2D_PHARMACOPHORE = auto()
    #FINGERPRINT_EXTENDED_REDUCED_GRAPHS = auto()

    @property
    def meta_type(self: Self) -> SearchMetaType:
        if self is SearchType.SMARTS:
            return SearchMetaType.SMARTS
        if "FINGERPRINT" in self.value:
            return SearchMetaType.FINGERPRINT
        return SearchMetaType.STRUCTURAL

    @property
    def rdkit_fingerprint_name(self: Self) -> str | None:
        if self.meta_type is SearchMetaType.FINGERPRINT:
            return self.name.lower().replace("fingerprint_", "")
        return None


class Metric(_ApiEnum):
    DICE = auto()
    HAMMING = auto()
    JACCARD = auto()
    COSINE = auto()


class Item(BaseModel):
    rcsb_id: str
    type: str
    smiles: str
    inchi: str
    inchikey: str
    formula: str
    charge: int
    n_atoms: int
    n_heavy: int
    n_chiral: int
    n_aromatic: int

    class Config:
        frozen = True
        allow_inf_nan = False

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return f"<{self.rcsb_id}>:{self.inchikey}:({self.smiles})"

    def __lt__(self, other: "Item") -> bool:
        return self.rcsb_id < other.rcsb_id


class ElementReq(BaseModel):
    element: str
    min: int | None
    max: int | None


class Formula(BaseModel):
    elements: list[ElementReq]
    allows_additional: bool


class Use(BaseModel):
    optical: bool = True
    cis_trans: bool = True
    tautomer: bool = True
    protonation: bool = True
    hydrogen: bool = True
    isotope: bool = True


class Fingerprint(BaseModel):
    type: SearchType

    def __call__(self: Self, mol: Mol) -> bytes:
        raise NotImplementedError("Subclasses should implement this method.")


class RdkitFingerprint(Fingerprint):
    def __call__(self: Self, mol: Mol) -> bytes:
        fp = AllChem.GetRDKitFPGenerator(
            minPath=1,
            maxPath=7,
            useHs=True,
            branchedPaths=True,
            useBondOrder=True,
            includeChirality=True,
            countSimulation=True,
            countBounds=[2, 4, 8, 16, 32, 64],
            fpSize=4096,
            numBitsPerFeature=2,
        ).GetFingerprint(mol)
        return DataStructs.BitVectToBytes(fp)


class _MorganFingerprint(Fingerprint):
    radius: ClassVar[int]
    def __call__(self: Self, mol: Mol) -> bytes:
        fp = AllChem.GetMorganGenerator(
            radius=self.radius,
            includeChirality=True,
            numBitsPerFeature=2,
            fpSize=4096,
        ).GetFingerprint(mol)
        return DataStructs.BitVectToBytes(fp)


class Morgan2Fingerprint(_MorganFingerprint):
    radius: ClassVar[int] = 2


class Morgan3Fingerprint(_MorganFingerprint):
    radius: ClassVar[int] = 3


class Morgan4Fingerprint(_MorganFingerprint):
    radius: ClassVar[int] = 4


class TopologicalTorsionFingerprint(Fingerprint):
    def __call__(self: Self, mol: Mol) -> bytes:
        fp = AllChem.GetTopologicalTorsionGenerator(
            torsionAtomCount=4,
            includeChirality=True,
            countSimulation=True,
            countBounds=[2, 4, 8, 16, 32, 64],
            fpSize=4096,
        ).GetFingerprint(mol)
        return DataStructs.BitVectToBytes(fp)


JsonPrimitive = str | int | float | bool  # we disallow None!
JsonArray = list["JsonValue"]
JsonObject = dict[str, "JsonValue"]
JsonValue = JsonPrimitive | JsonArray | JsonObject
