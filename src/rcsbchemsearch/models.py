import enum
import re
from collections.abc import Callable
from dataclasses import dataclass
from enum import Enum

from rdkit import Chem
from rdkit.Chem import Mol

from rcsbchemsearch.core import StructureData
from rcsbchemsearch.core.errors import MoleculeBuildError


class InputType(Enum):
    INCHI = "InChI"
    INCHIKEY = "InChI Key"
    SMILES = "SMILES"
    RCSB_ID = "rcsb id"


class OpStatus(Enum):
    """
    The success, failure, or other outcome of an attempted operation.
    Note that this may not reflect the final state of the operation, or whether results are available.

    Attributes:
    - COMPLETED: The operation completed successfully.
    - HALTED: The operation was halted before completion, resulting in partial or no results.
    - HAD_ERRORS: The operation completed, but some or all results may be unavailable, incorrect, or incomplete.
    - FAILED: The operation failed to complete (e.g. due to a general error); no results are available.

    `HAD_ERRORS` pertains especially to operations that are applied to multiple items.
    It is distinct from `FAILED` in that it guarantees the operation was attempted on all items.
    Suppose an operation is applied to 100 items, and each is interrupted by an error.
    Then, even though no results are available, `HAD_ERRORS` is the correct status.
    """

    COMPLETED = enum.auto()
    HALTED = enum.auto()
    HAD_ERRORS = enum.auto()
    FAILED = enum.auto()

    @property
    def can_contain_results(self) -> bool:
        return self in (OpStatus.COMPLETED, OpStatus.HAD_ERRORS, OpStatus.HALTED)


@dataclass(slots=True, frozen=True)
class Structure(StructureData):
    mol: Mol


# @D.describe("{<key>} <{<smiles>}>, score={<score>}")
@dataclass(slots=True, frozen=True)
class Tautomer(Structure):
    score: int


@dataclass(slots=True, frozen=True, eq=False)
class MolFactory:
    def new_mol(self, input_str: str) -> Mol:
        raise NotImplementedError()


@dataclass(slots=True, frozen=True, eq=False)
class SimpleMolFactory(MolFactory):
    def new_mol(self, input_str: str) -> Mol:
        try:
            if input_str.startswith("InChI="):
                return Chem.MolFromInchi(input_str, treatWarningAsError=False)
            return Chem.MolFromSmiles(input_str)
        except Exception:
            raise MoleculeBuildError(input_str)


@dataclass(slots=True, frozen=True, eq=False)
class MoleculeFactory:
    def of(self, mol: Mol | str, key: str | None = None) -> Structure:
        raise NotImplementedError()


@dataclass(slots=True, frozen=True, eq=False)
class SimpleMoleculeFactory(MoleculeFactory):
    mol_factory: MolFactory

    def of(self, mol: Mol | str, key: str | None = None) -> Structure:
        if isinstance(mol, str):
            mol = self.mol_factory.new_mol(mol)
        return self._to_molecule(mol, None)

    def _to_molecule(self, mol: Mol, key: str | None) -> Structure:
        smiles = self._try(mol, InputType.SMILES, Chem.MolToSmiles)
        inchi = self._try(mol, InputType.INCHI, Chem.MolToInchi)
        inchikey = self._try(mol, InputType.INCHIKEY, Chem.MolToInchiKey)
        return Structure(key if key else inchikey, smiles, inchi, inchikey, mol)

    @staticmethod
    def _try(mol: Mol, op_name: InputType, op: Callable[[Mol], str]) -> str:
        try:
            return op(mol)
        except Exception:
            raise MoleculeBuildError(mol, detail=f"An error occurred while creating {op_name.value}.")


# Molecules = SimpleMoleculeFactory(SimpleMolFactory())


_INCHIKEY_PATTERN = re.compile(r"^[A-Z0-9]{14}-[A-Z0-9]{10}-[SN]$")
