from dataclasses import dataclass
from typing import TypedDict

__all__ = ["Molecule", "MoleculeDict"]


class MoleculeDict(TypedDict):
    key: str
    smiles: str
    inchi: str
    inchikey: str


@dataclass(slots=True, frozen=True)
class Molecule:
    key: str
    smiles: str
    inchi: str
    inchikey: str

    def __str__(self) -> str:
        return f"<{self.key}={self.inchikey}: {self.smiles}>"

    @property
    def as_molecule_dict(self) -> MoleculeDict:
        # Due to a limitation in Python
        return MoleculeDict(
            key=self.key,
            smiles=self.smiles,
            inchi=self.inchi,
            inchikey=self.inchikey,
        )
