import os
import sys
from dataclasses import dataclass
from typing import TypedDict

__all__ = ["StructureData", "StructureDict"]
StructureDict = TypedDict("StructureDict", {"key": str, "SMILES": str, "InChI": str, "InChI Key": str})
IS_STDOUT_TTY = os.isatty(sys.stdout.fileno())
IS_STDERR_TTY = os.isatty(sys.stderr.fileno())


@dataclass(slots=True, frozen=True)
class StructureData:
    key: str
    smiles: str
    inchi: str
    inchikey: str

    def __str__(self) -> str:
        return f"<{self.key}={self.inchikey}: {self.smiles}>"

    @property
    def as_dict(self) -> StructureDict:
        return StructureDict(
            **{"key": self.key, "SMILES": self.smiles, "InChI": self.inchi, "InChI Key": self.inchikey}
        )
