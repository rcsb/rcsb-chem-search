# SPDX-FileCopyrightText: Copyright 2024, Contributors to rcsb-chem-search
# SPDX-PackageHomePage: https://github.com/rcsb/rcsb-chem-search
# SPDX-License-Identifier: BSD-3-Clause
"""
Logging config.
"""

from collections.abc import Generator
from contextlib import contextmanager
from dataclasses import dataclass, field
from enum import Enum
from typing import Self

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import Mol

# from rdkit.rdBase import BlockLogs
# rdBase.LogToPythonStderr()

__all__ = [
    "RdkitLogLevel",
    "RdkitLogging",
]


class RdkitSeed:
    @classmethod
    def set(cls, n: int, /) -> None:
        rdBase.SeedRandomNumberGenerator(n)


@dataclass(frozen=True, slots=True)
class Problem:
    type: str
    atom_index: int
    message: str


class _ProblemCapture:
    def of(self, mol: Mol, /) -> Generator[Problem]:
        for p in Chem.DetectChemistryProblems(mol):
            yield Problem(p.GetType(), p.GetAtomIdx(), p.Message())


class RdkitLogLevel(Enum):
    """Log levels from rdkit and how they map to standard `logging`."""

    DEBUG = RDLogger.DEBUG
    INFO = RDLogger.INFO
    WARNING = RDLogger.WARNING
    ERROR = RDLogger.ERROR
    CRITICAL = RDLogger.CRITICAL

    @classmethod
    def from_rdkit(cls, level: Self | int | str) -> Self:
        """Returns a level from an rdkit-defined number or name."""
        if not isinstance(level, (cls, str)):
            msg = f"level must be a name; got '{level}' of type {type(level)}."
            raise TypeError(msg)
        if isinstance(level, cls):
            return cls
        return {m.name: m for m in cls}[level.upper()]

    @classmethod
    def from_std(cls, level: int | str) -> Self:
        """Returns a level from a name defined by standard `logging` (equiv. `loguru`)."""
        if not isinstance(level, (int, str)):
            msg = f"level must be a name; got '{level}' of type {type(level)}."
            raise TypeError(msg)
        return {m.name: m for m in cls}[level.upper()]


@dataclass(slots=True)
class _RdkitLogging:
    """Configuration for rdkit's logger."""

    level: RdkitLogLevel = field(default=RdkitLogLevel.INFO, init=False)

    def set(self, *, level: RdkitLogLevel | str) -> None:
        self.level = RdkitLogLevel.from_rdkit(level)
        RDLogger.logger().setLevel(self.level.value)

    @contextmanager
    def at_level(self, level: RdkitLogLevel | str) -> Generator[None]:
        """Modifies the rdkit level in-place for a context manager (**not** thread-safe)."""
        """Context manager for controlling the RDKit logger level."""
        rd_logger = RDLogger.logger()
        self.level = RdkitLogLevel.from_rdkit(level)
        # In case no logging is enabled, it will be set to critical.
        level_list = [*rdBase.LogStatus().split("\n"), "rdApp.critical:enabled"]
        rd_logger.setLevel(level.value)
        init_level = next(lev.split(":")[0] for lev in level_list if "enabled" in lev.lower())
        yield
        rd_logger.setLevel(getattr(RDLogger, init_level.split(".")[1].upper()))  # Restore the level


RdkitLogging = _RdkitLogging()
RdkitProblems = _ProblemCapture()
