# SPDX-FileCopyrightText: Copyright 2024, Contributors to rcsb-chem-search
# SPDX-PackageHomePage: https://github.com/rcsb/rcsb-chem-search
# SPDX-License-Identifier: BSD-3-Clause
"""
Logging config.
"""

import functools
import traceback
from collections.abc import Callable
from contextlib import contextmanager
from dataclasses import dataclass, field
from enum import Enum
from typing import Self, Generator
from traceback import StackSummary

import rdkit.rdBase as rdBase
from rdkit import Chem, RDLogger
from rdkit import Chem
from rdkit.Chem import Mol
# from rdkit.rdBase import BlockLogs

from rcsbchemsearch.errors import MoleculeOpError

__all__ = [
    "RdkitLogLevel",
    "RdkitLogging",
    "RdkitErrors",
]
#rdBase.LogToPythonStderr()


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
        level_list = rdBase.LogStatus().split("\n") + ["rdApp.critical:enabled"]
        rd_logger.setLevel(level.value)
        init_level = [lev.split(":")[0] for lev in level_list if "enabled" in lev.lower()][0]
        yield
        rd_logger.setLevel(getattr(RDLogger, init_level.split(".")[1].upper()))  # Restore the level


@dataclass(slots=True)
class _ErrorUtils:

    @contextmanager
    def wrapping[**P, T](self, *, op: str) -> Generator[Callable[[Callable[P, T]], Callable[P, T]]]:
        """
        Wraps a function to catch and re-raise RDKit errors as MoleculeOpErrors.

        Example:
            ```python
            with RdkitErrors.wrapping(op_name="calculate solubility"):
                process_molecule(molecule)
            ```
        """
        yield self.wrap(op=op)

    def wrap[**P, T](self, *, op: str = "process") -> Callable[[Callable[P, T]], Callable[P, T]]:
        """
        Wraps a function to catch and re-raise RDKit errors as MoleculeOpErrors.

        Example:
            ```python
            @wrap_rdkit_errors(op_name="calculate solubility")
            def process_molecule(molecule: Molecule) -> Molecule:
            ```
        """
        def wrapper(fn: Callable[P, T]) -> Callable[P, T]:

            @functools.wraps(fn)
            def wrapped(*args: P.args, **kwargs: P.kwargs) -> T:
                try:
                    return fn(*args, **kwargs)
                except Exception as e:
                    if self.is_rdkit_error(e):
                        raise MoleculeOpError(args[0], op, str(e)) from e
                    raise e

            return wrapped

        return wrapper

    def is_rdkit_error(self, exception: BaseException) -> bool:
        """
        Determines whether an exception is an RDKit error.

        See:
            - [is_package_in_traceback][]
        """
        return self.is_package_in_traceback(exception, "rdkit")

    def is_package_in_traceback(self, exception: BaseException, name: str) -> bool:
        """
        Determines whether a package is in the traceback of an exception.

        Example:
            ```python
            try:
                ...
            except Exception as e:
                if is_package_in_traceback(e, "rdkit"):
                    raise MoleculeOpError(molecule, "process", "RDKit error")
                raise e
            ```
        """
        summary: StackSummary = traceback.extract_tb(exception.__traceback__)
        return any(name in frame.filename for frame in summary)


RdkitErrors = _ErrorUtils()
RdkitLogging = _RdkitLogging()
RdkitProblems = _ProblemCapture()
