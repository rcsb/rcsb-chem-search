"""
Exceptions for rcsbchemsearch.
"""

from __future__ import annotations

import functools
import traceback
from collections.abc import Callable as Fn
from collections.abc import Generator
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from pathlib import Path
from traceback import StackSummary
from types import ModuleType
from typing import Concatenate as Ct

from rcsbchemsearch.core import StructureData
from rcsbchemsearch.core.json_utils import JsonObject

__all__ = ["AnyChemError", "MoleculeBuildError", "MoleculeOperationError"]


def decorator[**P, R](fn: Fn[P, R]) -> Fn[P, R]:
    @functools.wraps(fn)
    def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
        return fn(*args, **kwargs)

    wrapper.__is_decorator__ = True
    return wrapper


@dataclass(frozen=True, slots=True)
class AnyAppError(Exception):
    """An error from this app."""

    def properties(self) -> JsonObject:
        # noinspection PyTypeChecker
        return asdict(self)


@dataclass(frozen=True, slots=True)
class AnyChemError(AnyAppError):
    """An error relating to rdkit or chemistry."""


@dataclass(frozen=True, slots=True)
class MoleculeOperationError(AnyChemError):
    """An error caused by failure to perform some operation on an already-built molecule."""

    molecule: StructureData
    op_name: str = ""
    detail: str = ""

    def __str__(self) -> str:
        op_name = self.op_name or "process"
        detail = f": {self.detail}" if self.detail else ""
        return f"Failed to {op_name} molecule {self.molecule}{detail}."


@dataclass(frozen=True, slots=True)
class MoleculeBuildError(AnyChemError):
    """An error building a molecule from SMILES or InChI."""

    input_str: str
    detail: str = ""

    def __str__(self) -> str:
        detail = f": {self.detail}" if self.detail else ""
        return f"Could not build molecule from '{self.input_str}'{detail}."


class _ErrorUtils:
    """Utilities for handling rdkit errors."""

    def is_rdkit_error(self, exception: Exception) -> bool:
        """
        Determines whether an exception is an RDKit error.

        See:
            - [is_package_in_traceback][]
        """
        return self.is_package_in_traceback(exception, "rdkit")

    def is_package_in_traceback(self, exception: Exception, name: str) -> bool:
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

    @contextmanager
    def rdkit_op(self, molecule: StructureData, *, name: str = "") -> Generator[None]:
        """
        Wraps a Python block to catch any RDKit error and wrap it in a `MoleculeOperationError`.

        Example:
            ```python
            with ErrorUtils.rdkit_op(molecule, name="calculate solubility"):
                process_molecule(molecule)
            ```
        """
        try:
            yield
        except Exception as e:
            if self.is_rdkit_error(e):
                raise MoleculeOperationError(molecule, name, str(e))
            raise e

    @decorator
    def wrap_rdkit_op[**P, R](self, *, name: str = "") -> Fn[Fn[Ct[StructureData, P], R], Fn[Ct[StructureData, P], R]]:
        """
        Wraps a function to catch and wrap any RDKit error into an MoleculeOperationError.

        Example:
            ```python
            @ErrorUtils.wrap_rdkit_op(name="calculate solubility")
            def process_molecule(molecule: Molecule) -> Molecule:
            ```
        """
        # Python doesn't properly support higher-kinded types,
        # and -- although ParamSpec has a `bound` param -- its semantics are undefined.
        # `bound` is also only supported with the legacy `ParamSpec(...)` syntax.
        # So, we can't force `P` to start with or otherwise contain a `Molecule`.
        # But we can use `Concatenate` to prepend a `Molecule` param -- although slightly stupid, it works.
        # There's another silly Python problem:
        # `def wrapper(fn: PMR) -> PMR` can't see `PMR` unless we declare the latter `nonlocal` or `global`.
        # But `nonlocal X` requires that `X` is defined in the outer function.
        # (Even declaring `X` in an enclosing class doesn't work.)
        # And `global X` requires that `X` is defined at the module level, which wouldn't work.
        # There's no way around this, so we can't use `PMR` in the outer signature.
        # Sometimes you need to teach Python a lesson / bend it until it thinks it's Scala.
        # This isn't one of those times.
        type PMR = Fn[Ct[StructureData, P], R]

        def wrapper(fn: PMR) -> PMR:
            real_op_name: str = name or getattr(fn, "__name__", "")

            @functools.wraps(fn)
            def wrapped(molecule: StructureData, *args: P.args, **kwargs: P.kwargs) -> R:
                try:
                    return fn(molecule, *args, **kwargs)
                except Exception as e:
                    if self.is_rdkit_error(e):
                        raise MoleculeOperationError(molecule, real_op_name)
                    raise e

            return wrapped

        return wrapper


ErrorUtils = _ErrorUtils()


def inner():
    raise AnyAppError()


def outer():
    try:
        inner()
    except AnyAppError:
        raise AnyAppError()


import sysconfig


def run():
    pkgs_dir = Path(sysconfig.get_paths()["purelib"]).resolve(strict=True)
    print(pkgs_dir)
    import importlib.util
    import inspect

    spec = importlib.util.find_spec("rdkit")
    pkg_dir = Path(spec.origin).resolve(strict=True).parent.relative_to(pkgs_dir)
    pkg_mod = importlib.util.module_from_spec(spec)
    print(spec)
    print(pkg_dir)
    try:
        inner()
    except Exception as exception:
        mod: ModuleType = inspect.getmodule(exception)
        summary: StackSummary = traceback.extract_tb(exception.__traceback__)
        for frame in summary:
            path = Path(frame.filename).resolve(strict=True)
            if path.is_relative_to(pkg_dir):
                print(f"Is <{path}> in <{pkg_dir}>?")
            # print(frame.filename, frame)


run()
