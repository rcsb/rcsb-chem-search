import functools
import traceback
from collections.abc import Callable as Fn
from collections.abc import Generator
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from traceback import StackSummary
from typing import Concatenate as Cat

from . import Molecule
from .json_utils import JsonObject

__all__ = ["BaseChemError", "MoleculeBuildError", "MoleculeOperationError"]


def decorator[**P, R](fn: Fn[P, R]) -> Fn[P, R]:
    @functools.wraps(fn)
    def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
        return fn(*args, **kwargs)

    wrapper.__is_decorator__ = True
    return wrapper


@dataclass(frozen=True, slots=True)
class AppError(Exception):
    """An error from this app."""

    def properties(self) -> JsonObject:
        # noinspection PyTypeChecker
        return asdict(self)


@dataclass(frozen=True, slots=True)
class BaseChemError(AppError):
    """An error relating to rdkit or chemistry."""


@dataclass(frozen=True, slots=True)
class MoleculeOperationError(BaseChemError):
    """An error caused by failure to process"""

    molecule: Molecule
    op_name: str = ""
    detail: str = ""

    def __str__(self) -> str:
        op_name = self.op_name if self.op_name else "process"
        detail = f": {self.detail}" if self.detail else ""
        return f"Failed to {op_name} molecule {self.molecule}{detail}."


@dataclass(frozen=True, slots=True)
class MoleculeBuildError(BaseChemError):
    input_str: str
    detail: str = ""

    def __str__(self) -> str:
        detail = f": {self.detail}" if self.detail else ""
        return f"Could not build molecule from '{self.input_str}'{detail}."


class _ErrorUtils:
    """"""

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

    @contextmanager
    def rdkit_op(self, molecule: Molecule, *, name: str = "") -> Generator[None]:
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
    def wrap_rdkit_op[**P, R](self, *, name: str = "") -> Fn[Fn[Cat[Molecule, P], R], Fn[Cat[Molecule, P], R]]:
        """
        Wraps a function to catch and re-raise RDKit errors as MoleculeOpErrors.

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
        type PMR = Fn[Cat[Molecule, P], R]

        def wrapper(fn: PMR) -> PMR:
            real_op_name: str = name or getattr(fn, "__name__", "")

            @functools.wraps(fn)
            def wrapped(molecule: Molecule, *args: P.args, **kwargs: P.kwargs) -> R:
                try:
                    return fn(molecule, *args, **kwargs)
                except Exception as e:
                    if self.is_rdkit_error(e):
                        raise MoleculeOperationError(molecule, real_op_name)
                    raise e

            return wrapped

        return wrapper


ErrorUtils = _ErrorUtils()
