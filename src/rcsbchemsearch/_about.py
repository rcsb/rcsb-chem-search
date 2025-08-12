# SPDX-FileCopyrightText: Copyright 2020-2024, Contributors to Tyrannosaurus
# SPDX-PackageHomePage: https://github.com/dmyersturnbull/tyrannosaurus
# SPDX-License-Identifier: Apache-2.0
#
# SPDX-FileCopyrightText: Copyright 2024, Contributors to rcsb-chem-search
# SPDX-PackageHomePage: https://github.com/rcsb/rcsb-chem-search
# SPDX-License-Identifier: BSD-3-Clause
#
# Adapted from Tyrannosaurus <https://github.com/dmyersturnbull/tyrannosaurus>.
"""
A set of metadata about this package.
The metadata is determined dynamically.
In order:

1. Attempts to load `importlib.metadata(pkg)`, where `pkg` is this directory name.
2. Looks for a `../../pyproject.toml` file.
3. Uses an empty map.
"""

import logging
import tomllib
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, PackageMetadata
from importlib.metadata import metadata as __importlib_load
from pathlib import Path

__all__ = ["About", "about"]

from typing import Any, Self

__pkg = Path(__file__).parent.name
logger = logging.getLogger(__pkg)


class _FakeDict(dict):

    def get_all(self: Self, key: str, fallback: Any = None) -> list[Any]:
        return self.get(key, fallback)


def _load_metadata(pkg: str) -> dict[str, str] | PackageMetadata:
    try:
        return __importlib_load(pkg)
    except PackageNotFoundError:  # nocov
        pass
    logger.debug(f"Did not find importlib metadata for package `{pkg}`. Is it installed?")
    _pyproject = Path(__file__).parent.parent / "pyproject.toml"
    _data: dict[str, str] = {}
    if _pyproject.exists():
        try:
            _data = tomllib.loads(_pyproject.read_text(encoding="utf-8"))
            logger.debug(f"Using metadata for package `{pkg}` from pyproject.toml.")
        except tomllib.TOMLDecodeError as e:
            logger.debug(f"Encountered error while decoding `{_pyproject.resolve()}`.", e)
    if len(_data) == 0:
        logger.debug(f"Did not find metadata for package `{pkg}` in pyproject.toml.")
    return _FakeDict({k.capitalize(): v for k, v in _data.get("project", {}) if isinstance(k, str)})


@dataclass(frozen=True, slots=True)
class About:
    namespace: str
    homepage: str
    title: str
    summary: str
    description: str
    keywords: list[str]
    author: str
    maintainer: str
    hyperlinks: dict[str, str]
    license: str
    version: str


__metadata = _load_metadata(__pkg)
about = About(
    namespace=__pkg,
    homepage=__metadata.get("Home-page"),
    title=__metadata.get("Name"),
    summary=__metadata.get("Summary"),
    description=__metadata.get("Description", ""),
    keywords=__metadata.get("Keywords", "").split(","),
    author=__metadata.get("Author", ""),
    maintainer=__metadata.get("Maintainer", ""),
    hyperlinks={s[:s.index(",")]: s[s.index(",")+2:] for s in __metadata.get_all("Project-URL", [])},
    license=__metadata.get("License"),
    version=__metadata.get("Version"),
)
