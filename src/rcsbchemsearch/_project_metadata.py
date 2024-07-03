# SPDX-FileCopyrightText: Copyright 2020-2024, Contributors to CICD
# SPDX-PackageHomePage: https://github.com/dmyersturnbull/cicd
# SPDX-License-Identifier: Apache-2.0
#
# SPDX-FileCopyrightText: Copyright 2024, the RCSB PDB and contributors
# SPDX-PackageHomePage: https://github.com/rcsb/rcsb-chem-search
# SPDX-License-Identifier: BSD-3-Clause
#
# Original source code is from:
# CICD (https://github.com/dmyersturnbull/cicd).
# This file includes modifications.
#
"""
Metadata and environment variables.
"""

import logging
import tomllib
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError
from importlib.metadata import metadata as __load
from pathlib import Path

from semver import Version

__all__ = ["ProjectMetadata"]

_pkg = Path(__file__).parent.name
logger = logging.getLogger(_pkg)
_metadata = {}
try:
    _metadata = __load(_pkg)
except PackageNotFoundError:  # nocov
    _pyproject = Path(__file__).parent / "pyproject.toml"
    if _pyproject.exists():
        _data = tomllib.loads(_pyproject.read_text(encoding="utf-8"))
        _metadata = {k.capitalize(): v for k, v in _data["project"]}
    else:
        logger.error(f"Could not load metadata for package {_pkg}. Is it installed?")


@dataclass(frozen=True, slots=True)
class _Metadata:
    pkg: str
    homepage: str
    title: str
    summary: str
    license: str
    version: Version


ProjectMetadata = _Metadata(
    pkg=_pkg,
    homepage=_metadata.get("Home-page"),
    title=_metadata.get("Name"),
    summary=_metadata.get("Summary"),
    license=_metadata.get("License"),
    version=_metadata.get("Version"),
)
