# SPDX-FileCopyrightText: Copyright 2020-2025, Contributors to Tyrannosaurus
# SPDX-PackageHomePage: https://github.com/dmyersturnbull/tyrannosaurus
# SPDX-License-Identifier: Apache-2.0
#
# SPDX-FileCopyrightText: Copyright 2024, Contributors to rcsb-chem-search
# SPDX-PackageHomePage: https://github.com/rcsb/rcsb-chem-search
# SPDX-License-Identifier: BSD-3-Clause

# Modified from Tyrannosaurus <https://github.com/dmyersturnbull/tyrannosaurus>.

"""
Project metadata for use inside and outside `rcsb-chem-search`.
"""
from collections.abc import Sequence
from typing import overload, NotRequired, ReadOnly, TypedDict


class _FrozenList[T](Sequence[T]):

    def __init__(self, *items: T) -> None:
        self.__items = list(items)

    @overload
    def __getitem__(self, index: int) -> T:
        ...

    @overload
    def __getitem__(self, index: slice) -> Sequence[T]:
        ...

    def __getitem__(self, index: int | slice) -> T | Sequence[T]:
        return self.__items[index]

    def __len__(self):
        return len(self.__items)

    def __hash__(self) -> int:
        return hash(tuple(self.__items))

    def __eq__(self, other: Sequence[T]) -> bool:
        if isinstance(other, _FrozenList):
            return self.__items == other.__items
        if isinstance(other, Sequence):
            return self.__items == list(other)
        return NotImplemented


class UrlDict(TypedDict):
    """URLs for this project, per [PyPi Project Metadata](https://docs.pypi.org/project_metadata/)."""
    homepage: ReadOnly[str]
    changelog: ReadOnly[str]
    source: ReadOnly[str]
    documentation: NotRequired[ReadOnly[str]]
    download: NotRequired[ReadOnly[str]]
    bug: NotRequired[ReadOnly[str]]
    funding: NotRequired[ReadOnly[str]]


class About(TypedDict):
    """
    Metadata about this package.

    Attributes:
        namespace: name of the directory containing this module.
        name: pyproject `project.name`                   / importlib `Name`    .
        version: pyproject `project.version`             / importlib `version`.
        summary: pyproject `project.description`         / importlib `Summary`.
        license: pyproject `project.license.text`        / importlib `License`;
          an SPDX ID such as `Apache-2.0`.
        authors: pyproject `project.authors[*].name`      / importlib `Author` split by ` and `;
          e.g. `["Kerri Kerrigan", "Adam Addison"]`.
        maintainers: pyproject `project.authors[*].name` / importlib `Author` split by ` and `;
          e.g. `["Kerri Kerrigan", "Adam Addison"]`.
        keywords: pyproject `project.keywords`           / importlib `Keywords`.
        urls: pyproject  `project.urls` subset           / `Project-URL` subset.
          Only recognized general URLs from
          [PyPi Project Metadata](https://docs.pypi.org/project_metadata/).
          Attributes are lowercased versions of the "Name", from the "General URL" table.
    """

    name: ReadOnly[str]
    namespace: ReadOnly[str]
    version: ReadOnly[str]
    summary: ReadOnly[str]
    license: ReadOnly[str]
    authors: ReadOnly[Sequence[str]]
    maintainers: ReadOnly[Sequence[str]]
    keywords: ReadOnly[Sequence[str]]
    urls: ReadOnly[UrlDict]


__about__ = About(
    namespace="rcsbchemsearch",
    # :tyranno: name="${project.name}",
    name="rcsb-chem-search",
    # :tyranno: version="${project.version}",
    version="0.0.1-alpha0",
    # :tyranno: summary="${project.summary}",
    summary="Backend REST API supporting chemical similarity searches on RCSB PDB chemical components.",
    # :tyranno: authors=${project.authors[*].name|unpack(@)},
    authors=["Douglas Myers-Turnbull"],
    # :tyranno: maintainers=_FrozenList(${project.maintainers[*].name|unpack(@)}),
    maintainers=_FrozenList(["Douglas Myers-Turnbull"]),
    # :tyranno: name=${project.keywords|unpack(@)},
    keywords=_FrozenList("cheminformatics", "search", "rcsb", "pdb"),
    # :tyranno: license="${project.license.text}",
    license="Apache-2.0",
    urls=UrlDict(
        # :tyranno: homepage="${project.urls.Homepage}",
        homepage="https://github.com/rcsb/rcsb-chem-search",
        # :tyranno: source="${project.urls.Source}",
        source="https://github.com/rcsb/rcsb-chem-search",
        # :tyranno: documentation="${project.urls.Documentation}",
        documentation="https://github.com/rcsb/rcsb-chem-search",
        # :tyranno: changelog="${project.urls.Release Notes}",
        changelog="https://github.com/rcsb/rcsb-chem-search/releases",
        # :tyranno: download="${project.urls.Download}",
        download="https://pypi.org/rcsb/rcsb-chem-search/",
        # :tyranno: bug="${project.urls.Tracker}",
        bug="https://github.com/rcsb/rcsb-chem-search/issues",
    ),
)
