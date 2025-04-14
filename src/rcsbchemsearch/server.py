# SPDX-FileCopyrightText: Copyright 2024, Contributors to rcsb-chem
# SPDX-PackageHomePage: https://github.com/rcsb/rcsb-chem
# SPDX-License-Identifier: BSD-3-Clause

"""
FastAPI code.
"""

import os
import asyncio
import re
from collections import defaultdict
from functools import cached_property
from typing import Self, Annotated

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.httpsredirect import HTTPSRedirectMiddleware
from fastapi.middleware.trustedhost import TrustedHostMiddleware
from pydantic import BaseModel, ConfigDict, computed_field
from rdkit import Chem
from rdkit.Chem import Mol, rdMolDescriptors
from rfc9457 import StatusProblem
from starlette.middleware.errors import ServerErrorMiddleware
from starlette.middleware.exceptions import ExceptionMiddleware
from starlette_compress import CompressMiddleware
from fastapi_problem.error import (
    BadRequestProblem,
    Problem,
    ServerProblem,
    UnprocessableProblem,
)

from rcsbchemsearch.core.errors import MoleculeBuildError
from rcsbchem._about import about
from rcsbchemsearch.models import SearchType, Metric, Formula, Exact, Item, ExactFlag
from rcsbchemsearch.core.utils import JsonObject, Json

__all__ = ["app", "search"]

if os.name == "nt":
    # workaround for asyncio loop policy for Windows platforms
    asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())

app = FastAPI(
    title=about.title,
    summary=about.summary,
    description=about.description,
    version=about.version,
    license_info={"identifier": about.license},
)
app.add_middleware(HTTPSRedirectMiddleware)
app.add_middleware(ServerErrorMiddleware)
app.add_middleware(ExceptionMiddleware)
app.add_middleware(CompressMiddleware, minimum_size=1024)
#app.add_middleware(TrustedHostMiddleware, allowed_hosts=["rcsb.org"])


dataset = Dataset(by_inchikey={})
IS_READY = False


class UnhealthyProblem(StatusProblem):
    title: str = "Not healthy"


class NotReadyProblem(StatusProblem):
    title: str = "Not ready"


@app.get("/about")
def about() -> JsonObject:
    return about.as_dict


@app.get("/schema")
def schema() -> JsonObject:
    return app.openapi_schema


@app.get("/healthy", responses={204: {"model": None}})
def healthy() -> None:
    if not dataset.is_ok():
        raise UnhealthyProblem()


@app.get("/ready", responses={204: {"model": None}})
def ready() -> JsonObject:
    return IS_READY


@app.get("/ligand")
def get(v: Annotated[str, Query(min_length=1)]) -> JsonObject:
    return dict(dataset.convert(v))


@app.get("/ligands")
def get(v: Annotated[list[str], Query(min_length=1)]) -> JsonObject:
    return {"results": [dict(dataset.convert(s)) for s in v]}


@app.get("/search")
def search(
    match_type: SearchType = SearchType.FINGERPRINT_PUBCHEM,
    structure: Annotated[str, Query()],
    similarity_metric: Metric = Metric.DICE,
    formula: Formula | None = None,
    smarts: str | None = None,
    exact_isotope: bool = False,
    exact_stereoisomer: bool = False,
    exact_tautomer: bool = False,
    limit: int | None = None,
    offset: int | None = None,
) -> None:
    """
    Search for related structures.
    """
    raise ServerProblem(
        status_code=500,
        title="Internal Server Error.",
        detail="Nothing yet.",
    )
