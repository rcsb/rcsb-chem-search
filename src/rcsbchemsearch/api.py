# SPDX-FileCopyrightText: Copyright 2024, the RCSB PDB and contributors
# SPDX-PackageHomePage: https://github.com/rcsb/rcsb-chem-search
# SPDX-License-Identifier: BSD-3-Clause
#
"""
FastAPI code.
"""
import os

import asyncio
from fastapi import FastAPI, HTTPException
#from fastapi.responses import Response, ORJSONResponse
#from pydantic import BaseModel
#from starlette import status
from fastapi.middleware.httpsredirect import HTTPSRedirectMiddleware
from fastapi.middleware.trustedhost import TrustedHostMiddleware
from starlette_compress import CompressMiddleware

from rcsbchemsearch.models import Consider, Formula

__all__ = []

if os.name == "nt":
    # workaround for asyncio loop policy for Windows platforms
    asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())

app = FastAPI(
    title="rcsb-chem-search",
    description="rcsb-chem-search",
    version="0.0.1",
    # ...
)
app.add_middleware(ServerErrorMiddleware)
app.add_middleware(ExceptionMiddleware)
app.add_middleware(HTTPSRedirectMiddleware)
app.add_middleware(CompressMiddleware, minimum_size=1024)
app.add_middleware(TrustedHostMiddleware, allowed_hosts=["rcsb.org"])


DEFAULT_CONSIDER = Consider()


@app.get("/search/{search_type}")
def search(
    query: str,
    algorithm_params: dict | None = None,
    metric: str = "dice",
    formula: Formula | None = None,
    use: Consider = DEFAULT_CONSIDER,
    limit: int | None = None,
    offset: int | None = None,
) -> None:
    """
    Search for related structures.
    """
    raise HTTPException(
        status_code=500,
        title="Internal Server Error.",
        detail="Nothing yet.",
    )
