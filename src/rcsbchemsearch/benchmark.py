# SPDX-FileCopyrightText: Copyright 2024, Contributors to rcsb-chem-search
# SPDX-PackageHomePage: https://github.com/rcsb/rcsb-chem-search
# SPDX-License-Identifier: BSD-3-Clause

"""
Load.
"""
import json
import time
from dataclasses import dataclass, field
from urllib.parse import quote as encode_uri
from pathlib import Path
from typing import Self, Any
import httpx
from pydantic import BaseModel
from rcsbchemsearch.loader import _JSON_DUMPS_KWARGS, BASE_SEARCH_URI
from rcsbchemsearch.models import Item, JsonObject

_DATA_PATH = Path("~data/data.json")
_OUTPUT_FILE = Path("~data/rcsb-results.json")


def _default(d: Any) -> JsonObject:
    if hasattr(d, "model_dump"):
        return d.model_dump()
    raise TypeError(str(type(d)))


class Hit(BaseModel):
    identifier: str
    score: float

    def __str__(self: Self) -> str:
        return repr(self)

    def __repr__(self: Self) -> str:
        return f"{self.identifier}({self.score})"


@dataclass(slots=True)  # TODO: frozen=True with setattr
class Query:
    delay_sec: float = 0.0
    timeout_sec: float = 1.0
    _start: int = field(init=False, default=0)

    def next(self: Self, item: Item) -> list[Hit]:
        if self._start > 0:
            time.sleep(self.delay_sec)
        encoded_request = encode_uri(json.dumps(self._data_request(item), **_JSON_DUMPS_KWARGS))
        with httpx.Client(http1=False, http2=True, timeout=self.timeout_sec) as client:
            response = client.get(BASE_SEARCH_URI + encoded_request)
        response.raise_for_status()
        data = response.json().get("result_set")
        return [Hit.model_validate(d) for d in data]

    def _data_request(self: Self, item: Item) -> JsonObject:
        return {
            "query": {
                "type": "terminal",
                "service": "chemical",
                "parameters": {
                    "value": item.smiles,
                    "type": "descriptor",
                    "descriptor_type": "SMILES",
                    "match_type": "fingerprint-similarity"
                }
            },
            "return_type": "mol_definition",
            "request_options": {
                "paginate": {
                    "start": 0,
                    "rows": 100
                },
                "results_content_type": [
                    "experimental"
                ],
                "sort": [
                    {
                        "sort_by": "score",
                        "direction": "desc"
                    }
                ],
                #"scoring_strategy": "combined"
            }
        }


def write(all_hits: dict[str, list[Hit]]):
    _OUTPUT_FILE.write_text(json.dumps(all_hits, default=_default, **_JSON_DUMPS_KWARGS), encoding="utf-8")


def main():
    data = json.loads(_DATA_PATH.read_text(encoding="utf-8"))["data"]
    items = [
        Item.model_validate(d)
        for d in data
        if d["smiles"] and d["formula"] and d["inchikey"] and d["inchi"]
    ]
    all_hits = {}
    query = Query()
    for i, item in enumerate(items):
        all_hits[item.rcsb_id] = query.next(item)
        #print(f"Found {len(all_hits[item.rcsb_id])} hits for {item.rcsb_id}.")
        if i % 100 == 0 and i > 0:
            write(all_hits)
            print(f"Processed {i} items.")
    write(all_hits)
    print(f"Done. Finished {len(all_hits)} items.")


if __name__ == "__main__":
    main()
