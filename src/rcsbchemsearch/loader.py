# SPDX-FileCopyrightText: Copyright 2024, Contributors to rcsb-chem-search
# SPDX-PackageHomePage: https://github.com/rcsb/rcsb-chem-search
# SPDX-License-Identifier: BSD-3-Clause

"""
Load.
"""
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Self, Any

import httpx
import json
import time
from urllib.parse import quote as encode_uri

from rcsbchemsearch.models import JsonValue, JsonArray, JsonObject, JsonPrimitive, Item

BASE_SEARCH_URI = "https://search.rcsb.org/rcsbsearch/v2/query?json="
BASE_DATA_URI = "https://data.rcsb.org/graphql?query="
logger = logging.getLogger("query")
logger.setLevel("INFO")
_JSON_DUMPS_KWARGS = dict(ensure_ascii=False, allow_nan=False, sort_keys=True, separators=(",", ":"))


@dataclass(slots=True)  # TODO: frozen=True with setattr
class Loader:
    page_length: int = 1000
    delay_sec: float = 0.5
    timeout_sec: float = 1.0
    _start: int = field(init=False, default=0)

    def next(self: Self) -> list[Item] | None:
        if self._start > 0:
            time.sleep(self.delay_sec)
        with httpx.Client(http1=False, http2=True, timeout=self.timeout_sec) as client:
            encoded_request = encode_uri(json.dumps(self._search_request(), **_JSON_DUMPS_KWARGS))
            response = client.get(BASE_SEARCH_URI + encoded_request)
        response.raise_for_status()
        data = response.json().get("result_set")
        self._start += self.page_length
        if not data:
            return None
        ids = [d["identifier"] for d in data]
        encoded_request = encode_uri(self._data_request(ids))
        with httpx.Client(http1=False, http2=True, timeout=self.timeout_sec) as client:
            response = client.get(BASE_DATA_URI + encoded_request)
        response.raise_for_status()
        data = response.json()
        return [self._convert(d) for d in data["data"]["chem_comps"]]

    def _convert(self: Self, d: JsonObject) -> Item:
        c = d["chem_comp"]
        z = d["rcsb_chem_comp_descriptor"]
        i = d["rcsb_chem_comp_info"]
        return Item(
            rcsb_id=d["rcsb_id"],
            type=c["type"],
            inchi=z["InChI"],
            smiles=z["SMILES"],
            inchikey=z["InChIKey"],
            formula=c["formula"].replace(" ", "") if c["formula"] else None,  # TODO: shouldn't happen
            charge=int(c["pdbx_formal_charge"]),
            n_atoms=int(i["atom_count"]),
            n_heavy=int(i["atom_count_heavy"]),
            n_chiral=int(i["atom_count_chiral"]),
            n_aromatic=int(i["bond_count_aromatic"]),
        )

    def _search_request(self: Self) -> JsonValue:
        return {
            "query": {
                "type": "terminal",
                "label": "text_chem",
                "service": "text_chem",
                "parameters": {
                    "attribute": "rcsb_chem_comp_descriptor.InChIKey",
                    "operator": "exists",
                    "negation": False
                }
            },
            "return_type": "mol_definition",
            "request_options": {
                "paginate": {
                    "start": self._start,
                    "rows": self.page_length,
                },
                "results_content_type": [
                    "experimental"
                ]
            }
        }

    def _data_request(self: Self, ids: list[str]) -> str:
        return r"""
        {
          chem_comps(comp_ids:<<ids>>) {
            rcsb_id
            chem_comp {
              type
              formula
              pdbx_formal_charge
            }
            rcsb_chem_comp_descriptor {
              InChI
              InChIKey
              SMILES
            }
            rcsb_chem_comp_info {
              atom_count
              atom_count_heavy
              atom_count_chiral
              bond_count_aromatic
            }
          }
        }
        """.replace("<<ids>>", "[" + ",".join([f'"{x}"' for x in ids]) + "]")


def _default(d: Any) -> JsonObject:
    if hasattr(d, "model_dump"):
        return d.model_dump()
    raise TypeError(str(type(d)))


def main():
    finder = Loader()
    data = []
    while page := finder.next():
        data += page
        logger.info(f"Got {len(data)} results")
        #zz = {d["rcsb_id"]: {x: y for x, y in d.items() if x != "rcsb_id"} for d in data}
        jj = json.dumps({"data": data}, **_JSON_DUMPS_KWARGS, default=_default)
        Path("data.json").write_text(jj, encoding="utf-8")
        print(f"Wrote {len(data)} items")
    #Path("ids.txt").write_text("\n".join(data), encoding="utf-8")


if __name__ == "__main__":
    main()

