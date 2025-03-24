
from dataclasses import dataclass
from pathlib import Path, PurePath
from typing import Any, Self
from pydantic import BaseModel
from pydantic_core import from_json, to_json

__all__ = [
    "JsonPrimitive",
    "JsonArray",
    "JsonLeaf",
    "JsonObject",
    "JsonLimb",
    "Json",
]

type JsonPrimitive = str | int | float | bool  # Disallow None!
type JsonArray = list[Json]
type JsonLeaf = JsonPrimitive | JsonArray
type JsonObject = dict[str, Json]
type JsonLimb = dict[str, JsonLeaf]
type Json = JsonLeaf | JsonObject
type SupportsJsonIo = JsonObject | JsonArray | dataclass | BaseModel


class _JsonUtils:

    @classmethod
    def as_bytes(cls: type[Self], py: JsonObject | dataclass | BaseModel) -> bytes:
        return to_json(py, inf_nan_mode="strings", bytes_mode="base64")

    @classmethod
    def as_str(cls: type[Self], py: Any | dataclass | BaseModel) -> str:
        return cls.as_bytes(py).decode(encoding="utf-8")

    @classmethod
    def read(cls: type[Self], json: str | bytes | bytearray) -> JsonObject:
        return from_json(json, allow_inf_nan=False, cache_strings="none")

    @classmethod
    def read_file(cls: type[Self], path: PurePath | str) -> JsonObject:
        return cls.read(Path(path).read_bytes())

    @classmethod
    def write_file(cls: type[Self], py: SupportsJsonIo, path: PurePath | str) -> None:
        Path(path).write_bytes(cls.as_bytes(py))
