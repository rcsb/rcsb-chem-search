from dataclasses import dataclass
from os import PathLike
from pathlib import Path

from pydantic import BaseModel
from pydantic_core import from_json, to_json

__all__ = [
    "Json",
    "JsonArray",
    "JsonLeaf",
    "JsonLimb",
    "JsonObject",
    "JsonPrimitive",
]

type JsonPrimitive = str | int | float | bool  # Disallow None!
type JsonArray = list[Json]
type JsonObject = dict[str, Json]
type JsonLeaf = JsonPrimitive | JsonArray
type JsonLimb = dict[str, JsonLeaf]
type Json = JsonLeaf | JsonObject
type SupportsJsonIo = JsonObject | JsonArray | dataclass | BaseModel


class _JsonUtils:
    """Utilities to read and write JSON via Pydantic."""

    @classmethod
    def as_bytes(cls, py: JsonObject | dataclass | BaseModel) -> bytes:
        return to_json(py, inf_nan_mode="strings", bytes_mode="base64")

    @classmethod
    def as_str(cls, py: JsonObject | dataclass | BaseModel) -> str:
        return cls.as_bytes(py).decode(encoding="utf-8")

    @classmethod
    def read(cls, json: str | bytes | bytearray) -> JsonObject:
        return from_json(json, allow_inf_nan=False, cache_strings="none")

    @classmethod
    def read_file(cls, path: PathLike | str) -> JsonObject:
        return cls.read(Path(path).read_bytes())

    @classmethod
    def write_file(cls, py: SupportsJsonIo, path: PathLike | str) -> None:
        # noinspection PyTypeChecker
        Path(path).write_bytes(cls.as_bytes(py))
