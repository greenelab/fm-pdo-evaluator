"""Versioned data-bundle model: ``Tranche``."""

from __future__ import annotations

import re
from datetime import date

from pydantic import BaseModel, ConfigDict, Field, field_validator

_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


class Tranche(BaseModel):
    """A versioned, content-hashed bundle of data from one source dataset."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    tranche_id: str = Field(min_length=1)
    source: str = Field(min_length=1)
    version: str = Field(min_length=1)
    ingestion_date: date
    patient_count: int = Field(ge=0)
    sample_count: int = Field(ge=0)
    drug_count: int = Field(ge=0)
    subtypes: tuple[str, ...] = ()
    content_hash: str
    description: str | None = None

    @field_validator("content_hash")
    @classmethod
    def _hash_is_sha256(cls, v: str) -> str:
        if not _SHA256_RE.match(v):
            raise ValueError("content_hash must be a 64-character lowercase hex sha256")
        return v
