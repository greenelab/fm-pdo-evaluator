"""Data-subject models: ``Patient`` and ``Sample``."""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

SubtypeGranularity = Literal["fine", "coarse"]

MetadataValue = str | int | float | bool | None


class Patient(BaseModel):
    """A patient contributing one or more PDOs to a tranche."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    patient_id: str = Field(min_length=1)
    tranche_id: str = Field(min_length=1)
    tissue_of_origin: str = Field(min_length=1)
    subtype: str | None = None
    subtype_granularity: SubtypeGranularity | None = None
    metadata: dict[str, MetadataValue] = Field(default_factory=dict)


class Sample(BaseModel):
    """A single PDO derived from a patient, at a specific passage/replicate."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    sample_id: str = Field(min_length=1)
    patient_id: str = Field(min_length=1)
    tranche_id: str = Field(min_length=1)
    passage: int | None = Field(default=None, ge=0)
    replicate: int | None = Field(default=None, ge=0)
    metadata: dict[str, MetadataValue] = Field(default_factory=dict)
