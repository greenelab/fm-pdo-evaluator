"""Model-output model: ``Prediction``.

A single prediction record (one sample x one drug x one model x one split).
The aggregated ``PredictionRecord`` produced by the registry wraps
many of these.
"""

from __future__ import annotations

from datetime import datetime

from pydantic import BaseModel, ConfigDict, Field


class Prediction(BaseModel):
    """One predicted drug response for one sample by one model on one split."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    prediction_id: str = Field(min_length=1)
    model_version: str = Field(min_length=1)
    tranche_id: str = Field(min_length=1)
    sample_id: str = Field(min_length=1)
    drug_id: str = Field(min_length=1)
    split_name: str = Field(min_length=1)
    predicted_value: float
    predicted_responder: bool | None = None
    created_at: datetime
