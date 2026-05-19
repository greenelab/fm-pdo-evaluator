"""Pydantic schema models for the harness.

These models define the contract for every artifact the harness produces:
patients and samples (the data subjects), drug assays and baseline
expression (the measurements), tranches (versioned data bundles),
predictions (model outputs), and the provenance metadata (leakage profile,
environment snapshot) attached to every prediction record.

All models are immutable (``frozen=True``) and reject extra fields
(``extra="forbid"``).
"""

from fmharness.schema.assays import (
    BaselineExpression,
    DrugAssay,
    NormalizationMethod,
    ResponseMetric,
)
from fmharness.schema.entities import Patient, Sample, SubtypeGranularity
from fmharness.schema.predictions import Prediction
from fmharness.schema.provenance import EnvironmentSnapshot, LeakageProfile
from fmharness.schema.tranches import Tranche

__all__ = [
    "BaselineExpression",
    "DrugAssay",
    "EnvironmentSnapshot",
    "LeakageProfile",
    "NormalizationMethod",
    "Patient",
    "Prediction",
    "ResponseMetric",
    "Sample",
    "SubtypeGranularity",
    "Tranche",
]
