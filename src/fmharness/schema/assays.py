"""Measurement models: ``DrugAssay`` and ``BaselineExpression``."""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

ResponseMetric = Literal["viability", "ic50", "auc", "responder_binary"]
NormalizationMethod = Literal[
    "raw",
    "median_of_ratios",
    "log1p_median_of_ratios",
    "tpm",
]
GeneIdNamespace = Literal["ensembl", "symbol", "entrez"]


class DrugAssay(BaseModel):
    """One drug-response measurement on one sample."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    assay_id: str = Field(min_length=1)
    sample_id: str = Field(min_length=1)
    drug_id: str = Field(min_length=1)
    drug_name: str = Field(min_length=1)
    response_metric: ResponseMetric
    response_value: float
    responder: bool | None = None
    inchikey: str | None = None
    drugbank_id: str | None = None
    pubchem_cid: int | None = Field(default=None, ge=0)


class BaselineExpression(BaseModel):
    """Reference to a per-sample baseline expression artifact (AnnData on disk).

    The actual matrix lives at ``expression_matrix_uri``; this model captures
    the metadata that consumers need to interpret it without re-deriving
    upstream choices (reference genome, annotation version, normalization).
    """

    model_config = ConfigDict(extra="forbid", frozen=True)

    sample_id: str = Field(min_length=1)
    expression_matrix_uri: str = Field(min_length=1)
    gene_count: int = Field(gt=0)
    gene_id_namespace: GeneIdNamespace = "ensembl"
    normalization: NormalizationMethod
    reference_genome: str = "GRCh38"
    reference_annotation: str | None = None
