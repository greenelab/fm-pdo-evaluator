"""Property tests for the harness schema."""

from __future__ import annotations

import datetime as dt

import pytest
from hypothesis import given
from hypothesis import strategies as st
from pydantic import ValidationError

from fmharness.schema import (
    BaselineExpression,
    DrugAssay,
    EnvironmentSnapshot,
    LeakageProfile,
    Patient,
    Prediction,
    Sample,
    Tranche,
)

identifier = st.text(min_size=1, max_size=64).filter(lambda s: s.strip() != "")
content_hash = st.text(
    min_size=64,
    max_size=64,
    alphabet="0123456789abcdef",
)
finite_float = st.floats(allow_nan=False, allow_infinity=False)


@given(patient_id=identifier, tranche_id=identifier, tissue=identifier)
def test_patient_round_trip(patient_id: str, tranche_id: str, tissue: str) -> None:
    p = Patient(patient_id=patient_id, tranche_id=tranche_id, tissue_of_origin=tissue)
    assert Patient.model_validate_json(p.model_dump_json()) == p


@given(sample_id=identifier, patient_id=identifier, tranche_id=identifier)
def test_sample_round_trip(sample_id: str, patient_id: str, tranche_id: str) -> None:
    s = Sample(sample_id=sample_id, patient_id=patient_id, tranche_id=tranche_id)
    assert Sample.model_validate_json(s.model_dump_json()) == s


def test_patient_frozen() -> None:
    p = Patient(patient_id="p1", tranche_id="t1", tissue_of_origin="liver")
    with pytest.raises(ValidationError):
        p.subtype = "HCC"  # type: ignore[misc]


def test_patient_rejects_unknown_granularity() -> None:
    with pytest.raises(ValidationError):
        Patient(
            patient_id="p1",
            tranche_id="t1",
            tissue_of_origin="liver",
            subtype_granularity="medium",  # type: ignore[arg-type]
        )


def test_drug_assay_rejects_unknown_metric() -> None:
    with pytest.raises(ValidationError):
        DrugAssay(
            assay_id="a1",
            sample_id="s1",
            drug_id="d1",
            drug_name="lenvatinib",
            response_metric="bogus",  # type: ignore[arg-type]
            response_value=0.5,
        )


@given(value=finite_float, responder=st.booleans())
def test_drug_assay_round_trip(value: float, responder: bool) -> None:
    a = DrugAssay(
        assay_id="a1",
        sample_id="s1",
        drug_id="d1",
        drug_name="lenvatinib",
        response_metric="viability",
        response_value=value,
        responder=responder,
    )
    assert DrugAssay.model_validate_json(a.model_dump_json()) == a


def test_baseline_expression_requires_positive_gene_count() -> None:
    with pytest.raises(ValidationError):
        BaselineExpression(
            sample_id="s1",
            expression_matrix_uri="/tmp/x.h5ad",
            gene_count=0,
            normalization="median_of_ratios",
        )


def test_baseline_expression_rejects_unknown_normalization() -> None:
    with pytest.raises(ValidationError):
        BaselineExpression(
            sample_id="s1",
            expression_matrix_uri="/tmp/x.h5ad",
            gene_count=20000,
            normalization="zscore",  # type: ignore[arg-type]
        )


@given(
    tranche_id=identifier,
    source=identifier,
    version=identifier,
    h=content_hash,
)
def test_tranche_round_trip(tranche_id: str, source: str, version: str, h: str) -> None:
    t = Tranche(
        tranche_id=tranche_id,
        source=source,
        version=version,
        ingestion_date=dt.date(2026, 5, 19),
        patient_count=10,
        sample_count=12,
        drug_count=7,
        content_hash=h,
    )
    assert Tranche.model_validate_json(t.model_dump_json()) == t


def test_tranche_rejects_negative_counts() -> None:
    with pytest.raises(ValidationError):
        Tranche(
            tranche_id="t1",
            source="yang",
            version="2024.05",
            ingestion_date=dt.date(2026, 5, 19),
            patient_count=-1,
            sample_count=0,
            drug_count=0,
            content_hash="0" * 64,
        )


def test_tranche_rejects_non_sha256_hash() -> None:
    with pytest.raises(ValidationError):
        Tranche(
            tranche_id="t1",
            source="yang",
            version="2024.05",
            ingestion_date=dt.date(2026, 5, 19),
            patient_count=0,
            sample_count=0,
            drug_count=0,
            content_hash="not-a-hash",
        )


@given(value=finite_float)
def test_prediction_round_trip(value: float) -> None:
    p = Prediction(
        prediction_id="pred-1",
        model_version="tahoe_x1@v1.0",
        tranche_id="yang_pdo_liver_2024",
        sample_id="s1",
        drug_id="lenvatinib",
        split_name="leave_patient_out:fold-3",
        predicted_value=value,
        predicted_responder=None,
        created_at=dt.datetime(2026, 5, 19, 12, 0, 0),
    )
    assert Prediction.model_validate_json(p.model_dump_json()) == p


def test_leakage_profile_fraction_bounds() -> None:
    base = dict(
        tranche_id="yang_pdo_liver_2024",
        model_version="tahoe_x1@v1.0",
        generated_at=dt.datetime(2026, 5, 19, 12, 0, 0),
    )
    with pytest.raises(ValidationError):
        LeakageProfile(drug_overlap_fraction=1.5, **base)  # type: ignore[arg-type]
    with pytest.raises(ValidationError):
        LeakageProfile(drug_overlap_fraction=-0.1, **base)  # type: ignore[arg-type]


def test_environment_snapshot_round_trip() -> None:
    snap = EnvironmentSnapshot(
        code_commit="a" * 40,
        python_version="3.11.9",
        seed=42,
        cuda_deterministic=True,
        data_commit="b" * 64,
    )
    assert EnvironmentSnapshot.model_validate_json(snap.model_dump_json()) == snap


def test_environment_snapshot_rejects_short_commit() -> None:
    with pytest.raises(ValidationError):
        EnvironmentSnapshot(
            code_commit="abc",
            python_version="3.11.9",
            seed=0,
            cuda_deterministic=False,
            data_commit="b" * 64,
        )
