"""Drug crosswalk loader.

Reads ``data/static/drug_xref.parquet`` (built by ``scripts/build/build_drug_xref.py``)
and exposes lookup helpers. Verifies sha256 against ``data/static/manifest.json`` on
load and refuses to return if the hash does not match -- guards against accidental
asset drift across machines.

Canonical drug identifier is the PubChem CID (nullable ``Int64``).
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable
from pathlib import Path

import pandas as pd

XREF_REL_PATH = "drug_xref.parquet"
MANIFEST_REL_PATH = "manifest.json"


class DrugXrefManifestError(RuntimeError):
    """Raised when the on-disk parquet sha256 does not match the manifest record."""


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_drug_xref(static_dir: Path) -> pd.DataFrame:
    """Read drug_xref.parquet from ``static_dir`` after verifying its sha256.

    ``static_dir`` is the directory containing both ``drug_xref.parquet`` and
    ``manifest.json``. Raises ``FileNotFoundError`` if either is missing and
    ``DrugXrefManifestError`` on a sha256 mismatch.
    """
    parquet_path = static_dir / XREF_REL_PATH
    manifest_path = static_dir / MANIFEST_REL_PATH
    if not parquet_path.exists():
        raise FileNotFoundError(f"drug xref parquet missing: {parquet_path}")
    if not manifest_path.exists():
        raise FileNotFoundError(f"static manifest missing: {manifest_path}")
    manifest = json.loads(manifest_path.read_text())
    recorded = manifest.get("files", {}).get(XREF_REL_PATH, {}).get("sha256")
    if not recorded:
        raise DrugXrefManifestError(f"no sha256 recorded for {XREF_REL_PATH} in {manifest_path}")
    actual = _sha256(parquet_path)
    if recorded != actual:
        raise DrugXrefManifestError(
            f"{XREF_REL_PATH} sha256 mismatch -- expected {recorded}, got {actual}"
        )
    return pd.read_parquet(parquet_path)


def resolve_cid(xref: pd.DataFrame, name: str, source: str) -> int | None:
    """Return the canonical PubChem CID for ``(name, source)``, or None if unresolved.

    Lookup is case-insensitive on ``name``. ``source`` is matched exactly
    (typically ``"gdsc2"`` or ``"soragni"``).
    """
    key = name.lower()
    hit = xref[(xref["source"] == source) & (xref["input_name"].str.lower() == key)]
    if hit.empty:
        return None
    cid = hit.iloc[0]["pubchem_cid"]
    return None if pd.isna(cid) else int(cid)


def canonical_cids(xref: pd.DataFrame, names: Iterable[str], source: str) -> list[int | None]:
    """Vectorized ``resolve_cid``: return one CID (or None) per input name."""
    lower_names = pd.Series(list(names)).str.lower()
    sub = xref[xref["source"] == source].copy()
    sub["_key"] = sub["input_name"].str.lower()
    mapping = dict(zip(sub["_key"], sub["pubchem_cid"], strict=False))
    return [None if pd.isna(c := mapping.get(k)) else int(c) for k in lower_names]


def overlap_report(xref: pd.DataFrame) -> pd.DataFrame:
    """Return one row per PubChem CID present in BOTH sources.

    Columns: ``pubchem_cid``, ``gdsc2_names`` (comma-joined), ``soragni_names``
    (comma-joined), ``inchikey``, ``drugbank_id``. Useful for the panel-overlap
    discussion in ``docs/datasets.md``.
    """
    resolved = xref.dropna(subset=["pubchem_cid"])
    by_cid = resolved.groupby("pubchem_cid")["source"].nunique()
    both = by_cid[by_cid > 1].index
    sub = resolved[resolved["pubchem_cid"].isin(both)]
    out = (
        sub.groupby(["pubchem_cid", "source"])["input_name"]
        .apply(lambda s: ", ".join(sorted(set(s))))
        .unstack("source")
        .reset_index()
        .rename(columns={"gdsc2": "gdsc2_names", "soragni": "soragni_names"})
    )
    # Attach the shared inchikey + drugbank (first non-null per CID)
    extras = (
        sub.groupby("pubchem_cid")[["inchikey", "drugbank_id"]]
        .agg(lambda s: next((v for v in s if pd.notna(v)), None))
        .reset_index()
    )
    return out.merge(extras, on="pubchem_cid", how="left")
