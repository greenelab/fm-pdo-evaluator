#!/usr/bin/env python3
"""Build the drug crosswalk for the Soragni + GDSC2 sarcoma panels.

Resolves each input drug name to PubChem CID (canonical) + InChIKey + DrugBank ID,
writes one row per (source, input_name) tuple to data/static/drug_xref.parquet, and
updates data/static/manifest.json (same schema as scripts/download/_utils.py).

Resolution chain:

  1. GDSC2: use the PUBCHEM column where populated (parseable as int).
  2. GDSC2 / Soragni name not yet resolved: query PubChem PUG REST
     /compound/name/{name}/cids/JSON. Pick the smallest CID when multiple
     hits are returned (PubChem convention for the "parent" / canonical
     compound vs salt forms).
  3. Soragni Drug_Name == any GDSC2 DRUG_NAME or SYNONYM: reuse the
     GDSC2 row's CID without re-querying PubChem. Recorded as
     resolution_method = "soragni_via_gdsc2_synonym".
  4. CID -> InChIKey: PUG REST /compound/cid/{cid}/property/InChIKey/JSON.
  5. InChIKey -> DrugBank: UniChem cross-reference API
     https://www.ebi.ac.uk/unichem/rest/inchikey/{key} -> filter src_id=2.

Public APIs only; no auth. PubChem rate limit is 5 req/s; this script sleeps
~0.25 s between calls to stay well under. ~329 drugs total -> ~6-7 min.

Usage:
    python scripts/build/build_drug_xref.py
    python scripts/build/build_drug_xref.py --refresh  # rebuild from scratch

The crosswalk is committed to the repo. Re-run only when the input drug lists
change (new GDSC2 release, Soragni drug-screen update).
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path
from urllib.error import HTTPError, URLError

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
GDSC2_COMPOUNDS = (
    REPO_ROOT / "data" / "raw" / "gdsc2_sarcoma" / "gdsc2" / "screened_compounds_rel_8.5.csv"
)
SORAGNI_DRUG_SCREEN = REPO_ROOT / "data" / "raw" / "soragni" / "tables" / "drug_screen.parquet"
OUTPUT_DIR = REPO_ROOT / "data" / "static"
OUTPUT_PARQUET = OUTPUT_DIR / "drug_xref.parquet"

# scripts/download/_utils.py is the source of truth for the manifest schema; reuse it.
sys.path.insert(0, str(REPO_ROOT / "scripts" / "download"))
from _utils import sha256_file, write_manifest  # noqa: E402

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
UNICHEM_BASE = "https://www.ebi.ac.uk/unichem/rest"
UNICHEM_DRUGBANK_SRC_ID = "2"  # UniChem source ID for DrugBank
REQUEST_SLEEP = 0.25  # seconds between API calls; PubChem allows 5 req/s
TIMEOUT = 15


def http_get_json(url: str) -> dict | list | None:
    req = urllib.request.Request(url, headers={"User-Agent": "fm-pdo-evaluator/0.1"})
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
            return json.loads(resp.read())
    except HTTPError as e:
        if e.code == 404:
            return None
        print(f"  [http] {url} -> {e.code} {e.reason}", file=sys.stderr)
        return None
    except (URLError, json.JSONDecodeError, TimeoutError) as e:
        print(f"  [http] {url} -> {type(e).__name__}: {e}", file=sys.stderr)
        return None


def pubchem_name_to_cid(name: str) -> int | None:
    encoded = urllib.parse.quote(name, safe="")
    data = http_get_json(f"{PUBCHEM_BASE}/compound/name/{encoded}/cids/JSON")
    time.sleep(REQUEST_SLEEP)
    if not data or "IdentifierList" not in data:
        return None
    cids = data["IdentifierList"].get("CID") or []
    if not cids:
        return None
    return int(min(cids))  # smallest = canonical parent (PubChem convention)


def pubchem_cid_to_inchikey(cid: int) -> str | None:
    data = http_get_json(f"{PUBCHEM_BASE}/compound/cid/{cid}/property/InChIKey/JSON")
    time.sleep(REQUEST_SLEEP)
    if not data or "PropertyTable" not in data:
        return None
    props = data["PropertyTable"].get("Properties") or []
    if not props:
        return None
    return props[0].get("InChIKey")


def unichem_inchikey_to_drugbank(inchikey: str) -> str | None:
    data = http_get_json(f"{UNICHEM_BASE}/inchikey/{inchikey}")
    time.sleep(REQUEST_SLEEP)
    if not data or not isinstance(data, list):
        return None
    for entry in data:
        if str(entry.get("src_id")) == UNICHEM_DRUGBANK_SRC_ID:
            return entry.get("src_compound_id")
    return None


def load_gdsc2_drugs() -> pd.DataFrame:
    """Return GDSC2 compounds with cols: drug_id, drug_name, synonyms, pubchem_cid (nullable)."""
    df = pd.read_csv(GDSC2_COMPOUNDS)

    # PUBCHEM column has values like "5291", "5291,176870", "none", "-", or NaN.
    def parse_cid(val: object) -> int | None:
        if pd.isna(val):
            return None
        s = str(val).strip().lower()
        if s in ("none", "-", "", "nan"):
            return None
        first = s.split(",")[0].strip()
        try:
            return int(first)
        except ValueError:
            return None

    return pd.DataFrame(
        {
            "drug_id": df["DRUG_ID"].astype(str),
            "drug_name": df["DRUG_NAME"].astype(str),
            "synonyms": df.get("SYNONYMS", pd.Series([""] * len(df))).fillna("").astype(str),
            "pubchem_cid": df.get("PUBCHEM", pd.Series([None] * len(df))).apply(parse_cid),
        }
    )


def load_soragni_drugs() -> list[str]:
    df = pd.read_parquet(SORAGNI_DRUG_SCREEN)
    return sorted(df["Drug_Name"].dropna().astype(str).unique().tolist())


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="ignore any existing parquet and rebuild from scratch",
    )
    args = parser.parse_args()

    if not GDSC2_COMPOUNDS.exists():
        sys.exit(
            f"[fail] missing {GDSC2_COMPOUNDS.relative_to(REPO_ROOT)} "
            "-- run scripts/download/download_gdsc2_sarcoma.py first"
        )
    if not SORAGNI_DRUG_SCREEN.exists():
        sys.exit(
            f"[fail] missing {SORAGNI_DRUG_SCREEN.relative_to(REPO_ROOT)} "
            "-- run scripts/download/download_soragni.py first"
        )

    if OUTPUT_PARQUET.exists() and not args.refresh:
        sys.exit(
            f"[skip] {OUTPUT_PARQUET.relative_to(REPO_ROOT)} already exists; "
            "pass --refresh to rebuild"
        )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    gdsc2 = load_gdsc2_drugs()
    soragni_names = load_soragni_drugs()
    print(f"[load] {len(gdsc2)} GDSC2 compounds, {len(soragni_names)} unique Soragni drugs")

    rows: list[dict] = []
    cid_cache: dict[str, int | None] = {}  # name (lower) -> CID
    inchikey_cache: dict[int, str | None] = {}
    drugbank_cache: dict[str, str | None] = {}

    # Pass 1: GDSC2 compounds
    for _i, row in gdsc2.iterrows():
        name = row["drug_name"]
        cid = row["pubchem_cid"]
        method = "gdsc2_pubchem_column" if cid else None
        if not cid:
            print(f"[gdsc2] resolving by name: {name}")
            cid = pubchem_name_to_cid(name)
            method = "gdsc2_name_lookup" if cid else "unresolved"
        cid_cache[name.lower()] = cid
        for syn in [s.strip() for s in row["synonyms"].split(",") if s.strip()]:
            cid_cache.setdefault(syn.lower(), cid)

        inchikey = inchikey_cache.get(cid) if cid else None
        if cid and inchikey is None and cid not in inchikey_cache:
            inchikey = pubchem_cid_to_inchikey(cid)
            inchikey_cache[cid] = inchikey
        drugbank = drugbank_cache.get(inchikey) if inchikey else None
        if inchikey and drugbank is None and inchikey not in drugbank_cache:
            drugbank = unichem_inchikey_to_drugbank(inchikey)
            drugbank_cache[inchikey] = drugbank
        rows.append(
            {
                "input_name": name,
                "source": "gdsc2",
                "source_drug_id": row["drug_id"],
                "pubchem_cid": cid,
                "inchikey": inchikey,
                "drugbank_id": drugbank,
                "resolution_method": method,
                "notes": None,
            }
        )

    # Pass 2: Soragni drugs (re-use GDSC2 resolution where possible)
    for name in soragni_names:
        key = name.lower()
        cached = cid_cache.get(key)
        if cached:
            cid = cached
            method = "soragni_via_gdsc2_synonym"
        else:
            print(f"[soragni] resolving by name: {name}")
            cid = pubchem_name_to_cid(name)
            method = "soragni_name_lookup" if cid else "unresolved"
            cid_cache[key] = cid

        inchikey = inchikey_cache.get(cid) if cid else None
        if cid and inchikey is None and cid not in inchikey_cache:
            inchikey = pubchem_cid_to_inchikey(cid)
            inchikey_cache[cid] = inchikey
        drugbank = drugbank_cache.get(inchikey) if inchikey else None
        if inchikey and drugbank is None and inchikey not in drugbank_cache:
            drugbank = unichem_inchikey_to_drugbank(inchikey)
            drugbank_cache[inchikey] = drugbank
        rows.append(
            {
                "input_name": name,
                "source": "soragni",
                "source_drug_id": None,
                "pubchem_cid": cid,
                "inchikey": inchikey,
                "drugbank_id": drugbank,
                "resolution_method": method,
                "notes": None,
            }
        )

    out = pd.DataFrame(rows)
    # Use nullable Int64 for CIDs so missing values survive parquet round-trip.
    out["pubchem_cid"] = out["pubchem_cid"].astype("Int64")
    out.to_parquet(OUTPUT_PARQUET, index=False)

    n_gdsc2 = (out["source"] == "gdsc2").sum()
    n_soragni = (out["source"] == "soragni").sum()
    n_cid = out["pubchem_cid"].notna().sum()
    n_ikey = out["inchikey"].notna().sum()
    n_db = out["drugbank_id"].notna().sum()
    overlap = out.dropna(subset=["pubchem_cid"]).groupby("pubchem_cid")["source"].nunique()
    n_overlap = (overlap > 1).sum()
    print()
    print(f"[summary] {len(out)} rows ({n_gdsc2} gdsc2 + {n_soragni} soragni)")
    print(f"          pubchem_cid resolved: {n_cid}/{len(out)}")
    print(f"          inchikey resolved:    {n_ikey}/{len(out)}")
    print(f"          drugbank_id resolved: {n_db}/{len(out)}")
    print(f"          unique CIDs in both Soragni AND GDSC2: {n_overlap}")

    digest = sha256_file(OUTPUT_PARQUET)
    manifest_path = OUTPUT_DIR / "manifest.json"
    manifest = {
        "dataset": "fmharness_static_assets",
        "release": {"asset_set": "drug_xref_v1"},
        "files": {
            "drug_xref.parquet": {
                "sha256": digest,
                "bytes": OUTPUT_PARQUET.stat().st_size,
                "source_uri": "built://scripts/build/build_drug_xref.py",
                "rows": int(out.shape[0]),
                "cols": int(out.shape[1]),
                "columns": list(out.columns),
            }
        },
    }
    write_manifest(manifest_path, manifest)
    print(f"[done] {OUTPUT_PARQUET.relative_to(REPO_ROOT)} ({digest[:12]}...)")


if __name__ == "__main__":
    main()
