#!/usr/bin/env python3
"""Download GDSC2 drug-response + CCLE/DepMap RNA-seq for the sarcoma subset.

Two open data sources, paired on COSMIC ID:

1. GDSC2 (Sanger Institute, release 8.5, dated 27Oct23). Fitted dose-response
   (AUC, IC50), compound metadata, cell-line metadata. Hosted on the public
   Sanger CDN; no auth required.

2. DepMap (Broad Institute, latest public release, ~quarterly). RNA-seq raw
   read counts (per protein-coding gene, RSEM-style; fed through DESeq2
   median-of-ratios in the loader to match Soragni's normalization scheme)
   plus model metadata mapping DepMap ACH-XXXXX IDs to COSMIC IDs. DepMap
   distributes through figshare collections that change per release (26Q1,
   25Q4, ...). To swap to a newer release:

   Go to https://depmap.org/portal/download/all/, pick the latest release,
   right-click the file and copy the figshare ndownloader link, then update
   the corresponding entry in DEPMAP_FILES below.

Outputs land in data/raw/gdsc2_sarcoma/ with sha256 of each file recorded in
manifest.json. The manifest schema is shared with download_soragni.py via
scripts/download/_utils.py.

Usage:
    python scripts/download/download_gdsc2_sarcoma.py
    python scripts/download/download_gdsc2_sarcoma.py --verify-only
"""

from __future__ import annotations

import argparse
import sys
import urllib.request
from pathlib import Path
from urllib.error import HTTPError, URLError

from _utils import (
    MANIFEST_NAME,
    load_manifest,
    sha256_file,
    skip_or_fail_on_hash,
    verify_manifest,
    write_manifest,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = REPO_ROOT / "data" / "raw" / "gdsc2_sarcoma"
MANIFEST_PATH = OUTPUT_DIR / MANIFEST_NAME

# GDSC2 release 8.5 (27Oct23) -- verified live on the Sanger CDN.
# Bump these when Sanger publishes a new release.
GDSC2_RELEASE = "8.5"
GDSC2_BASE = f"https://cog.sanger.ac.uk/cancerrxgene/GDSC_release{GDSC2_RELEASE}"
GDSC2_FILES: list[tuple[str, str]] = [
    (
        "GDSC2_fitted_dose_response_27Oct23.xlsx",
        f"{GDSC2_BASE}/GDSC2_fitted_dose_response_27Oct23.xlsx",
    ),
    ("screened_compounds_rel_8.5.csv", f"{GDSC2_BASE}/screened_compounds_rel_8.5.csv"),
    ("Cell_Lines_Details.xlsx", f"{GDSC2_BASE}/Cell_Lines_Details.xlsx"),
]

# DepMap files -- release 26Q1 (canonical-id snapshot public-26q1-5bbf).
# URLs hit depmap.org which 302s to a signed Google Cloud Storage URL that's
# regenerated per request, so the depmap.org URL itself is stable across runs.
#
# RawReadCount is RSEM expected counts per protein-coding gene (linear,
# integer-ish); the loader runs DESeq2 median-of-ratios on this to match
# Soragni's pre-computed normalized counts. The log2(TPM+1) variant at the
# same canonical-id (.27) is available if a length-normalized view is ever
# needed for a sensitivity row; not pulled by default.
DEPMAP_RELEASE = "26Q1"
DEPMAP_BASE = "https://depmap.org/portal/download/api/download"
DEPMAP_FILES: list[tuple[str, str]] = [
    (
        "OmicsExpressionRawReadCountHumanProteinCodingGenes.csv",
        f"{DEPMAP_BASE}?file_name=downloads-by-canonical-id%2Fpublic-26q1-5bbf.27%2FOmicsExpressionRawReadCountHumanProteinCodingGenes.csv&dl_name=OmicsExpressionRawReadCountHumanProteinCodingGenes.csv&bucket=depmap-external-downloads",
    ),
    (
        "Model.csv",
        f"{DEPMAP_BASE}?file_name=downloads-by-canonical-id%2Fpublic-26q1-5bbf.37%2FModel.csv&dl_name=Model.csv&bucket=depmap-external-downloads",
    ),
]


def default_manifest() -> dict:
    return {
        "dataset": "gdsc2_sarcoma",
        "release": {"gdsc": GDSC2_RELEASE, "depmap": DEPMAP_RELEASE},
        "files": {},
    }


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": "fm-pdo-evaluator/0.1"})
    with urllib.request.urlopen(req) as response, tmp.open("wb") as fh:
        while True:
            chunk = response.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
    tmp.replace(dest)


def fetch_with_sha(name: str, url: str, manifest: dict) -> None:
    dest = OUTPUT_DIR / name
    if skip_or_fail_on_hash(name, dest, manifest["files"]):
        print(f"[skip] {name} present with matching sha256")
        return
    print(f"[get ] {name} <- {url}")
    try:
        download(url, dest)
    except (HTTPError, URLError) as e:
        sys.exit(f"[fail] {name} download error: {e}")
    digest = sha256_file(dest)
    manifest["files"][name] = {
        "sha256": digest,
        "bytes": dest.stat().st_size,
        "source_uri": url,
    }
    print(f"       sha256 {digest}")


def fetch_gdsc2(manifest: dict) -> None:
    for name, url in GDSC2_FILES:
        fetch_with_sha(f"gdsc2/{name}", url, manifest)


def fetch_depmap(manifest: dict) -> None:
    for name, url in DEPMAP_FILES:
        fetch_with_sha(f"depmap/{name}", url, manifest)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--verify-only", action="store_true", help="re-verify existing files; no download"
    )
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = load_manifest(MANIFEST_PATH, default_manifest())

    if args.verify_only:
        errors = verify_manifest(OUTPUT_DIR, MANIFEST_PATH)
        sys.exit(1 if errors else 0)

    fetch_gdsc2(manifest)
    fetch_depmap(manifest)
    expected = {f"gdsc2/{n}" for n, _ in GDSC2_FILES} | {f"depmap/{n}" for n, _ in DEPMAP_FILES}
    manifest["files"] = {k: v for k, v in manifest["files"].items() if k in expected}
    write_manifest(MANIFEST_PATH, manifest)
    print(f"[done] manifest written to {MANIFEST_PATH.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
