# fm-pdo-harness — 3-Week MVP Plan

**Owner:** Lucas Gillenwater (Greene Lab, CU Anschutz)
**Window:** May 18 – Jun 5, 2026 (15 working days; Memorial Day Mon May 25 worked)
**Repo (to create):** `greenelab/fm-pdo-evaluator` forked to dev repo`lagillenwater/fm-pdo-evaluator`
**Compute platform:** CU Anschutz Alpine HPC
**Status:** Day 1 (Mon May 18) missed; catching up Day 1 + Day 2 in a compressed session today.

---

## 1. Overview

Build the minimum viable harness needed to produce  **preliminary result** evaluating three models × three OOD splits × two PDTO datasets. Include positive and negative controls and a pretraining-leakage exposure profile per result. The analysis tests whether **in-distribution model rankings predict out-of-distribution model rankings** on real PDTO drug response data.

## 2. Thesis

Biological foundation models and world models are heading toward commoditization. The clinical-utility bottleneck is shifting from training to *trust* — which predictions should researchers and clinicians believe? The empirical claim the harness exists to test: **in-distribution metrics do not predict OOD performance, so OOD generalization tests on prospective and rare-subtype data are the only meaningful trust signal.**


## 3. Deliverable

A registry-backed Markdown report with:

- **3 × 3 × 2 = 18 metric cells**, each with bootstrap CIs and per-subtype stratification:
  - Models: linear baseline (on raw expression), Tahoe-x1, STATE
  - Splits: stratified in-distribution, leave-patient-out, leave-subtype-out
  - Datasets: Soragni 2024 sarcoma PDTOs, Yang 2024 primary liver cancer PDOs (399 organoids from 144 patients; 7-drug panel)
- **Controls:** negative (random, population-average, patient-ID-shuffled, gene-shuffled) and positive (canonical biomarker per dataset)
- **Diagnostics:** patient-ID-shuffle null on every cell; determinism check
- **Leakage exposure profile** per cell: drug overlap with Tahoe-100M, declared-corpus dataset overlap, subtype prevalence
- **Rank-correlation table:** Spearman across (model rankings) for every split-pair, per dataset and cross-dataset

## 4. Checklist

| # | Item | Status |
|---|---|---|
| 1 | Create `greenelab/fm-pdo-evaluator` and fork to dev repo`lagillenwater/fm-pdo-evaluator` repo on GitHub | Done |
| 2 | Confirm Synapse access to `synapse.org/PDTOSarcoma` (Soragni 2024) | Confirmed |
| 3 | Confirm Yang 2024 data access path (GSA accession for RNA-seq/WES + supplementary table for drug response) | Confirmed |
| 4 | HuggingFace personal access token (for Tahoe-x1 weights) | To do |
| 5 | Pre-commit, pre-pyright, pre-pytest configuration baked into template | To do |
| 6 | Apptainer module identified on Alpine; first test image builds end-to-end | To do |
| 7 | HF token + Synapse PAT staged at `~/.fmharness/secrets` on Alpine; verified readable by sbatch | To do |
| 8 | STAR 2.7.9a genome index built against GRCh38 + GENCODE on Alpine; `quantify_rnaseq.sbatch` tested end-to-end (STAR → GeneCounts → median-of-ratios via `pydeseq2`) on 1 Soragni FASTQ pair | To do |

## 5. Tech stack

| Concern | Choice |
|---|---|
| Language | Python 3.11+ |
| Env / deps | `uv` + `pyproject.toml` |
| Data container | `anndata`; per-sample drug-response stored as a structured `pandas` frame (sparse panels, not matrix-shaped) |
| Schema | `pydantic` v2 |
| Preprocessing | `scanpy` |
| RNA quantification | STAR 2.7.9a (full alignment) with `--quantMode GeneCounts` → per-sample gene counts → median-of-ratios normalization via `pydeseq2`. Matches Soragni 2024 protocol (STAR + GRCh38; median-ratio normalization per Anders & Huber 2010 / Love 2014). Reference: GRCh38 + GENCODE (FASTA + GTF listed in `data/static/manifest.json` with sha256; output tranche content hash captured by `EnvironmentSnapshot.data_commit`). Run as a Slurm array job on Alpine; one cached expression matrix per dataset under `data/tranches/{tranche}/expression.parquet`. Applied identically to both Soragni and Yang for embedding comparability. |
| Probe | `scikit-learn` (StandardScaler → ElasticNetCV / LogisticRegressionCV) |
| Bootstrap | `scipy.stats.bootstrap`, parallel via `joblib` |
| Foundation models | `transformers` (HF), `arc-state` (PyPI) |
| CLI | `typer` |
| Config | `pydantic-settings` + YAML |
| Registry | Parquet (rows) + JSON manifest (runs) |
| Tests | `pytest` + `hypothesis` |
| Lint / type | `ruff`, `pyright` |
| Container | Apptainer (Singularity-compatible) for HPC inference; one image per model wrapper if torch conflicts |
| Determinism | `torch.use_deterministic_algorithms(True)`, `CUBLAS_WORKSPACE_CONFIG=:4096:8`, seeded python/numpy/torch RNG via `fmharness.utils.determinism` |
| Secrets | `.env.example` in repo; real `.env` gitignored; Alpine sbatch sources `~/.fmharness/secrets` |
| Static assets | `data/static/manifest.json` with sha256 per file; loaders verify on read |
| CI | GitHub Actions |
| License | BSD-2-Clause Plus Patent |

## 6. Compute environment (Alpine HPC)

- **Development:** local Mac for code + CI; push to Alpine for inference and bootstrap runs.
- **Inference partitions:** Alpine GPU partition with one A100 40GB or A10 24GB per inference job. Tahoe-x1 and STATE each fit comfortably in 24GB at fp16 with batch size 1.
- **Storage:** Alpine scratch for tranches and model weights (~30GB MVP, ~60GB headroom). Permanent storage in Alpine project space for the registry; Parquet artifacts are small (<100MB total).
- **Batch jobs:** Slurm sbatch scripts under `scripts/alpine/` for embedding generation and bootstrap runs. Interactive development for everything else.
- **Compliance:** Synapse data is staged only on Alpine, never to laptop disk or cloud. DUA-permitted use case documented in `docs/data_governance.md`.
- **RNA quantification (FASTQ → expression matrix):** Both Soragni (Synapse) and Yang (GSA, expected) deposit raw FASTQ, not expression matrices. Pipeline matches the Soragni 2024 protocol: STAR 2.7.9a alignment to GRCh38 + GENCODE → `--quantMode GeneCounts` for per-sample gene counts → median-of-ratios normalization via `pydeseq2`. Runs as a Slurm array job (`scripts/alpine/quantify_rnaseq.sbatch`), parallelized across samples. STAR is RAM-heavy (~30GB / sample) and substantially slower than pseudoalignment — budget Alpine large-memory partition time accordingly. Output gene-level normalized counts cached under `data/tranches/{tranche}/expression.parquet` and content-hashed into the tranche manifest. Reference (GRCh38 + GENCODE version) pinned via `data/static/manifest.json` sha256. **Kick off at start of Week 1** so completion precedes the Day 4 / Day 5 loaders.
- **Container strategy:** Apptainer images under `containers/` for foundation-model inference. Tahoe-x1 image built Day 8; STATE either reuses it or gets its own image if torch/CUDA conflicts force a split (see risk register). Container digest is captured in every prediction record.
- **GPU determinism:** All CLI entrypoints call `fmharness.utils.determinism.fix_seeds(seed)`, which seeds python/numpy/torch and sets `torch.use_deterministic_algorithms(True)` + `CUBLAS_WORKSPACE_CONFIG=:4096:8`. Without this, success criterion #6 cannot hold on GPU.
- **Secrets management:** HuggingFace token (Tahoe weights) and Synapse PAT live in `~/.fmharness/secrets` on Alpine, sourced by sbatch headers. `.env.example` in repo documents the variable names; real `.env` is gitignored.
- **Static asset versioning:** Tahoe-100M drug list, drug crosswalk tables, gene-panel reconciliation tables tracked in `data/static/manifest.json` with sha256 per file. Loaders verify on read and refuse to proceed on mismatch.
- **Provenance via `EnvironmentSnapshot`:** Every `PredictionRecord` embeds `code_commit`, `container_digest`, `python_version`, `torch_version`, `cuda_version`, `model_weights_hash`, `data_commit`, `seed`, `cuda_deterministic`. Without this the determinism check is half-blind to node-to-node drift.

## 7. Week-by-week day-by-day


### Week 1 — schema, splits, two-dataset ingestion

| Day | Date | Goal | Concrete deliverables |
|---|---|---|---|
| 1+2 | Tue May 19 | Scaffolding **and** schema v1 | Repo cloned, `pyproject.toml`, `uv.lock` committed, license, README, `.gitignore`, `.pre-commit-config.yaml`, GitHub Actions CI (pytest + ruff + pyright). `.env.example` with documented variable names. `containers/fmharness.def` skeleton (full image built Day 8). Pydantic models in `src/fmharness/schema/`: `Patient`, `Sample`, `DrugAssay`, `BaselineExpression`, `Prediction`, `Tranche`, `LeakageProfile`, `EnvironmentSnapshot`. `src/fmharness/utils/determinism.py` with `fix_seeds(seed)` seeding python/numpy/torch + setting deterministic-algorithms and CUBLAS workspace env. `docs/environment.md` covering container strategy, secrets layout, and determinism contract. `tests/test_schema.py` with hypothesis property tests. First PR merged. |
| 3 | Wed May 20 | Drug crosswalk | `src/fmharness/data/drug_xref.py`: name → InChIKey → DrugBank ID → PubChem CID. Mapping tables for the union of Soragni and Yang panels stored under `data/static/` with sha256 entries in `data/static/manifest.json`; loader verifies hashes on read. Tests for ambiguous names and for manifest-mismatch refusal. |
| 4 | Thu May 21 | Soragni loader | Soragni RNA-seq on Synapse is raw FASTQ; `quantify_rnaseq.sbatch` (kicked off start-of-week) must complete by today and produce `data/tranches/soragni/expression.parquet`. `src/fmharness/data/loaders/soragni.py`: read cached expression matrix + Synapse drug-response tables into AnnData + `DrugAssay` frame, emit `Tranche(tranche_id='soragni_pdo_sarcoma_2024', ...)`. QC: gene filtering, response sanity, drug-ID resolution. Report per-subtype patient counts and panel-coverage density to `docs/datasets.md`. |
| 5 | Fri May 22 | Yang loader + Week-1 review | Yang RNA-seq on GSA is a .rds file. `src/fmharness/data/loaders/yang.py`: read cached `data/tranches/yang/expression.parquet`, emit `yang_pdo_liver_2024` tranche. Handle multiple organoids per patient (399 / 144 ≈ 2.77 per patient on average; leave-patient-out must hold out all PDOs from a held-out patient). Document HCC / ICC / cHCC-ICC subtype taxonomy and missing-drug handling across the 7-drug panel. End-of-week PR. CI green. |

**End-of-Week-1 success criterion:** `fmharness ingest --source soragni` and `fmharness ingest --source yang` produce validated tranches; schema locked.

### Week 2 — splits, adapter, first two models (5 days; Memorial Day worked)

| Day | Date | Goal | Concrete deliverables |
|---|---|---|---|
| 6 | Mon May 25 | OOD splits | `src/fmharness/splits/`: `StratifiedInDistribution`, `LeavePatientOut`, `LeaveSubtypeOut` (configurable granularity — fine / coarse). Property tests: no (patient, subtype) overlap across folds for either dataset. Splitters accept and persist a seed; seed flows into the `EnvironmentSnapshot` attached to downstream prediction records. `Evaluator` refuses to run without a named split. |
| 7 | Tue May 26 | Adapter contract + linear baseline | `src/fmharness/models/adapter.py`: `ModelAdapter` Protocol. Required `metadata()`: `pretraining_corpus`, `pretraining_cutoff_date`, `task_signal_in_pretrain`. `MockAdapter` for tests. `src/fmharness/models/wrappers/linear_baseline.py`: sklearn pipeline on normalized expression. First end-to-end model run on Soragni + Yang. |
| 8 | Wed May 27 | Tahoe-x1 wrapper (part 1) | **Build `containers/tahoe.def`** — pinned torch + CUDA matching Alpine GPU driver; image digest pinned in `containers/digests.json`. HuggingFace loading on Alpine GPU partition (token sourced from `~/.fmharness/secrets`). Gene-panel reconciliation. Bulk-input handling decision documented (treat sample as single observation; pseudobulk path as plan B). Pretraining metadata: declare Tahoe-100M corpus + ingest ~1,100-drug list as static asset under `src/fmharness/leakage/corpora/tahoe_100m_drugs.parquet` (sha256 registered in `data/static/manifest.json`). Smoke test on 5 Soragni samples runs inside the container. |
| 9 | Thu May 28 | Tahoe-x1 wrapper (part 2) + probe | Finish `embed()` with disk caching keyed by input + model version. `predict_native()` if Tx1 exposes a drug-aware head. `src/fmharness/probe/linear.py`: fixed probe architecture (StandardScaler → ElasticNetCV / LogisticRegressionCV) applied identically across all models. |
| 10 | Fri May 29 | First end-to-end matrix slice | `src/fmharness/tasks/sensitivity.py`: binary responder / non-responder. Run **(linear baseline, Tahoe-x1) × (in-distribution, LPO) × (Soragni, Yang)** = 8 cells. Predictions written to disk under the registry layout. End-of-week PR. |

**End-of-Week-2 success criterion:** 8 of 18 matrix cells completed and reproducible from CLI; embeddings cached on Alpine.

### Week 3 — STATE, metrics, registry, leakage scan, report

| Day | Date | Goal | Concrete deliverables |
|---|---|---|---|
| 11 | Mon Jun 1 | STATE wrapper | `src/fmharness/models/wrappers/state.py`: `arc-state` PyPI integration on Alpine GPU. Decide whether STATE reuses the Tahoe container or needs its own `containers/state.def`; document the rationale in `docs/environment.md`. Image digest pinned in `containers/digests.json` either way. SE embedding + ST transition. Pretraining metadata declaration. Smoke test. Run STATE × 3 splits × 2 datasets → 6 cells, completing the matrix to 14/18. Add leave-subtype-out for all three models → final 4 cells. **All 18 cells produced by end of day.** |
| 12 | Tue Jun 2 | Decision metrics + controls | `src/fmharness/metrics/`: `top_k_hit_rate`, `regret`, `brier_score`, `expected_calibration_error`. `src/fmharness/metrics/stratified.py`: bootstrap CIs (n=1000) per subtype, per dataset, parallel via joblib on Alpine CPU node. `src/fmharness/controls/negative.py`: random, population-average, patient-ID-shuffled, gene-shuffled. `src/fmharness/controls/positive.py`: HCC/lenvatinib (and any other multikinase-inhibitor / standard-of-care pair present in Yang's 7-drug panel — e.g. HCC/sorafenib, HCC/regorafenib, ICC/gemcitabine if present) for Yang; NTRK/larotrectinib, EZH2/tazemetostat, mTOR/sirolimus for Soragni — only those where both biomarker label and drug are present in the dataset. Final pairings finalized at loader time on Day 5 once the Yang panel is resolved. |
| 13 | Wed Jun 3 | Tranche-aware registry + leakage scan | `src/fmharness/registry/predictions.py`: `PredictionRecord(model_version, tranche_id, split_spec, prediction_hash, metric_values, leakage_exposure_profile, environment_snapshot)`. Content hash. Immutability test extended to cover `environment_snapshot`. `src/fmharness/leakage/scan.py`: drug overlap against Tahoe-100M list, declared-corpus overlap, suspected-subtype prevalence. Re-run all 18 cells through the registry to attach `LeakageProfile`s and `EnvironmentSnapshot`s. |
| 14 | Thu Jun 4 | Controls + rank-correlation analysis | Execute 4 negative-control types × 2 datasets = 8 control rows. Execute positive-control oracles per dataset. Run patient-ID-shuffle on every cell (18 shuffle rows). **Determinism check:** re-run one cell from inside the pinned container; confirm prediction hash identical to first run; fail loud on mismatch. Compute Spearman rank correlation of model rankings across split-pairs (in-distribution ↔ LPO, in-distribution ↔ LSO, LPO ↔ LSO) per dataset and cross-dataset. |
| 15 | Fri Jun 5 | Report + write-up + tag | `src/fmharness/report/markdown.py`: render `reports/preliminary.md` from the registry. Sections: matrix table with CIs, rank-correlation table, control results, leakage exposure summary, narrative on what flipped vs. what only weakened. Internal preliminary-results write-up (3–5 pages) for the proposal. README finalized. Tag `v0.1.0`. |

**End-of-Week-3 success criterion:** `git clone && uv sync && fmharness run --all-configs && fmharness report --output reports/preliminary.md` produces the headline analysis on Alpine. Report is shareable with collaborators on Day 15.

## 8. Success criteria

| # | Criterion | Measurable |
|---|---|---|
| 1 | Both datasets ingested and validated | `fmharness ingest` exits 0 for both sources |
| 2 | All 18 matrix cells produce predictions | 18 `PredictionRecord`s in the registry, each with a `LeakageProfile` |
| 3 | Controls run | ≥8 negative-control rows, ≥4 positive-control rows |
| 4 | Bootstrap CIs reported | Every metric row has a 95% CI |
| 5 | Rank correlation across splits computed | `rank_correlation` section in `preliminary.md` with values + CIs |
| 6 | Determinism | Re-running the same config produces identical prediction hashes |
| 7 | Leakage scan non-trivial | Tahoe-x1 drug-overlap fraction reported per dataset |
| 8 | CI green on `main` | GitHub Actions passing |
| 9 | Coverage | ≥80% on `schema/`, `splits/`, `registry/`, `metrics/` |
| 10 | Environment provenance | Every prediction record carries an `EnvironmentSnapshot` with non-null `container_digest`, `code_commit`, `torch_version`, `cuda_version`, `model_weights_hash`, and `seed` |

## 9. Repository structure

```
fm-pdo-harness/
├── pyproject.toml
├── LICENSE
├── README.md
├── .github/workflows/ci.yml
├── .pre-commit-config.yaml
├── .gitignore
├── .env.example
├── containers/
│   ├── fmharness.def
│   ├── tahoe.def
│   ├── state.def              # built only if torch conflicts force a split
│   └── digests.json           # pinned image digests
├── configs/
│   ├── linear_baseline_{soragni,yang}_{id,lpo,lso}.yaml    # 6 files
│   ├── tahoe_x1_{soragni,yang}_{id,lpo,lso}.yaml           # 6 files
│   └── state_{soragni,yang}_{id,lpo,lso}.yaml              # 6 files
├── scripts/alpine/                                          # sbatch scripts
│   ├── quantify_rnaseq.sbatch     # FASTQ → gene counts (STAR 2.7.9a + --quantMode GeneCounts) → median-of-ratios (pydeseq2)
│   ├── embed_tahoe_x1.sbatch
│   ├── embed_state.sbatch
│   └── bootstrap_metrics.sbatch
├── src/fmharness/
│   ├── schema/                # pydantic models
│   ├── data/loaders/          # soragni.py, yang.py
│   ├── data/drug_xref.py
│   ├── splits/                # stratified, lpo, lso
│   ├── models/adapter.py
│   ├── models/wrappers/       # linear_baseline.py, tahoe_x1.py, state.py
│   ├── probe/linear.py
│   ├── tasks/sensitivity.py
│   ├── metrics/               # decision, calibration, stratified
│   ├── controls/              # negative.py, positive.py
│   ├── registry/predictions.py
│   ├── leakage/               # scan.py, corpora/ (static assets)
│   ├── utils/determinism.py   # fix_seeds(); deterministic algos + CUBLAS env
│   ├── report/markdown.py
│   └── cli.py
├── tests/                     # mirrors src/ structure
├── notebooks/
│   ├── 00_smoke_test.ipynb
│   └── 10_split_flip_analysis.ipynb
├── data/tranches/             # cached tranche artifacts (gitignored)
├── data/static/               # versioned static assets + manifest.json (sha256)
├── reports/                   # generated reports (gitignored)
└── docs/
    ├── datasets.md            # subtype taxonomies, panel coverage, decisions
    ├── adapter_contract.md    # how to add a new model
    ├── leakage_methodology.md # what the scan catches and what it does not
    ├── environment.md         # container strategy, secrets, determinism contract
    └── data_governance.md     # DUA / Synapse / Alpine storage rules
```

## 10. Risks

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Tahoe-x1 input format requires heavy preprocessing | High | Medium | Allocate Wed–Thu Week 2. Fall back to scFoundation (bulk-friendly) if blocked >1 day |
| STATE installation has CUDA/torch conflicts on Alpine | Medium | Medium | Isolated optional dependency group `fm-pdo-harness[state]`. Fall back to scGPT or Geneformer |
| PDTO bulk RNA-seq vs single-cell FM training mismatch produces noisy embeddings | Medium | High | Document explicitly. The finding is itself publishable. Sensitivity analysis: pseudobulk for one model as a comparison row |
| Yang subtype taxonomy (HCC / ICC / cHCC-ICC) needs harmonization across 144 patients | Medium | Low | One-day buffer in Week 1; harmonize during loader work |
| RNA quantification pipeline (STAR alignment on ~500+ FASTQ pairs across Soragni + Yang; ~30GB RAM per sample) overruns the Day 4–5 window | High | High | Kick off Slurm array job at the start of Week 1, in parallel with scaffolding/schema work. STAR is substantially heavier than pseudoalignment — request Alpine large-memory nodes and parallelize aggressively. Real risk is gating on STAR genome index build (~1hr one-time) or queue delays under RAM contention. If a co-author/collaborator has a published Partek-Flow count matrix for Soragni, use that as a fallback while the pipeline runs |
| Yang's 7-drug panel collapses the Tahoe-100M leakage-overlap analysis (likely 7/7 or 6/7 in Tahoe corpus) | High | Medium | Report the overlap honestly; treat leakage scan on Yang as a binary disclosure rather than a signal axis; rely on Soragni's larger panel for the per-drug leakage gradient |
| Liver biomarker positive controls don't have matching drugs in Yang's 7-drug panel | Medium | Low | Pick what's testable (lenvatinib for HCC is confirmed); report what wasn't |
| Sarcoma biomarker positive controls don't have matching drugs in Soragni panel | Medium | Low | Pick what's testable; report what wasn't |
| Bootstrap CI computation too slow on a single node | Low | Low | Parallelize across Alpine CPU cores via joblib. Drop to n=500 during dev |
| Split-flip doesn't happen biologically | Medium | Medium | Report rank-correlation values quantitatively. Spearman dropping from 0.9 → 0.4 is the same thesis |
| Compressed Day 1 + Day 2 overruns into Wed | Medium | Low | Drug crosswalk is small and can defer to Day 4 if needed |
| GPU non-determinism breaks success criterion #6 | Medium | High | `torch.use_deterministic_algorithms(True)` + `CUBLAS_WORKSPACE_CONFIG=:4096:8` + seeded RNG enforced via `fmharness.utils.determinism`. Determinism check on Day 14 catches regressions |
| Alpine node-to-node embedding drift inside the same model | Medium | High | Apptainer image per model wrapper, pinned by digest in `containers/digests.json`. `EnvironmentSnapshot` records the digest on every prediction |
| HF token or Synapse PAT leaks into repo or laptop disk | Low | High | Secrets live only at `~/.fmharness/secrets` on Alpine; `.env` gitignored; pre-commit hook scans for token patterns |

## 11. Deferrals (explicit MVP non-goals)

| Item | Why deferred |
|---|---|
| Substrate-gap result (cell-line → PDTO) | Requires GDSC/PRISM ingestion + dose/readout harmonization |
| Drug-ranking task | Sensitivity is enough to demonstrate split-flip |
| Transcriptional-response task | Not central to clinical-decision framing |
| Drug-resistance task | Needs CancerFoundation + resistance-specific data shape |
| Cell-line-proxy negative control | Requires GDSC/CCLE ingestion |
| AetherCell, PharmaFormer wrappers | Wait until split-flip result is in hand and competitor code is verified |
| Public dashboard / leaderboard UI | Registry is the source of truth; UI is a thin derivation, build later |
| Active wet-lab prospective tranche | Architecture supports it; first prospective run lands when refinedscience tranche arrives |

## 12. Positive control questions to answer during build. 

1. Which 2–3 sarcoma biomarker–drug pairs are present in the Soragni panel for positive controls
2. Which liver biomarker–drug pairs (beyond HCC/lenvatinib) are covered by Yang's 7-drug panel

## 13. References

- Al Shihabi et al. 2024, *Cell Stem Cell* — landscape of drug sensitivity in sarcoma (Soragni)
  - <https://www.cell.com/cell-stem-cell/fulltext/S1934-5909(24)00296-0>
  - Data: <https://www.synapse.org/PDTOSarcoma>
- Yang et al. 2024, *Cancer Cell* — pharmacogenomic profiling of intra-tumor heterogeneity using a large organoid biobank of liver cancer (399 PDOs from 144 patients; 7-drug screen)
  - <https://www.cell.com/cancer-cell/fulltext/S1535-6108(24)00089-8>
  - DOI: 10.1016/j.ccell.2024.03.004
- Tahoe-x1 / Tahoe-100M (Tahoe Bio) — <https://www.tahoebio.ai/>
- STATE (Arc Institute) — <https://arcinstitute.org/manuscripts/State> ; `arc-state` PyPI
- Wei et al. 2026, *Nature Methods* — generalizable single-cell perturbation benchmarking
- Vlachogiannis et al. 2018, *Science* — PDTOs predict clinical response
- Ooft et al. 2019, *Science Translational Medicine* — PDTOs predict chemotherapy response
