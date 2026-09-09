# Datasets

**As of** 2026-09-01.

The registry of every dataset this project reads: what each one is, where it came from, the
script and date of its download, and what was done to it afterwards. A task's `design.md` links
to the entry here. A dataset enters this registry with the first rung that reads it,[^1]
and its registration is a schema-validated `Tranche` record under `data/tranches/` — this
document is the index over those records.

[^1]: The archived predecessor of this registry — dataset notes from the pre-rebuild lineage —
is preserved verbatim at
[`docs/archive/2026-05-25-datasets-old-lineage.md`](archive/2026-05-25-datasets-old-lineage.md),
and gets reconciled into the entries of the rungs that read those datasets, when those rungs
arrive.

Corpus-level facts below are asserted from the source's own documentation and cited. Facts about
the *restricted pool a run actually consumed* — replicate depth per pair, dose levels, control
wells — are emitted by that run into its task folder as a pool-description table.

---

## Tahoe-100M

**What it is.** A giga-scale single-cell perturbation atlas from Vevo Therapeutics: over 100
million single-cell transcriptomes covering **50 cancer cell lines** exposed to **~1,100
small-molecule perturbations** (drug–dose conditions), generated on Vevo's pooled *Mosaic*
platform across **14 96-well plates**. Vehicle control is DMSO (`DMSO_TF`), plate-matched: every
plate carries its own DMSO wells, so treated-vs-control comparisons are within-plate.

**Source and citation.** Hugging Face [`tahoebio/Tahoe-100M`](https://huggingface.co/datasets/tahoebio/Tahoe-100M),
license CC0-1.0. Paper: Zhang et al., *Tahoe-100M: A Giga-Scale Single-Cell Perturbation Atlas
for Context-Dependent Gene Function and Cellular Modeling*, bioRxiv 2025,
[10.1101/2025.02.20.639398](https://doi.org/10.1101/2025.02.20.639398).

**What this project reads.** Not the 100M raw cells but two derived configurations:

- `pseudobulk_differential_expression` — the table every delta rung is built on. ~4.1 billion
  rows (89 GB, 1,026 parquet shards): one row per (cell line, drug, dose, plate) per gene,
  carrying the full DESeq2 result for the treated pseudobulk against the plate-matched DMSO
  pseudobulk. Sixteen columns: `gene_name`; `baseMean`, `log2FoldChange`, `lfcSE`, `stat`,
  `pvalue`, `padj`; `plate`, `n_cells_trt`, `n_cells_ctrl`; `Cell_ID_Cellosaur`,
  `Cell_ID_DepMap`, `Cell_Name_Vevo`; `drug`, `concentration`, `concentration_unit`. A gene
  DESeq2 could not test carries `baseMean` 0 and null in the five statistic columns, so
  untestable genes are absent by finiteness rather than by any filter of ours. `lfcSE` is the
  within-plate standard error of that one contrast and does not see plate-to-plate variation.
  These statistics are **accepted upstream** — this project does no re-quantification and no
  re-normalization of them. The registry records that fact; the decision behind it — reasons,
  risks, and what would revisit it — is made in the rung-0 design
  ([docs/tasks/rung0-assay-reliability/design.md](tasks/rung0-assay-reliability/design.md)).
- `drug_metadata` — used only to resolve drug names to PubChem CIDs.

**Download.** 2026-07-24, by `scripts/alpine/01_pseudobulk_shortcut.sbatch` (download logic
  extracted as `scripts/download_tahoe_pseudobulk_de.py`, based on the archived lineage's
  2026-07-24 pull; authenticated, single-process; shard set has no drug partition, so
  concurrent readers trip rate limit) to Alpine scratch at `/scratch/alpine/$USER/tahoe_pseudobulk_de`
  (verified present 2026-08-27, 83G). Scratch is purgeable; the same script rebuilds the copy (~12h).
Ingested as tranche `tahoe100m-pseudobulk-de.v1` (`data/tranches/`): the record carries the
Hugging Face dataset revision as its version and a content hash over the sorted per-shard
manifest, and the provenance record of any promoted result names that content hash as its
`data_commit`, per the environment contract (`docs/environment.md` §6).

**Measured shape of the table — smaller than the corpus above, and not a contradiction.** The
50 lines and ~1,100 perturbations stated under "What it is" are the atlas's, from the paper. The
derived `pseudobulk_differential_expression` table this project reads covers less of it, and the
difference matters for every count a rung reports. Counted on the cluster (job 31979673) rather
than inferred. **4,089,820,780 rows** and **3** dose levels (0.05, 0.5, 5.0 uM), both exact.
**Fifty** cell lines, exact; **379** drugs with a plate id, exact (the split-assignment scan,
job 32341610, 2026-09-09 — the HyperLogLog sketch had read 391, which is the couple of percent
it is good to); approximately **49,040** genes, still from the sketch. An earlier draft of this
entry reported 45 cell lines from that sketch, which the exact count contradicts.

**18,950** (line, drug) conditions and **56,827** (line, drug, dose) triples, all exact.

**Dose is very nearly confounded with plate, and it governs what rung 0 can measure.** 86.6% of
dose-conditions (49,186 of 56,827) sit on a **single plate**, so only **7,641** have the two
plates a split needs — and those are not spread evenly:

| dose | dose-conditions | with >= 2 plates | |
|---|---|---|---|
| 0.05 uM | 18,948 | 1,245 | 6.6% |
| 0.5 uM | 18,929 | 1,000 | 5.3% |
| 5.0 uM | 18,950 | **5,396** | **28.5%** |

So the replicated base is dominated by the top dose: 71% of it is 5.0 uM. It spans **121 drugs
and 50 cell lines**, so it is a broad slice of the screen rather than a few compounds.

DESeq2 could not test much of the table — **59%** of rows carry `baseMean` zero and **80%** a
null `padj` — so the gene set any statistic actually scores is far smaller than the gene count
above, and the per-condition counts in a result are the number to read.

**Restriction (later rungs, not rung 0).** A delta rung that must line up with GDSC2 reads only
the drugs shared between Tahoe and GDSC2, matched by PubChem CID: **32 drugs**
(`data/static/tahoe_target_cids.txt`, on Alpine in the repository checkout; not tracked on this
branch — generated by `scripts/alpine/00_target_cids.sbatch`). That list and the pool arithmetic
around it belong to the rung that needs them. Rung 0 does not use it: it measures at the assay's
full extent, and the superseded rung's "1,600 scored pairs" describes that restricted pool, not
this one.

**Processing after download, per use.**

- *Rung 0 (assay reliability)*: drop rows with no plate identifier, then work **within one
  (line, drug, dose)** — dose is part of the key, because pooling it would put different doses in
  the two halves for 99.7% of conditions. Split that triple's plates by `hash(plate) % 2`,
  average `log2FoldChange` per gene in each half, and correlate the halves, twice: over every
  gene, and over the genes the first half called differentially expressed by `padj`. Also, per
  (line, drug, dose, gene), the variance of `log2FoldChange` across plates against the mean
  `lfcSE²`, splitting replicate noise into its between-plate and within-plate parts. No drug
  filter and no gene panel: the 32-drug restriction above belongs to the rungs that need it.
- *Delta bundle* (arrives with rung 1; its download half is landed as
  `scripts/download_tahoe_pseudobulk_de.py`): the same aggregation without the split —
  mean `log2FoldChange` and `baseMean` per (line, drug, gene) over all plates and doses.
  Consumed by rung 1 ([docs/SPEC.md](SPEC.md)) when that rung arrives.
