# Data access and input preparation

This repository contains **no patient-level or donor-level data**. It contains code, configuration,
frozen comparison targets and figure renderers. You must obtain the inputs yourself.

`Required input inventory.tsv` lists all **40** configured inputs with their official source, recorded
version and the preprocessing step that produced them. The configuration keeps the literal
`INPUT_ROOT` token as a placeholder for your own input directory.

## Public data sources

| Key | Source | Notes |
|---|---|---|
| GSE272769 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272769 | discovery; processed matrix + sample metadata |
| GSE65682 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65682 | discovery; MARS cohort |
| GSE185263 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185263 | discovery |
| GSE236713 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE236713 | discovery |
| GSE63042 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63042 | external |
| GSE95233 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95233 | external |
| VANISH | https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-7581 | external; ArrayExpress E-MTAB-7581 |
| GSE154918 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154918 | ordered-severity task |
| GSE216009 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216009 | neutrophil maturity, paired immature/mature pseudobulk, module localization |
| GSE279452 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279452 | PBMC donor-level primary analysis |
| GSE290679 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290679 | independent PBMC replication cohort (Figure 6e, Figure S9, Table T9); 19 sepsis ICU, 19 non-septic ICU, 9 healthy donors |
| CTS reference files | https://github.com/bpscicluna/Consensus-Transcriptomic-Subtypes-of-Sepsis/tree/467fa897ba33b8ec355a8f33a1d3c4f4486846bd | pinned commit, three files |
| MSigDB Hallmark | https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#H | `h.all.v2026.1.Hs.symbols.gmt`, exact release |

## What each stage expects

- **Processed patient-level tables** (`GSE272769.tsv`, …, `VANISH.tsv`) — one row per patient with
  `sample_id`, `event_death` and exactly 50 `HALLMARK_*` score columns. These define the analysis
  subset; the expression matrices may contain more samples than the analysis table.
- **Normalised expression matrices** — one per cohort, gene symbols as row names, sample accessions as
  column names. Used for the within-sample rank score and the portable-scoring stage.
- **Single-cell and deconvolution intermediates** — frozen donor-level objects for `pbmc_score`,
  `cellular`, `neutrophil_paired` and `granulopoiesis`.
- **GSE290679 donor × lineage pseudobulk counts, donor metadata and the predefined direction
  registry** — for `pbmc_replication`.
- **Frozen subtype label joins** — VANISH SRS/Hi-DEF and MARS.
- **MSigDB Hallmark GMT** — the exact release, obtained from MSigDB under its terms.

## Optional acquisition helper

`scripts/acquire public prerequisites.py` retrieves GEO SOFT series metadata and the three pinned CTS
author files. It does **not** download analysis matrices and does not run any analysis:

```powershell
python "scripts/acquire public prerequisites.py" --destination "<your-input-directory>" --plan
python "scripts/acquire public prerequisites.py" --destination "<your-input-directory>" --cts
python "scripts/acquire public prerequisites.py" --destination "<your-input-directory>" --geo-metadata
```

The CTS command fetches three unchanged author reference/target expression files at the pinned commit;
it does not train a model. GEO metadata acquisition retrieves series SOFT metadata, not normalised
matrices or reads. The original author model reference has 460 samples and 18 features; do not
substitute the separate 421-sample package reference.

## GSE290679: the replication cohort

GSE290679 is a separate public scRNA-seq cohort used only as an independent PBMC replication
comparator. It is **not** part of the Score16 discovery or external mortality set, and it does not
contribute to any Score16 effect estimate.

The package starts from a **donor × lineage pseudobulk count matrix**, not from the 10.9 GB
author-annotated `.h5ad`. Regenerating that matrix from the h5ad is possible and the script is included
as `scripts/upstream/pseudobulk_from_h5ad.py` (requires `h5py`, `numpy`, `scipy` and the raw file); it
is **not** part of any analysis stage and is not run by the driver. The lineage map is frozen
(Monocyte / CD4 T / CD8 T / B / NK) and plasmablasts are deliberately kept out of B, matching the
GSE279452 convention.

```bash
python scripts/upstream/pseudobulk_from_h5ad.py \
    --input-h5ad path/to/GSE290679_raw_labeled_filtered.h5ad \
    --output-dir path/to/output
```

This writes `GSE290679_pseudobulk_counts.tsv.gz` and `GSE290679_pseudobulk_metadata.csv` in the
output directory.

| Key | Content |
|---|---|
| `gse290679_pseudobulk_counts` | gene × donor-lineage pseudobulk integer counts |
| `gse290679_pseudobulk_metadata` | donor, lineage, group (`CI-Sep` / `CI-NS` / `NHC`), cell count |
| `gse290679_direction_registry` | frozen GSE279452 Cliff's deltas used as the replication target |

## MSigDB licensing

The Hallmark GMT is distributed by MSigDB under its own terms. It is **not** bundled here. Obtain the
exact release yourself and keep it locally. Do not commit it to this repository.

## Upstream steps that are not included

- Raw sequencing/array intensities (FASTQ, CEL, IDAT) and raw single-cell `.h5ad`.
- Original normalisation / probe-to-gene collapse wrappers and the frozen expression builds.
- QC exclusion maps and the exact clinical sample-ID joins.
- Single-cell object construction and MuSiC/DWLS projection code; the frozen patient-level estimates
  are used instead.
- The published-117 CTS exact sample-ID join.

These are named per input in `Required input inventory.tsv`. They are the reason this package is a
processed-input reproduction rather than a raw-to-paper rebuild.

## Preparing a run

1. Download the processed matrices for the seven adult cohorts and the single-cell objects you need.
2. Compute 50 Hallmark ssGSEA scores per cohort under GSVA 2.4.9 with the pinned GMT, producing one
   `<cohort>.tsv` per cohort with `patient_id`, `sample_id`, `event_death` and 50 `HALLMARK_*`.
3. Obtain or regenerate the GSE290679 donor × lineage pseudobulk matrix and metadata.
4. Fill in `config/config.local.json` and run the stages (see `README.md`).
