# Sepsis Score16 — reproducibility code

Analysis code, configuration template and frozen numerical comparison targets for the study
**"A blood transcriptomic pathway signature associated with mortality in sepsis"**.

Score16 is a fixed 16-member Hallmark transcriptomic score:

```
Score16 = mean(14 positive Hallmarks) − mean(2 negative Hallmarks)
```

The 14 positive and 2 negative members, and their directions, are fixed. They were selected once in
four adult discovery cohorts and are never re-selected, re-weighted or sign-flipped.

---

## 1. Cohorts

| Role | Cohort | Accession | n | Deaths |
|---|---|---|---|---|
| Discovery | GSE272769 | GEO | 161 | 60 |
| Discovery | GSE65682 (MARS) | GEO | 478 | 114 |
| Discovery | GSE185263 | GEO | 342 | 50 |
| Discovery | GSE236713 | GEO | 124 | 30 |
| External | GSE63042 | GEO | 106 | 28 |
| External | GSE95233 | GEO | 51 | 17 |
| External | VANISH (E-MTAB-7581) | ArrayExpress | 176 | 48 |
| PBMC, donor-level | GSE279452 | GEO | 108 donors | — |
| PBMC replication | GSE290679 | GEO | 19 sepsis ICU + 19 non-septic ICU + 9 healthy donors | — |
| Neutrophil states | GSE216009 | GEO | 39 acute donors | — |
| Ordered severity | GSE154918 | GEO | 91 | — |

Frozen external benchmark: **OR 1.879 (95% CI 1.209–2.920)** across the three external cohorts
(REML with HKSJ adjustment).

---

## 2. What the package reproduces

Each stage reads patient- or donor-level **processed** inputs, recomputes the quantities reported in
the paper, and compares them against the frozen targets in `checks/`. A stage exits non-zero if any
comparison fails. Bootstrap interval endpoints are handled as described in §7.

| Stage | Checks | What it recomputes |
|---|---:|---|
| `core` | 1934 | 7 cohorts; fixed 16 members; 4 × 50 Hallmark univariable models; REML meta with HKSJ and modified-KH; LOCO/LOSO; cell-cycle collapse; QC exclusion |
| `pbmc` | 1445 | donor-level Wilcoxon and Cliff's δ, 80 pathways × 2 contrasts, BH families |
| `pbmc_replication` | 658 | GSE290679 donor-level replication: 50 effect rows plus the 14/25 and 0/25 summary |
| `clinical_task` | 550 | SOFA-50 (3 model forms) and ordered-severity Spearman |
| `cellular` | 352 | 17-donor maturity variance; deconvolution correlations and conditioned models |
| `neutrophil_paired` | 306 | GSE216009 paired immature vs mature, 51 metrics, 17 donors |
| `benchmarks` | 195 | CTS / SRS / Hi-DEF / MARS, 2000 paired stratified bootstraps, apparent AUC |
| `granulopoiesis` | 104 | five-gene module (AZU1/CTSG/ELANE/MPO/PRTN3): correlation with Score16 in the seven adult cohorts, MARS-group comparison in GSE65682, immature-vs-mature localization in GSE216009 |
| `pbmc_score` | 80 | GSVA 2.4.9 ssGSEA recomputed from donor logCPM |
| `scoring_portable` | 63 | ssGSEA re-derived for 7 cohorts; rank-score portability |
| `cross_task` | 28 | diagnostic-only nested CV (20×5 outer, 5 inner); unchanged OOF evaluated for mortality |
| `cts_reconstruction` | 14 | author randomForest workflow; 116/117 published CTS labels |

**Total: 5729 compared quantities.**

### Run requirements and scope

- **40 processed inputs are required and are not bundled.** `Required input inventory.tsv` lists each
  one with its official source and the preprocessing step behind it.
- **No analysis matrices are bundled or downloaded automatically.** The optional helper retrieves GEO
  series metadata and the pinned CTS reference files only.
- **Not a raw-to-paper rebuild.** Normalisation, probe collapse, QC exclusion maps and the exact
  clinical sample-ID joins are outside the package. The starting point is a per-cohort table of 50
  Hallmark scores plus the expression matrices.
- **Table T7 is only partly recomputed.** The 51-row paired immature-vs-mature block is recomputed;
  the 714 unpaired cell-subtype rows, the 48 paired subtype-vs-mature rows and the frozen
  `fdr_global` / `fdr_within_metric` families are frozen targets only.

---

## 3. Requirements

- **R 4.5.1**; observed package versions in `environment/R packages snapshot.tsv`
- Always required: `data.table`, `metafor`, `jsonlite`
- Stage-specific: `GSVA` **exactly 2.4.9**, `BiocParallel`, `glmnet`, `randomForest`, `edgeR`, and
  `R.utils` only if you want `fread()` to decompress `.gz` natively
- PowerShell 5.1+ for the Windows entry point; any shell plus `Rscript` elsewhere
- Python 3.12 with NumPy / pandas / Matplotlib and Arial, only for the figure renderers

There is no automatic installer. Point `r_libraries` in your config at your own library paths.

---

## 4. Running

```powershell
# 1. copy the example config and replace every INPUT_ROOT placeholder
Copy-Item config/config.example.json config/config.local.json
#    - inputs.*      -> your processed inputs (40 keys)
#    - output        -> a new, empty directory
#    - r_libraries   -> your R library paths (optional)

# 2. run a fresh analysis into a new empty output directory
.\scripts\run release.ps1 `
    -Config config\config.local.json `
    -Rscript "<path-to-Rscript>" `
    -OutputDirectory "<your-new-empty-directory>" `
    -Stages core,benchmarks,pbmc,pbmc_score,cellular,cross_task,clinical_task,cts_reconstruction,scoring_portable,pbmc_replication,neutrophil_paired,granulopoiesis

# 3. rebuild the summary from the real ledgers and assert it
Rscript tests/summarize_and_assert_ledgers.R --out TEST_RESULTS.tsv "<your-new-empty-directory>"
```

Linux / macOS:

```bash
export SEPSIS_CODE_ROOT="$PWD"
Rscript "scripts/run release.R" config/config.local.json core
```

Each stage writes into `<output>/results/` the `*_comparison.tsv` ledger, the analysis tables and a
`*_sessionInfo.txt`; into `<output>/logs/` the stage log; and `<output>/execution.json` with exit
codes and timings.

---

## 5. Figures and tables

`RESULTS_TO_FIGURES.md` maps every main figure, supplementary figure and table to the stage and the
files behind it.

The repository contains drawing scripts for **Figure 1, Figure 3, Figure S3, Figure S4 and
Figure S10**:

| Script | Output |
|---|---|
| `figures/Figure1_S4/draw_figure1_s4.py` | Figure 1, Figure S4 |
| `figures/Figure3_S10/draw_figure3_s10.py` | Figure 3, Figure S10 |
| `figures/Figure_S3/draw_figure_s3.py` | Figure S3 |

Each reads bundled frozen aggregates and redraws them; none refits a model or reruns a bootstrap.
The remaining artefacts are not generated by scripts in this repository.

---

## 6. Multiple-testing families

P values are adjusted within defined testing families using Benjamini–Hochberg.

- **GSE290679 replication (Table T9).** For each comparator, BH correction was applied across the
  25 prespecified lineage–pathway comparisons (5 metabolic/stress Hallmarks × 5 PBMC lineages).
  The frozen T9 table is verified for internal consistency (row counts, the 5 × 5 grid,
  within-comparator BH, the support rule, and the 14/25 and 0/25 counts). That verification reads
  the frozen table only: it is **not** a re-run of the `pbmc_replication` stage, which requires
  donor-level inputs not bundled here.
- **T7.** The frozen `fdr_global` and `fdr_within_metric` families are not re-derived here.

---

## 7. Limitations

1. **Processed inputs only.** The upstream steps listed in `DATA_ACCESS.md` are not included.
2. **MSigDB GMT not bundled** for licensing reasons; the exact release string is pinned.
3. **Single-cell object construction not included.** Five stages start from frozen donor-level objects.
4. **T7 partial** — see §2.
5. **Bootstrap intervals are recomputed; their endpoints are not asserted.** The `granulopoiesis`
   stage recomputes bootstrap percentile intervals from a single stage seed and writes them to the
   stage outputs. Point estimates and non-interval quantities are asserted against `checks/`. A
   subset of interval endpoints has not been shown to reproduce the manuscript values exactly in
   this environment, so interval endpoints are not asserted. The endpoints in `checks/` are the
   manuscript-reported targets.
6. **No dependency lock file.** `environment/R packages snapshot.tsv` is an observed snapshot.
7. **Windows-only execution was tested.** The R code is platform-neutral and the Linux invocation is
   documented, but only Windows with R 4.5.1 was exercised.

---

## 8. Licence and citation

**Licence: not yet chosen.** No reuse right is granted by this repository until the authors add a
licence file. Do not redistribute.

`CITATION.cff` lists the manuscript authors. ORCID iDs, version, release date and DOI will be added
when they exist.
