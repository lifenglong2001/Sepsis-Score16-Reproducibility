# scripts/upstream — GSE290679 pseudobulk construction

**These files are reference documentation of the upstream step. They are NOT run by
`scripts/run release.ps1` and are NOT part of any analysis stage.**

| File | Purpose |
|---|---|
| `pseudobulk_from_h5ad.py` | aggregates the author-annotated GSE290679 `.h5ad` into a donor × lineage pseudobulk count matrix plus donor metadata. Requires `h5py`, `numpy`, `scipy` and the ~10.9 GB raw file. |

## Why they are here

The release starts at the donor × lineage pseudobulk matrix, because that is the smallest
donor-level input that supports a real recomputation of Table T9. Shipping the construction script
makes the boundary explicit and auditable instead of hiding it.

## What `pseudobulk_from_h5ad.py` does

1. reads `obs` from the h5ad and maps the author's fine cell labels onto the five frozen lineages
   (Monocyte, CD4 T, CD8 T, B, NK). Plasmablasts are deliberately **not** folded into B, matching
   the GSE279452 convention;
2. asserts the donor design: 19 `CI-Sep`, 19 `CI-NS`, 9 `NHC`, 47 unique donor IDs, all adults;
3. sums raw counts per donor × lineage;
4. writes `GSE290679_pseudobulk_counts.tsv.gz` and `GSE290679_pseudobulk_metadata.csv`.

## Differences from the shipped stage

The shipped `src/pbmc replication.R` reads the resulting pseudobulk matrix and does the
normalisation (TMM + logCPM), ssGSEA scoring, Cliff's δ, bootstraps, BH and the frozen-direction
comparison. It does not touch the h5ad and does not re-run the aggregation.

If you regenerate the pseudobulk matrix yourself, verify that `pbmc_replication` still reproduces
Table T9 (14/25 against healthy controls, 0/25 against non-septic ICU). A different lineage map or
a different cell-count floor will break that comparison — which is the point of checking it.
