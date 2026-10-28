# Figure S3

Redraws the three-panel Figure S3 from frozen aggregates. No model is refitted, no bootstrap is rerun
and no pathway is re-selected.

## Requirements

Python 3.12 with NumPy, pandas and Matplotlib, and the Arial font installed.

## Run

From this directory:

```text
python draw_figure_s3.py --output-dir NEW_OUTPUT_DIRECTORY
```

The script locates `source_data/` relative to itself.

## Inputs

| File | Used for |
|---|---|
| `source_data/discovery fixedscore meta.csv` | primary fixed Score16 and the four fixed-score discovery omissions (panel a) |
| `source_data/discovery cohort effects.csv` | primary cohort effects; reference points in panels b and c |
| `source_data/cell cycle collapsed.csv` | cell-cycle-collapsed cohort effects (panel b) |
| `source_data/qc exclusion.csv` | predefined QC-exclusion cohort effects (panel c) |
| `source_data/cell cycle collapsed meta.csv` | pooled collapsed-score effect, read and checked for identity and CI validity; not plotted |

## Outputs

`Figure_S3.pdf`, `Figure_S3.png` (300 dpi), `Figure_S3.svg` and `Figure_S3_render_QC.json`. The JSON
records the input hashes, the omission odds ratios, clipping checks and a no-fitting flag.
