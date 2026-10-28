# Figure 3 and Figure S10

Redraws Figure 3 (four panels) and Figure S10 from frozen aggregates. No model is fitted, no bootstrap
is rerun and no pathway is re-selected.

## Requirements

Python 3.12 with NumPy, pandas and Matplotlib, and the Arial font installed.

## Run

From this directory:

```text
python draw_figure3_s10.py --output-dir NEW_OUTPUT_DIRECTORY
```

The script reads `source_data/` relative to itself and writes the figure files into the directory you
pass.

## Notes

- Figure 3 covers adults only. Its panel d is a point-only AUC comparison.
- Figure S10 uses the unchanged B = 1,000 portability outputs.
- The rendered panels label the fixed score as **Score16**.
