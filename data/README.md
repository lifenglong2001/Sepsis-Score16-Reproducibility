# Data directory

**This directory intentionally contains no data.**

No patient-level table, expression matrix, single-cell object, deconvolution estimate or MSigDB
GMT is shipped with this repository. See `../DATA_ACCESS.md` for where each input comes from and
`../Required input inventory.tsv` for the full list of 40 configured input keys.

If you are adding local inputs for your own run, put them **outside** the repository, for example:

```
<your-input-directory>/
    GSE272769.tsv
    GSE65682.tsv
    ...
    h.all.v2026.1.Hs.symbols.gmt
```

and point `config/config.local.json` at them. Do not commit them.

`.gitignore` already excludes common local input and output patterns, but it is a safety net, not
a guarantee: check `git status` before every commit.
