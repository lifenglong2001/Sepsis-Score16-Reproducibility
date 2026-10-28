# Runtime

**This is a record of the environment the code was developed and executed in. It is not a dependency
lock, and it is not a verification certificate for every stage.**

## Host

| Item | Value |
|---|---|
| OS | Windows, x86_64 |
| R | 4.5.1 |
| PowerShell | 5.1 or later for `scripts/run release.ps1` |
| Python (figure renderers only) | 3.12 |

## Packages the stages load

| Package | Version | Needed by |
|---|---|---|
| data.table | 1.18.0 | every stage |
| metafor | 4.8-0 | `core`, `benchmarks` |
| jsonlite | 2.0.0 | `src/common.R` config parsing |
| GSVA | **2.4.9** | `pbmc_score`, `scoring_portable`, `pbmc_replication` |
| BiocParallel | 1.44.0 | `pbmc_score`, `scoring_portable`, `pbmc_replication` |
| edgeR | local build | `pbmc_replication` |
| glmnet | 5.0 | `cross_task` |
| randomForest | 4.7-1.2 | `cts_reconstruction` |
| R.utils | 2.13.0 | optional; only for `fread()` on `.gz` |

`GSVA` version is asserted in code (`stopifnot(packageVersion('GSVA') == '2.4.9')`). A different GSVA
changes ssGSEA values and the stage refuses to run.

`environment/R packages snapshot.tsv` is the broader observed library snapshot. Packages marked
"built under R 4.5.2 / 4.5.3" there load and run on R 4.5.1; R only warns about the build version.

## Requirements that affect whether a stage runs

- `R.utils` is **not** required: `src/common.R` falls back to a base-R `gzfile()` connection for `.gz`
  inputs. Without either, `fread()` fails with `fread() requires 'R.utils' package`.
- Windows PowerShell 5.1 turns native stderr into a terminating error under
  `$ErrorActionPreference='Stop'`; the launcher relaxes the preference around the native call and uses
  the exit code.
- A non-ASCII (Chinese) directory in an input path works if the config is read as UTF-8. R libraries on
  such a path do **not** work reliably (R resolves the real path and fails to convert the encoding
  under a C locale); keep R libraries on an ASCII path.
- There is no `renv.lock`. Producing one requires a clean-room library build.

## Not verified here

- Linux or macOS execution. The R code is platform-neutral and the invocation is documented, but only
  Windows was exercised.
- A clean-room dependency install.
- Whether every stage runs end to end on a fresh machine: that depends on the 40 processed inputs,
  which are not part of this package (see `Required input inventory.tsv`).
