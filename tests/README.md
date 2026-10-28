# Tests

Three suites. None of them is a substitute for a real run.

| Test | What it checks | Runtime |
|---|---|---|
| `test_static_integrity.R` | every required file exists, every R source parses, no absolute local paths, no credentials, all text is valid UTF-8, the config still uses `INPUT_ROOT` placeholders for all 40 inputs, the 12 stages are wired into both entry points, and no stale stage/check counts appear in the reader-facing documents | seconds |
| `test_checks_integrity.R` | every `checks/*.csv` loads, the 14+2 membership is unchanged, and the headline frozen values are still the ones the paper reports | seconds |
| `summarize_and_assert_ledgers.R` | rebuilds `TEST_RESULTS.tsv` from the real stage ledgers and asserts it against the real exit codes and stage logs | seconds |

Run the first two:

```powershell
.\tests\run tests.ps1
```

Run the summary/assertion against one or more run directories:

```powershell
Rscript tests/summarize_and_assert_ledgers.R --out TEST_RESULTS.tsv <run-directory>
Rscript tests/summarize_and_assert_ledgers.R --out TEST_RESULTS.tsv <run-directory-1> <run-directory-2>
```

All three exit non-zero on the first failure.

## Path requirement on Windows

R resolves the package root from the script's own path. If the package sits under a directory whose
name contains non-ASCII characters, R may fail before any test runs, with
`unable to translate '...' to UTF-8` or `file name conversion problem -- name too long?`.

**Workaround: copy the package to a path that is pure ASCII and run there.** If the script path cannot
be read but you know where the package is, set `SEPSIS_CODE_ROOT` to that ASCII location:

```powershell
$env:SEPSIS_CODE_ROOT = "<your-ASCII-path>\Sepsis-Score16-Reproducibility"
Rscript "<your-ASCII-path>\Sepsis-Score16-Reproducibility\tests\test_static_integrity.R"
```

This affects the test harness only; the analysis stages are unaffected because they take absolute
input paths from the config.

## Scope

`test_checks_integrity.R` fails if any headline value in `checks/` changes, before a stage is run.
`summarize_and_assert_ledgers.R` reads each ledger with a real TSV parser and cross-checks the row
count against the `execution.json` exit code and the `CHECKS n FAIL m` line in the stage log.
