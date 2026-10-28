#!/usr/bin/env Rscript
# Static integrity tests. No data needed. Exit code 0 = all pass.
# Resolves the package root from this file's own location, so it can be run from anywhere.
# Set SEPSIS_CODE_ROOT to override, or run from the package root if the path contains
# characters R cannot translate on this host (see the fallback below).

args <- commandArgs(trailingOnly = FALSE)
this <- sub("^--file=", "", args[grep("^--file=", args)])
root <- Sys.getenv("SEPSIS_CODE_ROOT", "")
if (!nzchar(root)) {
  root <- tryCatch(normalizePath(file.path(dirname(this), ".."), winslash = "/", mustWork = TRUE),
                   error = function(e) NA_character_)
}
if (is.na(root) || !nzchar(root)) {
  # R on Windows cannot always convert a non-ASCII script path to UTF-8. Fall back to the
  # working directory, but only accept it if it really looks like the package root.
  cand <- tryCatch(normalizePath(getwd(), winslash = "/", mustWork = TRUE), error = function(e) NA_character_)
  if (is.na(cand) || !file.exists(file.path(cand, "src", "common.R"))) {
    cat("FATAL: cannot resolve the package root from the script path or the working directory.\n")
    cat("       Run from the package root, or set SEPSIS_CODE_ROOT.\n")
    quit(status = 2L)
  }
  root <- cand
  cat("NOTE: resolved the package root from the working directory, not the script path.\n")
}

fails <- 0L
ok <- function(name, cond, detail = "") {
  if (isTRUE(cond)) {
    cat("PASS", name, "\n")
  } else {
    fails <<- fails + 1L
    cat("FAIL", name, ":", detail, "\n")
  }
}

# 1. required top-level files
required <- c(
  "README.md", "RESULTS_TO_FIGURES.md", "DATA_ACCESS.md", "CITATION.cff", ".gitignore",
  "config/config.example.json", "Required input inventory.tsv",
  "scripts/run release.ps1", "scripts/run release.R", "scripts/verify release.ps1",
  "scripts/acquire public prerequisites.py", "scripts/upstream/README.md",
  "tests/summarize_and_assert_ledgers.R", "tests/README.md",
  "environment/runtime.txt", "environment/R packages snapshot.tsv", "environment/RUNTIME.md",
  "data/README.md"
)
for (f in required) ok(paste0("exists/", f), file.exists(file.path(root, f)), f)

# 2. every R source parses
r_files <- c(list.files(file.path(root, "src"), pattern = "\\.R$", full.names = TRUE),
             list.files(file.path(root, "tests"), pattern = "\\.R$", full.names = TRUE),
             file.path(root, "scripts", "run release.R"))
r_files <- r_files[file.exists(r_files)]
# src/ holds the 12 analysis stages plus the shared helper common.R
ok("src_has_12_stages_plus_common", length(list.files(file.path(root, "src"), pattern = "\\.R$")) == 13L,
   length(list.files(file.path(root, "src"), pattern = "\\.R$")))
for (f in r_files) {
  res <- tryCatch({ parse(f, encoding = "UTF-8"); TRUE }, error = function(e) conditionMessage(e))
  ok(paste0("parse/", basename(f)), isTRUE(res), if (isTRUE(res)) "" else res)
}

# 3. no absolute local paths in shipped analysis code
# Scan shipped code only: the test files themselves legitimately contain these
# patterns, so including them would be a self-match.
# Keep the pattern pure ASCII: R running in the C locale cannot translate a
# non-ASCII regex, and this test must work on any host.
shipped <- c(list.files(file.path(root, "src"), pattern = "\\.R$", full.names = TRUE),
             file.path(root, "scripts", "run release.ps1"),
             file.path(root, "scripts", "verify release.ps1"),
             file.path(root, "scripts", "run release.R"),
             file.path(root, "scripts", "acquire public prerequisites.py"),
             file.path(root, "config", "config.example.json"))
hits <- character(0)
for (f in shipped) {
  if (!file.exists(f)) next
  ln <- readLines(f, warn = FALSE, encoding = "UTF-8")
  # A drive-letter path, not a URL scheme: the colon must not be preceded by a
  # letter or digit, so "https://" does not match but "C:/Users" does.
  m <- grep("(^|[^A-Za-z0-9])[A-Za-z]:[\\\\/]|/Users/|/home/", ln, value = TRUE)
  if (length(m)) hits <- c(hits, paste0(basename(f), ": ", trimws(head(m, 1))))
}
ok("no_absolute_local_paths", !length(hits), paste(hits, collapse = " | "))

# 4. no credentials in shipped code
cred <- character(0)
for (f in shipped) {
  if (!file.exists(f)) next
  ln <- readLines(f, warn = FALSE, encoding = "UTF-8")
  m <- grep("(api[_-]?key|secret|password|passwd|token\\s*=|BEGIN [A-Z ]*PRIVATE KEY)", ln,
            ignore.case = TRUE, value = TRUE)
  if (length(m)) cred <- c(cred, basename(f))
}
ok("no_credentials_in_shipped_code", !length(cred), paste(unique(cred), collapse = ", "))

# 5. UTF-8 validity of shipped text
txt <- c(r_files, file.path(root, "README.md"), file.path(root, "config", "config.example.json"))
bad <- character(0)
for (f in txt) {
  if (!file.exists(f)) next
  raw <- readBin(f, "raw", file.size(f))
  if (length(raw) && !validUTF8(rawToChar(raw))) bad <- c(bad, basename(f))
}
ok("utf8_valid", !length(bad), paste(bad, collapse = ", "))

# 6. config still has INPUT_ROOT placeholders (i.e. no personal path leaked in)
cfg <- paste(readLines(file.path(root, "config", "config.example.json"), warn = FALSE), collapse = "\n")
ok("config_uses_placeholders", grepl("INPUT_ROOT/", cfg, fixed = TRUE), "config.example.json")
ok("config_has_40_inputs",
   length(gregexpr('"INPUT_ROOT/', cfg, fixed = TRUE)[[1]]) == 40L,
   length(gregexpr('"INPUT_ROOT/', cfg, fixed = TRUE)[[1]]))

# 7. the twelve-stage stage list is wired end to end
rdrv <- paste(readLines(file.path(root, "scripts", "run release.R"), warn = FALSE), collapse = "\n")
psdrv <- paste(readLines(file.path(root, "scripts", "run release.ps1"), warn = FALSE), collapse = "\n")
srcfiles <- list.files(file.path(root, "src"), pattern = "\\.R$")
stage_files <- setdiff(sub("\\.R$", "", srcfiles), "common")
ok("granulopoiesis_stage_source", "granulopoiesis.R" %in% srcfiles, paste(srcfiles, collapse = ","))
ok("granulopoiesis_in_r_driver", grepl("'granulopoiesis'", rdrv, fixed = TRUE), "run release.R")
ok("granulopoiesis_in_ps_driver", grepl("'granulopoiesis'", psdrv, fixed = TRUE), "run release.ps1")
ok("granulopoiesis_inputs_wired",
   grepl("'gse216009_pseudobulk_expression','mars_joined'", psdrv, fixed = TRUE),
   "run release.ps1 input list")

# 8. no stale stage/check counts in reader-facing documents
# The release ships 12 stages and 5729 compared quantities. Older counts must not reappear as the
# current claim anywhere a reader can see them.
stale <- c("4661", "5,625", "5625", "11 250", "11250", "eleven stages", "All 11 stages",
           "All 11 analysis stages", "11 / 11", "14 cohorts", "14-cohort", "5971", "5731")
docs <- c("README.md", "DATA_ACCESS.md", "RESULTS_TO_FIGURES.md", "tests/README.md")
hits <- character(0)
for (f in docs) {
  p <- file.path(root, f)
  if (!file.exists(p)) { hits <- c(hits, paste0(f, ": missing")); next }
  txt <- paste(readLines(p, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  for (s in stale) {
    if (grepl(s, txt, fixed = TRUE)) hits <- c(hits, paste0(f, " :: ", s))
  }
}
ok("no_stale_release_counts", !length(hits), paste(hits, collapse = " | "))

# 9. the current counts are actually stated
ok("readme_states_5729", grepl("5729", paste(readLines(file.path(root, "README.md"), warn = FALSE), collapse = "\n"), fixed = TRUE), "README.md")
ok("readme_states_12_stages", grepl("granulopoiesis", paste(readLines(file.path(root, "README.md"), warn = FALSE), collapse = "\n"), fixed = TRUE), "README.md")

cat("\nSTATIC INTEGRITY:", if (fails == 0L) "ALL PASS" else paste0(fails, " FAILURES"), "\n")
quit(status = if (fails == 0L) 0L else 1L)
