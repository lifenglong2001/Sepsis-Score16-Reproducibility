#!/usr/bin/env Rscript
# Build TEST_RESULTS.tsv from real stage ledgers and assert it against the real
# process exit codes. No hand-written numbers: every count is re-derived by
# reading each ledger row by row.
#
# Usage:
#   Rscript tests/summarize_and_assert_ledgers.R --out TEST_RESULTS.tsv <run_dir> [<run_dir> ...]
#
# A <run_dir> is any directory produced by scripts/run release.ps1; it must
# contain results/*_comparison.tsv and may contain execution.json.
#
# Exit code 0 only if, for every ledger:
#   * the file parses,
#   * the `status` column exists and every value is PASS or FAIL,
#   * the number of FAIL rows is consistent with the recorded exit code
#     (exit 0 <=> zero FAIL rows),
#   * the ledger's own "CHECKS n FAIL m" line in logs/<stage>.log, when present,
#     agrees with the recount.
# Otherwise it exits non-zero and prints every contradiction.

suppressPackageStartupMessages({library(data.table);library(jsonlite)})

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || length(args) < 2L) {
  stop("usage: summarize_and_assert_ledgers.R --out <file.tsv> <run_dir> [<run_dir> ...]")
}
out_file <- NULL
if (args[1] == "--out") {
  out_file <- args[2]
  run_dirs <- args[-(1:2)]
} else {
  run_dirs <- args
}
if (!length(run_dirs)) stop("no run directories given")

rows <- list()
problems <- character(0)
add_problem <- function(...) problems <<- c(problems, paste0(...))

for (rd in run_dirs) {
  if (!dir.exists(rd)) { add_problem("missing run dir: ", rd); next }
  res_dir <- file.path(rd, "results")
  ledgers <- if (dir.exists(res_dir)) list.files(res_dir, pattern = "_comparison\\.tsv$", full.names = TRUE) else character(0)
  if (!length(ledgers)) { add_problem("no ledgers in ", rd); next }

  exec <- NULL
  ef <- file.path(rd, "execution.json")
  if (file.exists(ef)) {
    exec <- tryCatch({
      rawb <- readBin(ef, "raw", file.info(ef)$size)
      if (length(rawb) >= 3 && all(rawb[1:3] == as.raw(c(0xEF, 0xBB, 0xBF)))) rawb <- rawb[-(1:3)]
      txt <- rawToChar(rawb[rawb != as.raw(0)])
      j <- jsonlite::fromJSON(txt, simplifyDataFrame = TRUE)
      # A single-stage run writes a JSON object, not an array; wrap it.
      if (!is.data.frame(j) && is.list(j) && !is.null(j$stage)) j <- as.data.frame(j, stringsAsFactors = FALSE)
      j
    }, error = function(e) NULL)
  }

  for (lg in ledgers) {
    stage <- sub("_comparison\\.tsv$", "", basename(lg))
    d <- tryCatch(fread(lg), error = function(e) { add_problem("unreadable ledger: ", lg, " :: ", conditionMessage(e)); NULL })
    if (is.null(d)) next
    if (!"status" %in% names(d)) { add_problem("no status column: ", lg); next }
    bad_vals <- setdiff(unique(d$status), c("PASS", "FAIL"))
    if (length(bad_vals)) add_problem("unexpected status values in ", lg, ": ", paste(bad_vals, collapse = ","))
    n <- nrow(d)
    npass <- sum(d$status == "PASS")
    nfail <- sum(d$status == "FAIL")

    exit_code <- NA_integer_
    if (!is.null(exec) && all(c("stage", "exit_code") %in% names(exec))) {
      hit <- exec[exec$stage == stage, , drop = FALSE]
      if (nrow(hit) == 1L) exit_code <- as.integer(hit$exit_code)
    }
    if (!is.na(exit_code)) {
      if (exit_code == 0L && nfail > 0L)
        add_problem("CONTRADICTION ", stage, ": exit 0 but ", nfail, " FAIL rows")
      if (exit_code != 0L && nfail == 0L)
        add_problem("CONTRADICTION ", stage, ": exit ", exit_code, " but 0 FAIL rows")
    }

    logf <- file.path(rd, "logs", paste0(stage, ".log"))
    log_line <- ""
    if (file.exists(logf)) {
      # Stage logs are written by PowerShell redirection, which emits UTF-16LE
      # (with BOM) as soon as an error record is captured, and may also contain
      # localised non-UTF-8 text. Decode by BOM before matching, otherwise the
      # ASCII CHECKS line is invisible to grep.
      rawb <- readBin(logf, "raw", file.info(logf)$size)
      txt <- if (length(rawb) >= 2 && rawb[1] == as.raw(0xFF) && rawb[2] == as.raw(0xFE)) {
        iconv(list(rawb[-(1:2)]), from = "UTF-16LE", to = "UTF-8")
      } else if (length(rawb) >= 2 && rawb[1] == as.raw(0xFE) && rawb[2] == as.raw(0xFF)) {
        iconv(list(rawb[-(1:2)]), from = "UTF-16BE", to = "UTF-8")
      } else if (length(rawb) >= 3 && all(rawb[1:3] == as.raw(c(0xEF, 0xBB, 0xBF)))) {
        rawToChar(rawb[-(1:3)])
      } else {
        rawToChar(rawb[rawb != as.raw(0)])
      }
      if (is.na(txt)) txt <- ""
      ll <- strsplit(txt, "\r?\n", perl = TRUE)[[1]]
      # PowerShell redirection can leave a trailing space on each captured line,
      # so trim before matching the stage summary line.
      ll <- trimws(sub("\r$", "", ll))
      m <- grep(paste0("^", stage, " CHECKS [0-9]+ FAIL [0-9]+$"), ll, value = TRUE, useBytes = TRUE)
      if (length(m)) {
        log_line <- tail(m, 1)
        ln <- as.integer(sub(".*CHECKS ([0-9]+) FAIL ([0-9]+)$", "\\1", log_line))
        lf <- as.integer(sub(".*CHECKS ([0-9]+) FAIL ([0-9]+)$", "\\2", log_line))
        if (!is.na(ln) && ln != n) add_problem("CONTRADICTION ", stage, ": log CHECKS ", ln, " != ledger rows ", n)
        if (!is.na(lf) && lf != nfail) add_problem("CONTRADICTION ", stage, ": log FAIL ", lf, " != ledger FAIL ", nfail)
      } else {
        add_problem("no CHECKS line found in ", logf)
      }
    } else {
      add_problem("no stage log for ", stage, " in ", rd)
    }

    rows[[length(rows) + 1L]] <- data.table(
      stage = stage, run_dir = rd, checks = n, pass = npass, fail = nfail,
      exit_code = exit_code, ledger = file.path("results", basename(lg)),
      log_summary = log_line, status = if (nfail == 0L) "PASS" else "FAIL"
    )
  }
}

res <- rbindlist(rows, fill = TRUE)
if (!is.null(out_file)) {
  fwrite(res, out_file, sep = "\t", na = "NA")
  cat("wrote", out_file, "\n")
}
print(res[, .(stage, run_dir, checks, pass, fail, exit_code, status)])

cat("\nTOTAL rows:", sum(res$checks), " PASS:", sum(res$pass), " FAIL:", sum(res$fail), "\n")
if (length(problems)) {
  cat("\nLEDGER/SUMMARY ASSERTION FAILED:\n")
  cat(paste0("  - ", problems, collapse = "\n"), "\n")
  quit(status = 1L)
}
cat("LEDGER/SUMMARY ASSERTION: OK (summary matches every ledger and every exit code)\n")
quit(status = 0L)
