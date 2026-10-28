#!/usr/bin/env Rscript
# Frozen-target integrity tests. No patient data needed: only checks/*.csv.
# These guard the values the paper reports, so that a drifting target cannot
# make a failing analysis look like a pass. Exit code 0 = all pass.

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
  if (is.na(cand) || !file.exists(file.path(cand, "checks"))) {
    cat("FATAL: cannot resolve the package root from the script path or the working directory.\n")
    cat("       Run from the package root, or set SEPSIS_CODE_ROOT.\n")
    quit(status = 2L)
  }
  root <- cand
  cat("NOTE: resolved the package root from the working directory, not the script path.\n")
}
suppressPackageStartupMessages(library(data.table))

fails <- 0L
ok <- function(name, cond, detail = "") {
  if (isTRUE(cond)) cat("PASS", name, "\n") else {
    fails <<- fails + 1L
    cat("FAIL", name, ":", detail, "\n")
  }
}
near <- function(a, b, tol = 1e-9) isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol))
C <- function(f) fread(file.path(root, "checks", f))

# ---- 1. every check table loads -------------------------------------------------
files <- list.files(file.path(root, "checks"), pattern = "\\.csv$")
ok("checks_present", length(files) >= 20L, length(files))
for (f in files) {
  r <- tryCatch({ d <- C(f); nrow(d) > 0L }, error = function(e) conditionMessage(e))
  ok(paste0("loads/", f), isTRUE(r), if (isTRUE(r)) "" else r)
}

# ---- 2. fixed 14 + 2 membership -------------------------------------------------
m <- C("T2 Score16 membership.csv")
pos <- m[direction == "positive", pathway]
neg <- m[direction == "negative", pathway]
ok("membership_total_16", nrow(m) == 16L, nrow(m))
ok("membership_14_positive", length(pos) == 14L, length(pos))
ok("membership_2_negative", length(neg) == 2L, length(neg))
ok("membership_unique", !anyDuplicated(m$pathway), anyDuplicated(m$pathway))
ok("membership_negative_is_notch_wnt",
   setequal(neg, c("HALLMARK_NOTCH_SIGNALING", "HALLMARK_WNT_BETA_CATENIN_SIGNALING")),
   paste(neg, collapse = ","))
ok("membership_excludes_E2F", !("HALLMARK_E2F_TARGETS" %in% m$pathway), "E2F present")

# ---- 3. headline frozen values --------------------------------------------------
e <- C("T5 External meta LOSO.csv")[analysis == "primary_three_units"]
ok("external_meta_rows", nrow(e) == 1L, nrow(e))
ok("external_meta_k3", e$k == 3L, e$k)
ok("external_meta_n333", e$n == 333L, e$n)
ok("external_meta_deaths93", e$deaths == 93L, e$deaths)
ok("external_OR_1.879", near(e$OR, 1.87936767104486), e$OR)
ok("external_CI_low_1.209", near(e$CI_low, 1.20948912829893), e$CI_low)
ok("external_CI_high_2.920", near(e$CI_high, 2.92026010017648), e$CI_high)

d <- C("T4 Score16 discovery.csv")
ok("discovery_4_cohorts", nrow(d) == 4L, nrow(d))
ok("discovery_n_1105", sum(d$n) == 1105L, sum(d$n))
ok("discovery_deaths_254", sum(d$deaths) == 254L, sum(d$deaths))

v <- C("T12 VANISH CTS SRS HiDEF.csv")[analysis == "CTS_unadjusted"]
ok("vanish_CTS_deltaAUC_0.103", near(v$delta_AUC, 0.102783203125), v$delta_AUC)
ok("vanish_CTS_n176", v$n == 176L, v$n)
ok("vanish_CTS_B2000", v$B == 2000L, v$B)
ok("vanish_CTS_no_bootstrap_failures", v$bootstrap_failures == 0L, v$bootstrap_failures)

t10 <- C("T10 OOF diagnosis.csv")
# T10 carries four rows: diagnostic OOF -> diagnosis for two cohorts, plus the
# same unchanged OOF score evaluated for mortality in the same two cohorts.
ok("T10_four_rows", nrow(t10) == 4L, nrow(t10))
diag10 <- t10[grepl("-> diagnosis", analysis, fixed = TRUE)]
ok("T10_two_diagnosis_rows", nrow(diag10) == 2L, nrow(diag10))
ok("T10_GSE236713_AUC_0.8988",
   near(diag10[cohort == "GSE236713"]$AUC, 0.898769100169779),
   diag10[cohort == "GSE236713"]$AUC)
ok("T10_GSE63042_AUC_0.6317",
   near(diag10[cohort == "GSE63042"]$AUC, 0.631665299425759),
   diag10[cohort == "GSE63042"]$AUC)

p <- C("T13 Portable scoring.csv")
ok("T13_three_external", nrow(p) == 3L, nrow(p))
ok("T13_bootstrap_1000", all(p$AUC_CI_low <= p$AUC & p$AUC <= p$AUC_CI_high), "CI ordering")

c6 <- C("T6 Accepted CTS counts.csv")
ok("T6_three_classes", nrow(c6) == 3L, nrow(c6))
ok("T6_counts_sum_176", sum(c6$N) == 176L, sum(c6$N))

# ---- 4. no NaN / Inf in frozen numeric columns ----------------------------------
num_bad <- character(0)
for (f in c("T5 External cohorts.csv", "T4 Score16 discovery.csv",
            "T12 VANISH CTS SRS HiDEF.csv", "T13 Portable scoring.csv")) {
  dd <- C(f)
  for (nm in names(dd)) {
    if (is.numeric(dd[[nm]]) && any(is.infinite(dd[[nm]]))) num_bad <- c(num_bad, paste0(f, ":", nm))
  }
}
ok("no_infinite_frozen_values", !length(num_bad), paste(num_bad, collapse = ", "))

# ---- 5. T9 PBMC replication: testing family, support rule and support counts ---
# Scope: these checks validate the frozen T9 table and its internal consistency only.
# They do NOT re-run the pbmc_replication stage from its inputs; that stage is not
# executed by this test (see the run log for the stages actually executed).
t9 <- C("T9 PBMC replication.csv")
ok("T9_rowcount_50", nrow(t9) == 50L, nrow(t9))

T9_HALLMARKS <- c("HALLMARK_HYPOXIA", "HALLMARK_MTORC1_SIGNALING", "HALLMARK_GLYCOLYSIS",
                  "HALLMARK_CHOLESTEROL_HOMEOSTASIS", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY")
T9_LINEAGES <- c("Monocyte", "CD4 T", "CD8 T", "B", "NK")

ok("T9_two_contrasts", uniqueN(t9$contrast) == 2L, uniqueN(t9$contrast))
ok("T9_contrasts_named",
   setequal(unique(t9$contrast), c("Sepsis vs Healthy", "Sepsis vs ICU non-sepsis")),
   paste(unique(t9$contrast), collapse = ","))

for (ct in unique(t9$contrast)) {
  s <- t9[contrast == ct]
  tag <- gsub("[^A-Za-z]", "", ct)
  ok(paste0("T9_", tag, "_25_rows"), nrow(s) == 25L, nrow(s))
  ok(paste0("T9_", tag, "_hallmarks_exact"), setequal(unique(s$metric), T9_HALLMARKS),
     paste(setdiff(unique(s$metric), T9_HALLMARKS), collapse = ","))
  ok(paste0("T9_", tag, "_lineages_exact"), setequal(unique(s$cell_type), T9_LINEAGES),
     paste(setdiff(unique(s$cell_type), T9_LINEAGES), collapse = ","))
  ok(paste0("T9_", tag, "_combos_unique"),
     uniqueN(paste(s$metric, s$cell_type)) == 25L, uniqueN(paste(s$metric, s$cell_type)))
  # the 5 x 5 grid must be complete, not merely distinct
  ok(paste0("T9_", tag, "_grid_complete"),
     setequal(paste(s$metric, s$cell_type),
              as.vector(outer(T9_HALLMARKS, T9_LINEAGES, paste))), "grid mismatch")

  # BH over the 25 prespecified comparisons within this comparator
  bh <- p.adjust(s$P, "BH")
  ok(paste0("T9_", tag, "_BH_matches_FDR"), near(bh, s$FDR, tol = 1e-9), max(abs(bh - s$FDR)))

  # support rule: same direction, FDR < 0.05, and interval not crossing zero
  supp <- s$direction_concordant & s$FDR < 0.05 & (s$ci_lower > 0 | s$ci_upper < 0)
  ok(paste0("T9_", tag, "_support_rule_reproduces"),
     all(supp == s$same_direction_supported), sum(supp != s$same_direction_supported))
}

# formal support counts quoted in the manuscript
ok("T9_healthy_supported_14",
   sum(t9[contrast == "Sepsis vs Healthy"]$same_direction_supported) == 14L,
   sum(t9[contrast == "Sepsis vs Healthy"]$same_direction_supported))
ok("T9_icu_supported_0",
   sum(t9[contrast == "Sepsis vs ICU non-sepsis"]$same_direction_supported) == 0L,
   sum(t9[contrast == "Sepsis vs ICU non-sepsis"]$same_direction_supported))

# ---- 5b. T9 summary table must agree with the per-comparison table -------------
ts <- C("T9 Replication summary.csv")
ok("T9_summary_two_rows", nrow(ts) == 2L, nrow(ts))
for (i in seq_len(nrow(ts))) {
  s <- t9[contrast == ts$contrast[i]]
  tag <- gsub("[^A-Za-z]", "", ts$contrast[i])
  ok(paste0("T9_summary_", tag, "_targets_25"), ts$predefined_targets[i] == 25L,
     ts$predefined_targets[i])
  ok(paste0("T9_summary_", tag, "_same_direction"),
     ts$same_direction[i] == sum(s$direction_concordant),
     paste(ts$same_direction[i], "vs", sum(s$direction_concordant)))
  ok(paste0("T9_summary_", tag, "_supported"),
     ts$same_direction_supported[i] == sum(s$same_direction_supported),
     paste(ts$same_direction_supported[i], "vs", sum(s$same_direction_supported)))
  ok(paste0("T9_summary_", tag, "_supported_lineages"),
     ts$supported_lineages[i] == uniqueN(s$cell_type[s$same_direction_supported]),
     paste(ts$supported_lineages[i], "vs", uniqueN(s$cell_type[s$same_direction_supported])))
}

cat("\nFROZEN TARGET INTEGRITY:", if (fails == 0L) "ALL PASS" else paste0(fails, " FAILURES"), "\n")
quit(status = if (fails == 0L) 0L else 1L)
