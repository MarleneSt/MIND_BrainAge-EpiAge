###### run_all.R: run the full meta-analysis pipeline (scripts 01-21) ------
# Run from the Reproducibilty folder:   Rscript run_all.R
# To resume from a given script (e.g. after a failure):   Rscript run_all.R 02
#
# Expects:  data/     cohort-level result files (same layout as analysis_210526),
#                     the moderator files (moderators_*_R.csv, moderators_*_models.csv)
#                     cohort_age_sd.csv (SD of age per cohort x timepoint, used by R4.3)
#                     and descriptives_table_20260625.csv (cohort descriptives, used by Fig 1)
#           Nothing is written to data/.
# Writes:   results/  all model outputs, tables and figures
#           results/logs/  one log per script, run times and sessionInfo
#
# Each script runs in its own R session (as when run on its own), in order.
# Scripts 20 (Fig 4B heatmap) and 21 (Fig 1 descriptives) are Python and are run with python3.
# The pipeline stops at the first script that fails; see its log.


#### settings ----
scripts <- sort(list.files(".", pattern = "^[0-9]{2}_.*\\.(R|py)$"))

# optional: start from a given script number (earlier scripts are not rerun)
args  <- commandArgs(trailingOnly = TRUE)
start <- if (length(args) > 0) args[1] else "01"
scripts <- scripts[substr(scripts, 1, 2) >= start]

stopifnot(dir.exists("data"))
dir.create("results/logs", recursive = TRUE, showWarnings = FALSE)


#### run scripts in order ----
run_times <- data.frame(script = scripts, minutes = NA_real_, status = NA_character_)

# when resuming, keep the run times of earlier scripts from the previous run
if (start != "01" && file.exists("results/logs/run_times.csv")) {
  previous  <- read.csv("results/logs/run_times.csv")
  previous  <- previous[substr(previous$script, 1, 2) < start & previous$status %in% "ok", ]
} else {
  previous  <- run_times[0, ]
}

for (i in seq_along(scripts)) {
  f   <- scripts[i]
  log <- file.path("results/logs", sub("\\.(R|py)$", ".log", f))
  cat(sprintf("[%s] running %s ...", format(Sys.time(), "%H:%M:%S"), f))

  t0     <- Sys.time()
  cmd    <- if (grepl("\\.py$", f)) "python3" else "Rscript"
  status <- system2(cmd, f, stdout = log, stderr = log)
  run_times$minutes[i] <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  run_times$status[i]  <- ifelse(status == 0, "ok", "FAILED")

  cat(sprintf(" %s (%.1f min)\n", run_times$status[i], run_times$minutes[i]))
  if (status != 0) {
    write.csv(rbind(previous, run_times), "results/logs/run_times.csv", row.names = FALSE)
    stop("Pipeline stopped: ", f, " failed. See ", log)
  }
}

run_times <- rbind(previous, run_times)
write.csv(run_times, "results/logs/run_times.csv", row.names = FALSE)
cat(sprintf("\nAll %d scripts finished in %.1f minutes.\n", nrow(run_times), sum(run_times$minutes)))


#### key checks: do the results match the paper? ----
suppressMessages(library(metafor))

mod2_eff2 <- readRDS("results/mod2_eff2.rds")
global    <- read.csv("results/assoc_model2_PAR_global_final/global_estimates.csv")
omnibus   <- read.csv("results/revision/R2.5/assoc_model2_PAR_modelwise_final/omnibus_differences_main_text_R2.5.csv")
fdr_assoc <- read.csv("results/revision/fdr/moderator_fdr_assoc_revision.csv")

checks <- data.frame(
  check    = c("association effect sizes (k)",
               "global association beta",
               "model pair omnibus QM(135)",
               "model pair pFDR",
               "generation (4 level) pFDR",
               "proportion female pFDR",
               "time interval pFDR"),
  expected = c(2448, 0.017, 192.00, 0.004, 0.044, 0.004, 0.809),
  observed = c(nrow(mod2_eff2),
               round(global$estimate[global$term == "intrcpt"], 3),
               round(omnibus$QM[omnibus$test == "assoc_pairs"], 2),
               round(fdr_assoc$p_assoc_fdr[fdr_assoc$moderator == "model pair"], 3),
               round(fdr_assoc$p_assoc_fdr[fdr_assoc$moderator == "generation (4 level)"], 3),
               round(fdr_assoc$p_assoc_fdr[fdr_assoc$moderator == "proportion female"], 3),
               round(fdr_assoc$p_assoc_fdr[fdr_assoc$moderator == "time interval DNAm-MRI"], 3))
)
checks$match <- ifelse(abs(checks$observed - checks$expected) < 1e-3, "yes", "NO")
print(checks, row.names = FALSE)
write.csv(checks, "results/logs/key_checks.csv", row.names = FALSE)


#### software versions ----
writeLines(capture.output(sessionInfo()), "results/logs/sessionInfo.txt")
