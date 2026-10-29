# ---------------------------------------------------------------------------
# Table S7: study-clustering sensitivity analysis for k_deg.
#
# The REM enters each curated value as an independent study. For degradation,
# 6 of the 8 standard-error-reporting values come from a single report
# (Lee et al., 2021), so the pooled interval is narrower and I^2 higher than a
# study-level analysis would give. This mirrors the k_act analysis already in
# Figure S4, where the 10 M3R values from one report were collapsed into one
# entry via a within-cluster random-effects fit.
#
# Model specification, filters and the REML estimator are unchanged; only the
# clustering is new. Writes to its own results folder, so no committed figure
# or table from the published pipeline is overwritten.
# ---------------------------------------------------------------------------

if (!is.null(dev.list())) dev.off()
cat("\014")
rm(list = ls())

# Work from this script's folder (code/R), as the sibling scripts do.
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg)) {
  setwd(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
} else if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable() &&
           nzchar(rstudioapi::getSourceEditorContext()$path)) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}
stopifnot("run from code/R (see README)" = dir.exists("visualize") && dir.exists("etc"))

library(metafor)
library(readxl)
library(stringr)

results_path <- "../../results/kdeg_cluster_sensitivity"
dir.create(results_path, recursive = TRUE, showWarnings = FALSE)

# I^2 the way Tables S4 and S6 report it: the Q-based formula, not metafor's
# own $I2 (which uses tau^2 and a typical sampling variance and reads ~0.7
# points higher for this dataset).
i2_q <- function(fit) 100 * (fit$QE - (fit$k - 1)) / fit$QE

# Helpers copied verbatim from rem_OXTR_parameters_log_round2.R -------------
parse_data <- function(data) {
  split_data <- strsplit(data, " ± ")
  mean <- as.numeric(sapply(split_data, function(x) trimws(x[1])))
  sd <- as.numeric(sapply(split_data, function(x) trimws(x[2])))
  data.frame(Mean = mean, SD = sd)
}

create_log_data <- function(data) {
  epsilon <- 1e-10
  data$log_Mean <- log(data$Mean + epsilon)
  data$log_SD <- data$SD / data$Mean
  data
}

# Data: same file, same filter as the published k_deg fit -------------------
filename <- "../../data/OXT_All_Parameters_10162024.xlsx"
OXTR_data <- as.data.frame(read_excel(filename))

kdeg <- OXTR_data[str_detect(OXTR_data$"association rate",
                             pattern = "± (0\\.[0-9]*[1-9][0-9]*|[1-9][0-9]*)") &
                  str_detect(OXTR_data$Symbol, pattern = "Kdeg"), ]
kdeg <- kdeg[!is.na(kdeg$"association rate"), ]
kdeg <- cbind(kdeg, parse_data(kdeg$"association rate"))
kdeg <- create_log_data(kdeg)

cat("\n=== curated k_deg values entering the REM ===\n")
print(data.frame(Receptor = kdeg$Receptor, Reference = kdeg$Reference,
                 Mean = kdeg$Mean, SD = kdeg$SD))

# 1. Control: reproduce the published value-level fit -----------------------
rm_value <- rma(yi = log_Mean, sei = log_SD, data = kdeg)
published <- 4.04e-3
control <- exp(as.numeric(rm_value$b))
cat(sprintf(paste0("\n[control] value-level REM : %.3g min^-1 [%.3g, %.3g]  n=%d  ",
                   "I2(Q-based)=%.1f%%  I2(metafor)=%.1f%%  Q=%.4g  df=%d  tau2=%.4g\n"),
            control, exp(rm_value$ci.lb), exp(rm_value$ci.ub),
            rm_value$k, i2_q(rm_value), rm_value$I2,
            rm_value$QE, rm_value$k - 1, rm_value$tau2))
cat(sprintf("[control] published Table S1: %.3g min^-1  -> ratio %.4f\n",
            published, control / published))

# 2. Study-level fit: collapse each report to a single entry ---------------
# Each cluster is first pooled by its own random-effects fit on the log scale
# (a single-value cluster passes through unchanged), exactly as the k_act M3R
# collapse in Figure S4.
kdeg$cluster <- kdeg$Reference
clusters <- unique(kdeg$cluster)
cat(sprintf("\n%d values from %d distinct reports\n", nrow(kdeg), length(clusters)))
print(table(kdeg$cluster))

collapse <- function(df) {
  if (nrow(df) == 1L) {
    return(data.frame(cluster = df$cluster[1], Receptor = df$Receptor[1],
                      k = 1L, log_Mean = df$log_Mean, log_SD = df$log_SD))
  }
  fit <- rma(yi = log_Mean, sei = log_SD, data = df)
  data.frame(cluster = df$cluster[1],
             Receptor = paste(unique(df$Receptor), collapse = "/"),
             k = nrow(df), log_Mean = as.numeric(fit$b), log_SD = fit$se)
}

study_df <- do.call(rbind, lapply(clusters, function(cl) collapse(kdeg[kdeg$cluster == cl, ])))
cat("\n=== study-level entries ===\n")
print(data.frame(cluster = study_df$cluster, Receptor = study_df$Receptor,
                 values_collapsed = study_df$k,
                 Mean = exp(study_df$log_Mean), log_SE = study_df$log_SD))

rm_study <- rma(yi = log_Mean, sei = log_SD, data = study_df)
study <- exp(as.numeric(rm_study$b))
cat(sprintf(paste0("\n[study-level] REM : %.3g min^-1 [%.3g, %.3g]  n=%d  ",
                   "I2(Q-based)=%.1f%%  I2(metafor)=%.1f%%  Q=%.4g  df=%d  tau2=%.4g  Qp=%.3g\n"),
            study, exp(rm_study$ci.lb), exp(rm_study$ci.ub),
            rm_study$k, i2_q(rm_study), rm_study$I2,
            rm_study$QE, rm_study$k - 1, rm_study$tau2, rm_study$QEp))

# 3. Cluster-robust variance estimate on the value-level fit ---------------
# Keeps all 8 values but corrects the standard error for the clustering.
rm_cr <- rma.mv(yi = log_Mean, V = log_SD^2, random = ~ 1 | cluster/obs,
                data = transform(kdeg, obs = seq_len(nrow(kdeg))), method = "REML")
robust_fit <- robust(rm_cr, cluster = kdeg$cluster, clubSandwich = FALSE)
robust_est <- exp(as.numeric(robust_fit$b))
cat(sprintf("\n[cluster-robust] REM : %.3g min^-1 [%.3g, %.3g]  clusters=%d\n",
            robust_est, exp(robust_fit$ci.lb), exp(robust_fit$ci.ub),
            length(unique(kdeg$cluster))))

# Save ---------------------------------------------------------------------
out <- data.frame(
  Analysis = c("Value-level (published)", "Study-level (clusters collapsed)",
               "Cluster-robust (RVE)"),
  n = c(rm_value$k, rm_study$k, rm_value$k),
  Clusters = c(length(clusters), length(clusters), length(clusters)),
  Estimate_min1 = c(control, study, robust_est),
  CI_low = c(exp(rm_value$ci.lb), exp(rm_study$ci.lb), exp(robust_fit$ci.lb)),
  CI_high = c(exp(rm_value$ci.ub), exp(rm_study$ci.ub), exp(robust_fit$ci.ub)),
  I2_Qbased_percent = c(i2_q(rm_value), i2_q(rm_study), NA),
  I2_metafor_percent = c(rm_value$I2, rm_study$I2, NA),
  Q = c(rm_value$QE, rm_study$QE, NA),
  df = c(rm_value$k - 1, rm_study$k - 1, NA),
  tau2 = c(rm_value$tau2, rm_study$tau2, NA),
  Q_p = c(rm_value$QEp, rm_study$QEp, NA),
  Ratio_to_published = c(control, study, robust_est) / published
)
write.csv(out, file.path(results_path, "kdeg_cluster_sensitivity.csv"), row.names = FALSE)

sessionInfo_txt <- file.path(results_path, "session_info.txt")
writeLines(c(R.version.string,
             paste("metafor", as.character(packageVersion("metafor"))),
             paste("run", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))),
           sessionInfo_txt)

cat("\n=== summary ===\n")
print(out, row.names = FALSE)
cat(sprintf("\nwrote %s\n", file.path(results_path, "kdeg_cluster_sensitivity.csv")))
