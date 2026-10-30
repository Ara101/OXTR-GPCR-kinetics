# Table S5 - scope sensitivity, REM half ------------------
# Refit the random-effects model for every parameter whose curated values
# include non-Class A data (Class B or C GPCRs, or non-GPCR proteins), once on
# the full set and once on Class A GPCRs only, and compare.
#
# Receptor classes come from data/receptor_class_map.csv. Filters, the log
# transform and the delta-method SE are identical to rem_OXTR_parameters_log.R,
# so the "full" rows reproduce Table S1 exactly (checked below).
#
# Output: results/tables/Table_S5_scope_sensitivity_REM.csv, which
# code/python/scope_sensitivity_bhm.py merges with the BHM refits into
# results/tables/Table_S5_scope_sensitivity.xlsx / .csv.
#
# k_on,G, k_off,G and K_d have no REM estimate in Table S1 (too few values with
# a standard error), so they appear in the BHM half only. k_on,G and K_d are
# all Class A anyway.

if (!is.null(dev.list())) dev.off()
rm(list = ls())

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg)) {
  setwd(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
} else if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable() &&
           nzchar(rstudioapi::getSourceEditorContext()$path)) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}
stopifnot("run from code/R (see README)" = dir.exists("visualize") && dir.exists("etc"))

pkg_list <- c("readxl", "stringr", "metafor")
for (p in pkg_list) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

OXTR_data <- as.data.frame(read_excel("../../data/OXT_All_Parameters_10162024.xlsx"))
class_map <- read.csv("../../data/receptor_class_map.csv", stringsAsFactors = FALSE, encoding = "UTF-8")
receptor_class <- function(r) {
  cl <- class_map$class[match(trimws(r), class_map$receptor)]
  if (any(is.na(cl))) stop("receptor missing from data/receptor_class_map.csv: ",
                           paste(unique(trimws(r)[is.na(cl)]), collapse = ", "))
  cl
}

# Same SE filter as rem_OXTR_parameters_log.R: keep values whose reported +/- is non-zero
se_pattern <- "± (0\\.[0-9]*[1-9][0-9]*|[1-9][0-9]*)"
select_rows <- function(symbol, column) {
  d <- OXTR_data[str_detect(OXTR_data[[column]], pattern = se_pattern) &
                   str_detect(OXTR_data$Symbol, pattern = symbol), ]
  d <- d[!is.na(d[[column]]), ]
  parts <- strsplit(d[[column]], " ± ")
  d$Mean <- as.numeric(sapply(parts, function(x) trimws(x[1])))
  d$SD <- as.numeric(sapply(parts, function(x) trimws(x[2])))
  d$log_Mean <- log(d$Mean + 1e-10)
  d$log_SD <- d$SD / d$Mean
  d$class <- receptor_class(d$Receptor)
  d
}

# parameter -> (workbook Symbol, value column); same as the REM pipeline
params <- list(
  kact   = c("Kact, Kinact", "association rate"),
  kinact = c("Kact, Kinact", "dissociation rate"),
  krec   = c("Krec", "association rate"),
  kexpt  = c("Kins", "association rate"),
  krgt   = c("Krgt", "association rate"),
  kdeg   = c("Kdeg", "association rate")
)

fit_row <- function(name, subset, d) {
  base <- data.frame(parameter = name, subset = subset, n = nrow(d),
                     receptors = paste(d$Receptor, collapse = "; "),
                     mean = NA, ci_lb = NA, ci_ub = NA, tau2 = NA, Q = NA, I2_Q = NA, note = "")
  if (nrow(d) < 2) {
    base$note <- sprintf("not computable: %d SE-reporting value%s remain%s; REM pooling needs n >= 2",
                         nrow(d), ifelse(nrow(d) == 1, "", "s"), ifelse(nrow(d) == 1, "s", ""))
    return(base)
  }
  rm <- rma(yi = log_Mean, sei = log_SD, data = d)
  base$mean <- exp(as.numeric(rm$b)); base$ci_lb <- exp(rm$ci.lb); base$ci_ub <- exp(rm$ci.ub)
  base$tau2 <- rm$tau2; base$Q <- rm$QE
  base$I2_Q <- max(0, 100 * (rm$QE - (rm$k - 1)) / rm$QE)  # Q-based I^2, as in Table S4
  base
}

rows <- list(); classification <- list()
for (name in names(params)) {
  d <- select_rows(params[[name]][1], params[[name]][2])
  classification[[name]] <- data.frame(parameter = name, Reference = d$Reference, Receptor = d$Receptor,
                                       class = d$class, value = d$Mean, SD = d$SD)
  rows[[paste(name, "full")]] <- fit_row(name, "full", d)
  rows[[paste(name, "classA")]] <- fit_row(name, "Class A only", d[d$class == "Class A", ])
}
out <- do.call(rbind, rows)

# Sanity check: the full fits must reproduce Table S1 (REM column)
table_s1 <- c(kact = 3.04e1, kinact = 7.02e0, krec = 2.37e-2, kexpt = 3.75e-2, krgt = 1.31e-2, kdeg = 4.04e-3)
full <- out[out$subset == "full", ]
stopifnot(all(abs(signif(full$mean, 3) - table_s1[full$parameter]) / table_s1[full$parameter] < 1e-6))

print(out[, c("parameter", "subset", "n", "mean", "ci_lb", "ci_ub", "I2_Q", "note")], row.names = FALSE)
dir.create("../../results/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(out, "../../results/tables/Table_S5_scope_sensitivity_REM.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(do.call(rbind, classification), "../../results/tables/Table_S5_scope_sensitivity_REM_rows.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\nSaved results/tables/Table_S5_scope_sensitivity_REM.csv and _REM_rows.csv\n")
