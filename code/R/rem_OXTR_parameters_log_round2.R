# ---------------------------------------------------------------------------
# Plot revisions of rem_OXTR_parameters_log.R used in Figures S6, S11, S13, S16:
#   - phosphorylation forest plot labelled by receptor
#   - k_rec, k_rgt, k_deg and k_cons,int forest plots on a x10^2 axis so the
#     annotations do not collapse to 0.00 [0.00, 0.00]
# Model specification, filters and seed are unchanged; only the x10^2 datasets
# and the plot labels are new. Writes to a separate results folder.
# ---------------------------------------------------------------------------

# Clear all ---------------------------------------------------------------
# Clear plots
if (!is.null(dev.list())) dev.off()
# Clear console
cat("\014")
# Clear workspace
rm(list = ls())

# Change working directory ------------------------------------------------
# Work from this script's folder (code/R) so the relative paths below resolve on any machine:
# Rscript passes --file=, RStudio exposes the open file; otherwise start R in code/R.
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg)) {
  setwd(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
} else if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable() &&
           nzchar(rstudioapi::getSourceEditorContext()$path)) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}
stopifnot("run from code/R (see README)" = dir.exists("visualize") && dir.exists("etc"))

# Add path ----------------------------------------------------------------
subfolders = c("etc", "visualize")
for (i in 1:length(subfolders)) {
  a = list.files(path = subfolders[i], pattern = "[.]R$", full.names = TRUE)
  for (j in 1:length(a)) {
    source(a[j])
  }
}

# Check and generate a result folder --------------------------------------
results_path = "../../results/figures/R/geometric_mean_round2"

# Generate new results folder
dir.create(results_path, recursive = TRUE)

# Load libraries ----------------------------------------------------------
pkg_list = c("ggplot2", "metafor", "readxl", "weights", "latex2exp", "ggpubr", 
             "shades", "ggnewscale", "scales", "ggsignif", "colormap","stringi", 
             "lme4", "haven", "tidyverse","RColorBrewer", "lmerTest","stringr")
instant_pkgs(pkg_list)

# Load data ---------------------------------------------------------------
filename = '../../data/OXT_All_Parameters_10162024.xlsx'
OXTR_data <- as.data.frame(read_excel(filename))

# Clean the data and separate ---------------------------------------------

# KD data
kd_data <- OXTR_data[!str_detect(OXTR_data$Kd, pattern = "± 0.00") 
                     & str_detect(OXTR_data$Symbol, pattern = "Kon, Koff"), ]
kd_data <- kd_data[!is.na(kd_data$Kd), ]

# kong koffg data

koffg <- OXTR_data[str_detect(OXTR_data$"dissociation rate", pattern = "± (0\\.[0-9]*[1-9][0-9]*|[1-9][0-9]*)") 
                   & str_detect(OXTR_data$Symbol, pattern = "Kong, Koffg"), ]
koffg <- koffg[!is.na(koffg$"dissociation rate"), ]

# kp kdp data

kp <- OXTR_data[str_detect(OXTR_data$"association rate", pattern = "± (0\\.[0-9]*[1-9][0-9]*|[1-9][0-9]*)") 
                & str_detect(OXTR_data$Symbol, pattern = "kp, kdp"), ]
kp <- kp[!is.na(kp$"association rate"), ]

kdp <- OXTR_data[str_detect(OXTR_data$"dissociation rate", pattern = "± (0\\.[0-9]*[1-9][0-9]*|[1-9][0-9]*)")  
                 & str_detect(OXTR_data$Symbol, pattern = "kp, kdp"), ]
kdp <- kdp[!is.na(kdp$"dissociation rate"), ]

# kact kinact data
kact <- OXTR_data[str_detect(OXTR_data$"association rate", pattern = "± (0\\.[0-9]*[1-9][0-9]*|[1-9][0-9]*)") 
                  & str_detect(OXTR_data$Symbol, pattern = "Kact, Kinact"), ]
kact <- kact[!is.na(kact$"association rate"), ]

kinact <- OXTR_data[str_detect(OXTR_data$"dissociation rate", pattern = "± (0\\.[0-9]*[1-9][0-9]*|[1-9][0-9]*)")
                    & str_detect(OXTR_data$Symbol, pattern = "Kact, Kinact"), ]
kinact <- kinact[!is.na(kinact$"dissociation rate"), ]

# kba data
#kba <- OXTR_data[OXTR_data$Receptor == 'Oxytocin Receptor' & 
#                 str_detect(OXTR_data$"association rate", pattern = "± (0\\.[0-9]*[1-9][0-9]*|[1-9][0-9]*)") & 
#                 str_detect(OXTR_data$Symbol, pattern = "Kba"), ]
#kba <- kba[!is.na(kba$"association rate"), ]

kba <- OXTR_data[str_detect(OXTR_data$"association rate", pattern = "± (0\\.[0-9]*[1-9][0-9]*|[1-9][0-9]*)") 
                 & str_detect(OXTR_data$Symbol, pattern = "Kba"), ]
kba <- kba[!is.na(kba$"association rate"), ]

# For exact matching of beta arrestin 1
kba_1 <- OXTR_data[str_detect(OXTR_data$"association rate", pattern = "± (0\\.[0-9]*[1-9][0-9]*|[1-9][0-9]*)") 
                   & str_detect(OXTR_data$Symbol, pattern = "Kba")
                   & OXTR_data$"Binding Partner" == "Beta arrestin 1", ]
kba_1 <- kba_1[!is.na(kba_1$"association rate"), ]

# For exact matching of beta arrestin 2
kba_2 <- OXTR_data[str_detect(OXTR_data$"association rate", pattern = "± (0\\.[0-9]*[1-9][0-9]*|[1-9][0-9]*)") 
                   & str_detect(OXTR_data$Symbol, pattern = "Kba")
                   & OXTR_data$"Binding Partner" == "Beta arrestin 2", ]
kba_2 <- kba_2[!is.na(kba_2$"association rate"), ]

# Kint data

kint <- OXTR_data[str_detect(OXTR_data$"association rate", pattern = "± (0\\.[0-9]*[1-9][0-9]*|[1-9][0-9]*)") 
                  & str_detect(OXTR_data$Symbol, pattern = "Kint"), ]
kint <- kint[!is.na(kint$"association rate"), ]

# kcons_int

kcons_int <- OXTR_data[str_detect(OXTR_data$"association rate", pattern = "± (0\\.[0-9]*[1-9][0-9]*|[1-9][0-9]*)")
                       & str_detect(OXTR_data$Symbol, pattern = "Kcons_int"), ]
kcons_int <- kcons_int[!is.na(kcons_int$"association rate"), ]

# krec
krec <- OXTR_data[str_detect(OXTR_data$"association rate", pattern = "± (0\\.[0-9]*[1-9][0-9]*|[1-9][0-9]*)") 
                  & str_detect(OXTR_data$Symbol, pattern = "Krec"), ]
krec <- krec[!is.na(krec$"association rate"), ]

# Kins
Kins <- OXTR_data[str_detect(OXTR_data$"association rate", pattern = "± (0\\.[0-9]*[1-9][0-9]*|[1-9][0-9]*)") 
                  & str_detect(OXTR_data$Symbol, pattern = "Kins"), ]
Kins <- Kins[!is.na(Kins$"association rate"), ]

# Krgt
Krgt <- OXTR_data[str_detect(OXTR_data$"association rate", pattern = "± (0\\.[0-9]*[1-9][0-9]*|[1-9][0-9]*)") 
                  & str_detect(OXTR_data$Symbol, pattern = "Krgt"), ]
Krgt <- Krgt[!is.na(Krgt$"association rate"), ]

# kdeg

kdeg <- OXTR_data[str_detect(OXTR_data$"association rate", pattern = "± (0\\.[0-9]*[1-9][0-9]*|[1-9][0-9]*)")
                  & str_detect(OXTR_data$Symbol, pattern = "Kdeg"), ]
kdeg <- kdeg[!is.na(kdeg$"association rate"), ]

# Extract mean SE and SD data

parse_data <- function(data) {
  # Split the string at "±"
  split_data <- strsplit(data, " ± ")
  
  # Extract the mean and SD values
  mean <- as.numeric(sapply(split_data, function(x) trimws(x[1])))
  sd <- sapply(split_data, function(x) trimws(x[2]))
  
  # Remove "nM" from the SD values if it exists
  #sd <- sub(" nM", "", sd)
  
  # Convert the SD values to numeric
  sd <- as.numeric(sd)
  
  # Create a dataframe with the mean and SD values
  df <- data.frame(Mean = mean, SD = sd)
  
  return(df)
}

# Use the function to parse the data and find the mean and SD
kd_data <- cbind(kd_data, parse_data(kd_data$Kd))

koffg <- cbind(koffg, parse_data(koffg$"dissociation rate"))

kp <- cbind(kp, parse_data(kp$"association rate"))

kdp <- cbind(kdp, parse_data(kdp$"dissociation rate"))

kact <- cbind(kact, parse_data(kact$"association rate"))

kact_0  <-  kact
# Filter rows where Receptor is 'M3R'
M3R_df <- kact_0[kact_0$Receptor == 'M3R',]
# Calculate the mean of the 'Mean' and 'SD' columns
mean_mean <- mean(M3R_df$Mean, na.rm = TRUE)
mean_sd <- mean(M3R_df$SD, na.rm = TRUE)
# Extract a single row where Receptor is 'M3R' and Binding_Partner is 'Gq'
single_row <- kact_0[kact_0$Receptor == 'M3R' & kact_0$SD == '18',]
# Change the 'Mean' and 'SD' values of the single row
single_row$Mean <- mean_mean  # replace new_mean_value with your actual value
single_row$SD <- mean_sd  # replace new_sd_value with your actual value
# Remove the old rows with 'M3R' in the data frame
kact_0 <- kact_0[kact_0$Receptor != 'M3R',]
# Add the new row to the data frame
kact_0 <- rbind(kact_0, single_row)


kinact <- cbind(kinact, parse_data(kinact$"dissociation rate"))

kba <- cbind(kba, parse_data(kba$"association rate")*60)
kba_1 <- cbind(kba_1, parse_data(kba_1$"association rate")*60)
kba_2 <- cbind(kba_2, parse_data(kba_2$"association rate")*60)
kint <- cbind(kint, parse_data(kint$"association rate"))

kcons_intx10 <- cbind(kcons_int, parse_data(kcons_int$"association rate")*100)
kcons_int <- cbind(kcons_int, parse_data(kcons_int$"association rate"))

krecx10 <- cbind(krec, parse_data(krec$"association rate")*100)
krec <- cbind(krec, parse_data(krec$"association rate"))

Kins <- cbind(Kins, parse_data(Kins$"association rate"))
Krgtx10 <- cbind(Krgt, parse_data(Krgt$"association rate")*100)
Krgt <- cbind(Krgt, parse_data(Krgt$"association rate"))

kdegx10 <- cbind(kdeg, parse_data(kdeg$"association rate")*100)
kdeg <- cbind(kdeg, parse_data(kdeg$"association rate"))

# LOG TRANSFORMATION: Create log-transformed columns for meta-analysis
create_log_data <- function(data) {
  # Add small constant to avoid log(0)
  epsilon <- 1e-10
  data$log_Mean <- log(data$Mean + epsilon)
  # For log scale, we approximate the SE using delta method: SE_log ≈ SD/Mean
  data$log_SD <- data$SD / data$Mean
  return(data)
}

# Apply log transformation to all datasets
kd_data <- create_log_data(kd_data)
koffg <- create_log_data(koffg)
kp <- create_log_data(kp)
kdp <- create_log_data(kdp)
kact <- create_log_data(kact)
kact_0 <- create_log_data(kact_0)
# Collapse M3R on the log scale (Figure S4): the arithmetic mean of the 10 M3R values
# (1.44-354 min^-1) with their mean SD gave a CI that ignored the ~250-fold spread, so the
# single M3R row is replaced by a within-M3R random-effects estimate and its SE.
rm_M3R <- rma(yi = log_Mean, sei = log_SD, data = kact[kact$Receptor == 'M3R', ])
kact_0$log_Mean[kact_0$Receptor == 'M3R'] <- as.numeric(rm_M3R$b)
kact_0$log_SD[kact_0$Receptor == 'M3R'] <- rm_M3R$se
kact_0$Mean[kact_0$Receptor == 'M3R'] <- exp(as.numeric(rm_M3R$b))
kinact <- create_log_data(kinact)
kba <- create_log_data(kba)
kba_1 <- create_log_data(kba_1)
kba_2 <- create_log_data(kba_2)
kint <- create_log_data(kint)
kcons_int <- create_log_data(kcons_int)
kcons_intx10 <- create_log_data(kcons_intx10)
krec <- create_log_data(krec)
krecx10 <- create_log_data(krecx10)
Kins <- create_log_data(Kins)
Krgt <- create_log_data(Krgt)
Krgtx10 <- create_log_data(Krgtx10)
kdeg <- create_log_data(kdeg)
kdegx10 <- create_log_data(kdegx10)

# Meta-analysis -----------------------------------------------------------
# Compute weighted average and SD -----------------------------------------

# Kd data - using log-transformed data
rm_Kd <- rma(yi = log_Mean, sei = log_SD, data=kd_data)
summary(rm_Kd)
# Koffg
rm_koffg <- rma(yi = log_Mean, sei = log_SD, data=koffg)
summary(rm_koffg)
# kp
rm_kp <- rma(yi = log_Mean, sei = log_SD, data=kp)
summary(rm_kp)
# kdp
rm_kdp <- rma(yi = log_Mean, sei = log_SD, data=kdp)
# kact
rm_kact <- rma(yi = log_Mean, sei = log_SD, data=kact)
summary(rm_kact)
rm_kact_0 <- rma(yi = log_Mean, sei = log_SD, data=kact_0)
summary(rm_kact_0)

# kinact
rm_kinact <- rma(yi = log_Mean, sei = log_SD, data=kinact)
summary(rm_kdp)
# kba
rm_kba <- rma(yi = log_Mean, sei = log_SD, data=kba)
summary(rm_kba)
rm_kba_1 <- rma(yi = log_Mean, sei = log_SD, data=kba_1)
summary(rm_kba_1)
rm_kba_2 <- rma(yi = log_Mean, sei = log_SD, data=kba_2)
summary(rm_kba_2)

# kint
rm_kint <- rma(yi = log_Mean, sei = log_SD, data=kint)
summary(rm_kint)

# kcont_int
rm_kcons_int <- rma(yi = log_Mean, sei = log_SD, data=kcons_int)
summary(rm_kcons_int)
rm_kcons_intx10 <- rma(yi = log_Mean, sei = log_SD, data=kcons_intx10)
summary(rm_kcons_intx10)

# krec
rm_krec <- rma(yi = log_Mean, sei = log_SD, data=krec)
rm_krecx10 <- rma(yi = log_Mean, sei = log_SD, data=krecx10)
summary(rm_krec)

# Kins
rm_Kins <- rma(yi = log_Mean, sei = log_SD, data=Kins)
summary(rm_Kins)

# Krgt
rm_Krgt <- rma(yi = log_Mean, sei = log_SD, data=Krgt)
rm_Krgtx10 <- rma(yi = log_Mean, sei = log_SD, data=Krgtx10)
summary(rm_Krgt)

# kdeg
rm_kdeg <- rma(yi = log_Mean, sei = log_SD, data=kdeg)
summary(rm_kdeg)

rm_kdegx10 <- rma(yi = log_Mean, sei = log_SD, data=kdegx10)
summary(rm_kdegx10)

# TRANSFORM BACK: Create functions to transform meta-analysis results back to original scale
transform_meta_results <- function(rm_object) {
  # Transform the meta-analysis results back to original scale (exponentiate)
  # For log-normal distribution, the exponentiated mean is the geometric mean
  rm_orig <- rm_object
  # Transform the main estimates
  rm_orig$b <- exp(rm_object$b)  # Back-transform the estimate
  rm_orig$ci.lb <- exp(rm_object$ci.lb)  # Back-transform the lower CI
  rm_orig$ci.ub <- exp(rm_object$ci.ub)  # Back-transform the upper CI
  return(rm_orig)
}

# Transform all meta-analysis results back to original scale
rm_Kd_orig <- transform_meta_results(rm_Kd)
rm_koffg_orig <- transform_meta_results(rm_koffg)
rm_kp_orig <- transform_meta_results(rm_kp)
rm_kdp_orig <- transform_meta_results(rm_kdp)
rm_kact_orig <- transform_meta_results(rm_kact)
rm_kact_0_orig <- transform_meta_results(rm_kact_0)
rm_kinact_orig <- transform_meta_results(rm_kinact)
rm_kba_orig <- transform_meta_results(rm_kba)
rm_kba_1_orig <- transform_meta_results(rm_kba_1)
rm_kba_2_orig <- transform_meta_results(rm_kba_2)
rm_kint_orig <- transform_meta_results(rm_kint)
rm_kcons_int_orig <- transform_meta_results(rm_kcons_int)
rm_kcons_intx10_orig <- transform_meta_results(rm_kcons_intx10)
rm_krec_orig <- transform_meta_results(rm_krec)
rm_krecx10_orig <- transform_meta_results(rm_krecx10)
rm_Kins_orig <- transform_meta_results(rm_Kins)
rm_Krgt_orig <- transform_meta_results(rm_Krgt)
rm_Krgtx10_orig <- transform_meta_results(rm_Krgtx10)
rm_kdeg_orig <- transform_meta_results(rm_kdeg)
rm_kdegx10_orig <- transform_meta_results(rm_kdegx10)

# Student's t-test --------------------------------------------------------

kp_vs_kdp = wtd.t.test(x=kp$Mean, y=kdp$Mean,
                            weight=1/(kp$SD^2+rm_kp_orig$tau2), 
                            weighty=1/(kdp$SD^2+rm_kdp_orig$tau2),
                            alternative="greater", samedata=FALSE)

kcons_int_vs_kint = wtd.t.test(x=kint$Mean, y=kcons_int$Mean,
                               weight=1/(kint$SD^2+rm_kint_orig$tau2), 
                               weighty=1/(kcons_int$SD^2+rm_kcons_int_orig$tau2),
                               alternative="greater", samedata=FALSE)

krec_vs_kdeg = wtd.t.test(x=krec$Mean, y=kdeg$Mean,
                            weight=1/(krec$SD^2+rm_krec_orig$tau2), 
                            weighty=1/(kdeg$SD^2+rm_kdeg_orig$tau2),
                            alternative="greater", samedata=FALSE)

kact_vs_kinact = wtd.t.test(x=kact$Mean, y=kinact$Mean,
                            weight=1/(kact$SD^2+rm_kact_orig$tau2), 
                            weighty=1/(kinact$SD^2+rm_kinact_orig$tau2),
                            alternative="greater", samedata=FALSE)
# kba_1_vs_kba_2
kba_1_vs_kba_2 = wtd.t.test(x=kba_1$Mean, y=kba_2$Mean,
                            weight=1/(kba_1$SD^2+rm_kba_1_orig$tau2), 
                            weighty=1/(kba_2$SD^2+rm_kba_2_orig$tau2),
                            alternative="greater", samedata=FALSE)

# Forest plot -------------------------------------------------------------
# kd
png(file=sprintf("%s/forest_binding_affinity_source_geometric_mean.png", results_path), width=1300, height=700)
forest_ylee(data=kd_data, rm=rm_Kd_orig, atransf=exp, slab=kd_data$"Ligand/Receptor Source",
            header=c("Cell lines", ""),
            unit = stri_unescape_unicode("nM"),
            xlab=TeX("Binding affinity, (nM)"), xlim = c(-3.5, 4.5), alim = c(-1.3, 2.7), cex=2)
dev.off()

png(file=sprintf("%s/forest_binding_affinity_methods_geometric_mean.png", results_path), width=1300, height=700)
forest_ylee(data=kd_data, rm=rm_Kd_orig, atransf=exp, slab=kd_data$Reference,
            unit = stri_unescape_unicode("nM"),
            xlab=TeX("Binding affinity, (nM)"), xlim = c(-4, 7.5), alim = c(-2, 4.5), cex=2)
dev.off()

# kong/koffg # more work needed
png(file=sprintf("%s/forest_gprotein_dissociation_geometric_mean.png", results_path), width=1300, height=600)
forest_ylee(data=koffg, rm=rm_koffg_orig, atransf=exp, slab=koffg$"Binding Partner",
            header=c("Binding Partner", ""),
            unit = stri_unescape_unicode(gsub("<U\\+(....)>", "\\\\u\\1", "min<U+207B><U+00B9>")),
            xlab=TeX("G protein dissociation rate constant, $k_{offg}$ ($min^{-1}$)"), xlim = c(-10, 15), alim = c(-4, 8.5), cex=2)
dev.off()

# kba beta arrestin
# kP and kba are apparent first-order rate constants (min^-1): ln(2)/t1/2 from saturating-agonist assays,
# with no GRK or beta-arrestin concentration reported. Labels only; no value changed.
png(file=sprintf("%s/forest_beta_arrestin_geometric_mean.png", results_path), width=1300, height=700)
forest_ylee(data=kba, rm=rm_kba_orig, atransf=exp, slab=kba$Receptor,
            header=c("Receptor type", ""),
            unit = stri_unescape_unicode(gsub("<U\\+(....)>", "\\\\u\\1", "min<U+207B><U+00B9>")),
            xlab=TeX("Beta arrestin Rate constant, $k_{ba}$ ($min^{-1}$)"), xlim = c(-4, 5.5), alim = c(-2, 3.2), cex=2)
dev.off()
png(file=sprintf("%s/forest_beta_arrestin_Reference_geometric_mean.png", results_path), width=1300, height=700)
forest_ylee(data=kba, rm=rm_kba_orig, atransf=exp, slab=kba$Reference,
            #header=c("Receptor type", ""),
            #title = "Oxytocin receptor ",
            unit = stri_unescape_unicode(gsub("<U\\+(....)>", "\\\\u\\1", "min<U+207B><U+00B9>")),
            xlab=TeX("Beta arrestin rate constant, $k_{ba}$ ($min^{-1}$)"), xlim = c(-4, 5.5), alim = c(-2, 3.2), cex=2)
dev.off()
# kba beta arrestin 1 and 2
png(file=sprintf("%s/forest_beta_arrestin_1_geometric_mean.png", results_path), width=1300, height=500)
forest_ylee(data=kba_1, rm=rm_kba_1_orig, atransf=exp, slab=kba_1$Reference,
            #header=c("Receptor type", ""),
            #title = "β-arrestin 1",
            unit = stri_unescape_unicode(gsub("<U\\+(....)>", "\\\\u\\1", "min<U+207B><U+00B9>")),
            xlab=TeX("Beta arrestin 1 rate constant, $k_{ba}$ ($min^{-1}$)"), xlim = c(-6, 9), alim = c(-2, 5), cex=2)
dev.off()
png(file=sprintf("%s/forest_beta_arrestin_2_geometric_mean.png", results_path), width=1300, height=600)
forest_ylee(data=kba_2, rm=rm_kba_2_orig, atransf=exp, slab=kba_2$Reference,
            #header=c("Receptor type", ""),
            #title = "β-arrestin 2",
            unit = stri_unescape_unicode(gsub("<U\\+(....)>", "\\\\u\\1", "min<U+207B><U+00B9>")),
            xlab=TeX("Beta arrestin 2 rate constant, $k_{ba}$ ($min^{-1}$)"), xlim = c(-6, 8), alim = c(-3, 3.2), cex=2)
dev.off()
png(file=sprintf("%s/forest_beta_arrestin_1_receptor_type_geometric_mean.png", results_path), width=1300, height=500)
forest_ylee(data=kba_1, rm=rm_kba_1_orig, atransf=exp, slab=kba_1$Receptor,
            header=c("Receptor type", ""),
            title = "β-arrestin 1",
            unit = stri_unescape_unicode(gsub("<U\\+(....)>", "\\\\u\\1", "min<U+207B><U+00B9>")),
            xlab=TeX("Beta arrestin 1 rate constant, $k_{ba}$ ($min^{-1}$)"), xlim = c(-6, 9), alim = c(-2, 5), cex=2)
dev.off()
png(file=sprintf("%s/forest_beta_arrestin_2_receptor_type_geometric_mean.png", results_path), width=1300, height=600)
forest_ylee(data=kba_2, rm=rm_kba_2_orig, atransf=exp, slab=kba_2$Receptor,
            header=c("Receptor type", ""),
            title = "β-arrestin 2",
            unit = stri_unescape_unicode(gsub("<U\\+(....)>", "\\\\u\\1", "min<U+207B><U+00B9>")),
            xlab=TeX("Beta arrestin 2 rate constant, $k_{ba}$ ($min^{-1}$)"), xlim = c(-5, 8), alim = c(-3, 3.2), cex=2)
dev.off()


# kint
png(file=sprintf("%s/forest_Internalization_geometric_mean.png", results_path), width=1300, height=500)
forest_ylee(data=kint, rm=rm_kint_orig, atransf=exp, slab=kint$Receptor,
            header=c("Receptor type", ""),
            unit = stri_unescape_unicode(gsub("<U\\+(....)>", "\\\\u\\1", "min<U+207B><U+00B9>")),
            xlab=TeX("Internalization Rate constant, $k_{int}$ ($min^{-1}$)"), xlim = c(-6, 1.7), alim = c(-5, 0.4), cex=2)
dev.off()

# kcons_int
png(file=sprintf("%s/forest_Constitutive_Internalization_geometric_mean.png", results_path), width=1300, height=500)
forest_ylee(data=kcons_int, rm=rm_kcons_int_orig, atransf=exp, slab=kcons_int$Receptor,
            header=c("Receptor type", ""),
            unit = stri_unescape_unicode(gsub("<U\\+(....)>", "\\\\u\\1", "min<U+207B><U+00B9>")),
            xlab=TeX("Constitutive Internalization Rate constant, $k_{cons_int}$ ($min^{-1}$)"), xlim = c(-11, 2.5), alim = c(-9, 0.4), cex=2)
dev.off()

# kcons_intx10
png(file=sprintf("%s/forest_Constitutive_Internalizationx10_geometric_mean.png", results_path), width=1300, height=500)
forest_ylee(data=kcons_intx10, rm=rm_kcons_intx10_orig, atransf=exp, slab=kcons_intx10$Receptor,
            header=c("Receptor type", ""),
            unit = stri_unescape_unicode(gsub("<U\\+(....)>", "\\\\u\\1", "<U+00D7>10<U+207B><U+00B2> min<U+207B><U+00B9>")),
            xlab=TeX("Constitutive Internalization Rate constant, $k_{cons,int}$ (×$10^{-2}$ $min^{-1}$)"), xlim = c(-6, 5.6), alim = c(-3, 2.7), cex=2)
dev.off()

# kp
png(file=sprintf("%s/forest_Phosphorylation_geometric_mean.png", results_path), width=1300, height=500)
forest_ylee(data=kp, rm=rm_kp_orig, atransf=exp,
            slab=paste0(kp$Receptor, " (", kp$Reference, ")"),
            header=c("Receptor type (source)", ""),
            unit = stri_unescape_unicode(gsub("<U\\+(....)>", "\\\\u\\1", "min<U+207B><U+00B9>")),
            xlab=TeX("Phosphorylation Rate constant, $k_{p}$ ($min^{-1}$)"), xlim = c(-9, 9), alim = c(-4.5, 5), cex=2)
dev.off()

# kdp
png(file=sprintf("%s/forest_dephosphorylation_geometric_mean.png", results_path), width=1300, height=500)
forest_ylee(data=kdp, rm=rm_kdp_orig, atransf=exp,
            slab=paste0(kdp$Receptor, " (", kdp$Reference, ")"),
            header=c("Receptor type (source)", ""),
            unit = stri_unescape_unicode(gsub("<U\\+(....)>", "\\\\u\\1", "min<U+207B><U+00B9>")),
            xlab=TeX("Dephosphorylation Rate constant, $k_{dp}$ ($min^{-1}$)"), xlim = c(-5, 1.1), alim = c(-4, 0.12), cex=2)
dev.off()

# kact
png(file=sprintf("%s/forest_Gprotein_activation_geometric_mean.png", results_path), width=1300, height=800)
forest_ylee(data=kact, rm=rm_kact_orig, atransf=exp, slab=kact$Receptor,
            header=c("Receptor type", ""),
            unit = stri_unescape_unicode(gsub("<U\\+(....)>", "\\\\u\\1", "min<U+207B><U+00B9>")),
            xlab=TeX("G protein activation Rate constant, $k_{act}$ ($min^{-1}$)"), xlim = c(-5, 12), alim = c(-1, 7), cex=2)
dev.off()

# kact_0
png(file=sprintf("%s/forest_Gprotein_activation_reformated_geometric_mean.png", results_path), width=1300, height=600)
forest_ylee(data=kact_0, rm=rm_kact_0_orig, atransf=exp, slab=kact_0$Receptor,
            header=c("Receptor type", ""),
            unit = stri_unescape_unicode(gsub("<U\\+(....)>", "\\\\u\\1", "min<U+207B><U+00B9>")),
            xlab=TeX("G protein activation Rate constant, $k_{act}$ ($min^{-1}$)"), xlim = c(-5, 12), alim = c(-1, 7), cex=2)
dev.off()

# kinact
png(file=sprintf("%s/forest_Gprotein_inactivation_geometric_mean.png", results_path), width=1300, height=500)
forest_ylee(data=kinact, rm=rm_kinact_orig, atransf=exp, slab=kinact$Receptor,
            header=c("Receptor type", ""),
            unit = stri_unescape_unicode(gsub("<U\\+(....)>", "\\\\u\\1", "min<U+207B><U+00B9>")),
            xlab=TeX("G protein inactivation Rate constant, $k_{inact}$ ($min^{-1}$)"), xlim = c(-4.5, 8.5), alim = c(0, 3.5), cex=2)
dev.off()

# krec
png(file=sprintf("%s/forest_Recycling_geometric_mean.png", results_path), width=1300, height=450)
forest_ylee(data=krec, rm = rm_krec_orig, atransf=exp, slab=krec$Receptor,
            header=c("Receptor type", ""),
            unit = stri_unescape_unicode(gsub("<U\\+(....)>", "\\\\u\\1", "min<U+207B><U+00B9>")),
            xlab=TeX("Recycling Rate constant, $k_{rec}$ ($min^{-1}$)"), xlim = c(-7.8, 0.9), alim = c(-6, -0.4), cex=2)
dev.off()

# krecx10 (x10^2 axis so the annotations do not collapse to 0.00)
png(file=sprintf("%s/forest_Recyclingx10_geometric_mean.png", results_path), width=1300, height=450)
forest_ylee(data=krecx10, rm = rm_krecx10_orig, atransf=exp, slab=krecx10$Receptor,
            header=c("Receptor type", ""),
            unit = stri_unescape_unicode(gsub("<U\\+(....)>", "\\\\u\\1", "<U+00D7>10<U+207B><U+00B2> min<U+207B><U+00B9>")),
            xlab=TeX("Recycling Rate constant, $k_{rec}$ (\u00d7$10^{-2}$ $min^{-1}$)"), xlim = c(-3.4, 7.3), alim = c(-1.4, 4.2), cex=2)
dev.off()

# Kins
png(file=sprintf("%s/forest_Export_geometric_mean.png", results_path), width=1300, height=500)
forest_ylee(data=Kins, rm=rm_Kins_orig, atransf=exp, slab=Kins$Receptor,
            header=c("Receptor type", ""),
            unit = stri_unescape_unicode(gsub("<U\\+(....)>", "\\\\u\\1", "min<U+207B><U+00B9>")),
            xlab=TeX("Export Rate constant, $k_{expt}$ ($min^{-1}$)"), xlim = c(-9, 2), alim = c(-6, 0.1), cex=2)
dev.off()

# Krgt
png(file=sprintf("%s/forest_Retrograde_transport_geometric_mean.png", results_path), width=1300, height=500)
forest_ylee(data=Krgt, rm=rm_Krgt_orig, atransf=exp, slab=Krgt$Receptor,
            header=c("Receptor type", ""),
            unit = stri_unescape_unicode(gsub("<U\\+(....)>", "\\\\u\\1", "min<U+207B><U+00B9>")),
            xlab=TeX("Retrograde transport Rate constant, $k_{rgt}$ ($min^{-1}$)"), xlim = c(-9, 3), alim = c(-7, 1), cex=2)
dev.off()

# Krgtx10 (x10^2 axis so the annotations do not collapse to 0.00)
png(file=sprintf("%s/forest_Retrograde_transportx10_geometric_mean.png", results_path), width=1300, height=500)
forest_ylee(data=Krgtx10, rm=rm_Krgtx10_orig, atransf=exp, slab=Krgtx10$Receptor,
            header=c("Receptor type", ""),
            unit = stri_unescape_unicode(gsub("<U\\+(....)>", "\\\\u\\1", "<U+00D7>10<U+207B><U+00B2> min<U+207B><U+00B9>")),
            xlab=TeX("Retrograde transport Rate constant, $k_{rgt}$ (\u00d7$10^{-2}$ $min^{-1}$)"), xlim = c(-4.6, 9.6), alim = c(-2.4, 5.6), cex=2)
dev.off()


# kdeg
png(file=sprintf("%s/forest_degradation_geometric_mean.png", results_path), width=1300, height=550)
forest_ylee(data=kdeg, rm=rm_kdeg_orig, atransf=exp, slab=kdeg$Receptor,
            header=c("Receptor type", ""),
            unit = stri_unescape_unicode(gsub("<U\\+(....)>", "\\\\u\\1", "min<U+207B><U+00B9>")),
            xlab=TeX("Degradation Rate constant, $k_{deg}$ ($min^{-1}$)"), xlim = c(-15, 4 ), alim = c(-10, 0.5), cex=2)
dev.off()

# kdegx10
png(file=sprintf("%s/forest_degradation10^2_geometric_mean.png", results_path), width=1300, height=550)
forest_ylee(data=kdegx10, rm=rm_kdegx10_orig, atransf=exp, slab=kdegx10$Receptor,
            header=c("Receptor type", ""),
            unit = stri_unescape_unicode(gsub("<U\\+(....)>", "\\\\u\\1", "<U+00D7>10<U+207B><U+00B2> min<U+207B><U+00B9>")),
            xlab=TeX("Degradation Rate constant, $k_{deg}$ (×$10^{-2}$ $min^{-1}$)"), xlim = c(-5.5, 3.0), alim = c(-4, .9), cex=2)
dev.off()

# Scatter plot ------------------------------------------------------------

# krec_vs_kdeg
p = ggplot() +
  geom_point(data = krec, aes(x = "krec", y = Mean, colour = Receptor), size = 7) +
  geom_point(data = krec, aes(x = "krec", y=rm_krec_orig$b), shape = 95, size=20, colour = "darkblue") +
  labs(color="krec") +
  lightness(scale_color_brewer(palette="Blues"), scalefac(0.8)) +
  guides(color = guide_legend(order=1)) +
  new_scale_color() + 
  geom_point(data = kdeg, aes(x = "kdeg", y = Mean, colour = Receptor), size = 7) +
  geom_point(data = kdeg, aes(x = "kdeg", y=rm_kdeg_orig$b), shape = 95, size=20, colour = "darkgreen") +
  labs(color="kdeg") +
  lightness(scale_color_brewer(palette="Greens"),scalefac(0.8)) + 
  guides(color = guide_legend(order=2)) +
  xlab("") + ylab(TeX(" Rate constant, $(min^{-1})$")) +
  scale_y_continuous(trans= 'log10', breaks=trans_breaks('log10', function(x) 10^x,n = 4),
                     labels=trans_format('log10', math_format(10^.x)), limits = c(1e-6, 1e1),
                     sec.axis = sec_axis(trans=~.*1, name=TeX("$(min^{-1})$"),
                                         breaks=trans_breaks('log10', function(x) 10^x, n= 4),
                                         labels=trans_format('log10', math_format(10^.x)))) +
  scale_x_discrete(limits=c("krec", "kdeg")) +  
  theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, face="bold"))

show(p)
ggsave(sprintf("%s/Degradation and recycling rate constant_geometric_mean.png", results_path), width=4500, height=3000, units="px")
dev.off()

# kp_vs_kdp
p = ggplot() +
  geom_point(data = kp, aes(x = "kp", y = Mean, colour = Method), size = 7) +
  geom_point(data = kp, aes(x = "kp", y=rm_kp_orig$b), shape = 95, size=20, colour = "darkblue") +
  labs(color="kp") +
  lightness(scale_color_brewer(palette="Blues"), scalefac(0.8)) +
  guides(color = guide_legend(order=1)) +
  new_scale_color() + 
  geom_point(data = kdp, aes(x = "kdp", y = Mean, colour = Method), size = 7) +
  geom_point(data = kdp, aes(x = "kdp", y=rm_kdp_orig$b), shape = 95, size=20, colour = "darkgreen") +
  labs(color="kdp") +
  lightness(scale_color_brewer(palette="Greens"),scalefac(0.8)) + 
  guides(color = guide_legend(order=2)) +
  xlab("") + ylab(TeX(" Rate constant, $(min^{-1})$")) +
  scale_y_continuous(trans= 'log10', breaks=trans_breaks('log10', function(x) 10^x,n = 4),
                     labels=trans_format('log10', math_format(10^.x)), limits = c(1e-2, 1e1),
                     sec.axis = sec_axis(trans=~.*1, name=TeX("$(min^{-1})$"),
                                         breaks=trans_breaks('log10', function(x) 10^x, n= 4),
                                         labels=trans_format('log10', math_format(10^.x)))) +
  scale_x_discrete(limits=c("kp", "kdp")) +  
  theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, face="bold"))

show(p)
ggsave(sprintf("%s/Phosphorylation and Dephosphorylation rate constant_geometric_mean.png", results_path), width=4500, height=3000, units="px")
dev.off()

# kact_vs_kinact
p = ggplot() +
  geom_point(data = kact, aes(x = "kact", y = Mean, colour = Method), size = 7) +
  geom_point(data = kact, aes(x = "kact", y=rm_kact_orig$b), shape = 95, size=20, colour = "darkblue") +
  labs(color="kact") +
  lightness(scale_color_brewer(palette="Blues"), scalefac(0.8)) +
  guides(color = guide_legend(order=1)) +
  new_scale_color() +
  geom_point(data = kinact, aes(x = "kinact", y = Mean, colour = Method), size = 7) +
  geom_point(data = kinact, aes(x = "kinact", y=rm_kinact_orig$b), shape = 95, size=20, colour = "darkgreen") +
  labs(color="kinact") +
  lightness(scale_color_brewer(palette="Greens"),scalefac(0.8)) +
  guides(color = guide_legend(order=2)) +
  xlab("") + ylab(TeX(" Rate constant, $(min^{-1})$")) +
  xlab("") + ylab(TeX(" Rate constant, $(min^{-1})$")) +
  scale_y_continuous(trans= 'log10', breaks=trans_breaks('log10', function(x) 10^x,n = 4),
                     labels=trans_format('log10', math_format(10^.x)), limits = c(1e-1, 1e4),
                     sec.axis = sec_axis(trans=~.*1, name=TeX("$(min^{-1})$"),
                                         breaks=trans_breaks('log10', function(x) 10^x, n= 4),
                                         labels=trans_format('log10', math_format(10^.x)))) +
  scale_x_discrete(limits=c("kact", "kinact")) +  
  theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, face="bold"))

show(p)
ggsave(sprintf("%s/G protein activation and inactivation rate constant_geometric_mean.png", results_path), width=4500, height=3000, units="px")
dev.off()

# k_int_vs_kcons_int
p = ggplot() +
  geom_point(data = kint, aes(x = "kint", y = Mean, colour = Receptor), size = 7) +
  geom_point(data = kint, aes(x = "kint", y=rm_kint_orig$b), shape = 95, size=20, colour = "darkblue") +
  labs(color="kint") +
  lightness(scale_color_brewer(palette="Blues"), scalefac(0.8)) +
  guides(color = guide_legend(order=1)) +
  new_scale_color() + 
  geom_point(data = kcons_int, aes(x = "kcons_int", y = Mean, colour = Receptor), size = 7) +
  geom_point(data = kcons_int, aes(x = "kcons_int", y=rm_kcons_int_orig$b), shape = 95, size=20, colour = "darkgreen") +
  labs(color="kcons_int") +
  lightness(scale_color_brewer(palette="Greens"),scalefac(0.8)) + 
  guides(color = guide_legend(order=2)) +
  xlab("") + ylab(TeX(" Rate constant, $(min^{-1})$")) +
  scale_y_continuous(trans= 'log10', breaks=trans_breaks('log10', function(x) 10^x,n = 4),
                     labels=trans_format('log10', math_format(10^.x)), limits = c(1e-3, 1e0),
                     sec.axis = sec_axis(trans=~.*1, name=TeX("$(min^{-1})$"),
                                         breaks=trans_breaks('log10', function(x) 10^x, n= 4),
                                         labels=trans_format('log10', math_format(10^.x)))) +
  scale_x_discrete(limits=c("kint", "kcons_int")) +  
  theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, face="bold"))

show(p)
ggsave(sprintf("%s/Internalization and Constitutive Internalization rate constant_geometric_mean.png", results_path), width=4500, height=3000, units="px")
dev.off()

# kba_1_vs_kba_2
p = ggplot() +
  geom_point(data = kba_1, aes(x = "kba_1", y = Mean, colour = Receptor), size = 7) +
  geom_point(data = kba_1, aes(x = "kba_1", y=rm_kba_1_orig$b), shape = 95, size=20, colour = "darkblue") +
  labs(color="kba_1") +
  lightness(scale_color_brewer(palette="Blues"), scalefac(0.8)) +
  guides(color = guide_legend(order=1)) +
  new_scale_color() + 
  geom_point(data = kba_2, aes(x = "kba_2", y = Mean, colour = Receptor), size = 7) +
  geom_point(data = kba_2, aes(x = "kba_2", y=rm_kba_2_orig$b), shape = 95, size=20, colour = "darkgreen") +
  labs(color="kba_2") +
  lightness(scale_color_brewer(palette="Greens"),scalefac(0.8)) + 
  guides(color = guide_legend(order=2)) +
  xlab("") + ylab(TeX(" Rate constant, $(min^{-1})$")) +
  scale_y_continuous(trans= 'log10', breaks=trans_breaks('log10', function(x) 10^x,n = 4),
                     labels=trans_format('log10', math_format(10^.x)), limits = c(1e-2, 1e2),
                     sec.axis = sec_axis(trans=~.*1, name=TeX("$(min^{-1})$"),
                                         breaks=trans_breaks('log10', function(x) 10^x, n= 4),
                                         labels=trans_format('log10', math_format(10^.x)))) +
  scale_x_discrete(limits=c("kba_1", "kba_2")) +
  theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, face="bold"))

show(p)
ggsave(sprintf("%s/Comparison of kba_1 and kba_2_geometric_mean.png", results_path), width=4500, height=3000, units="px")
dev.off()

# Extracelular scatter plot -------------------------------------------------------
p = ggplot() +
  geom_point(data = kact, aes(x = "kact", y = Mean, colour = Receptor), size = 7) +
  geom_point(data = kact, aes(x = "kact", y=rm_kact_orig$b), shape = 95, size=20, colour = "darkblue") +
  labs(color="kact") +
  lightness(scale_color_brewer(palette="Blues"), scalefac(0.8)) +
  guides(color = guide_legend(order=1)) +
  new_scale_color() +
  geom_point(data = kinact, aes(x = "kinact", y = Mean, colour = Receptor), size = 7) +
  geom_point(data = kinact, aes(x = "kinact", y=rm_kinact_orig$b), shape = 95, size=20, colour = "darkgreen") +
  labs(color="kinact") +
  lightness(scale_color_brewer(palette="Greens"),scalefac(0.8)) +
  guides(color = guide_legend(order=2)) +
  new_scale_color() +
  geom_point(data = kp, aes(x = "kp", y = Mean, colour = Receptor), size = 7) +
  geom_point(data = kp, aes(x = "kp", y=rm_kp_orig$b), shape = 95, size=20, colour = "yellow") +
  labs(color="kp") +
  lightness(scale_color_brewer(palette="YlOrBr"),scalefac(0.8)) +
  guides(color = guide_legend(order=3)) +
  new_scale_color() +
  geom_point(data = kdp, aes(x = "kdp", y = Mean, colour = Receptor), size = 7) +
  geom_point(data = kdp, aes(x = "kdp", y=rm_kdp_orig$b), shape = 95, size=20, colour = "darkred") +
  labs(color="kdp") +
  lightness(scale_color_brewer(palette="Reds"),scalefac(0.8)) +
  guides(color = guide_legend(order=4)) +
  new_scale_color() +
  geom_point(data = kba, aes(x = "kba", y = Mean, colour = Receptor), size = 7) +
  geom_point(data = kba, aes(x = "kba", y=rm_kba_orig$b), shape = 95, size=20, colour = "darkgreen") +
  labs(color="kba") +
  lightness(scale_color_brewer(palette="Greens"),scalefac(0.8)) +
  guides(color = guide_legend(order=5)) +
  new_scale_color() +
  xlab("") + ylab(TeX(" Rate constant, $(min^{-1})$")) +
  scale_y_continuous(trans= 'log10', breaks=trans_breaks('log10', function(x) 10^x,n = 4),
                     labels=trans_format('log10', math_format(10^.x)), limits = c(1e-2, 1e4),
                     sec.axis = sec_axis(trans=~.*1, name=TeX("$(min^{-1})$"),
                                         breaks=trans_breaks('log10', function(x) 10^x, n= 4),
                                         labels=trans_format('log10', math_format(10^.x)))) +
  scale_x_discrete(limits=c("kact", "kinact", "kp", "kdp", "kba")) +
  theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, face="bold"))
      
show(p)
ggsave(sprintf("%s/Membrane Dynamics plot_geometric_mean.png", results_path), width=4500, height=3000, units="px")
dev.off()

# Intracelular scatter plot -------------------------------------------------------
p = ggplot() +
  geom_point(data = kint, aes(x = "kint", y = Mean, colour = Receptor), size = 7) +
  geom_point(data = kint, aes(x = "kint", y=rm_kint_orig$b), shape = 95, size=20, colour = "darkblue") +
  labs(color="kint") +
  lightness(scale_color_brewer(palette="Blues"), scalefac(0.8)) +
  guides(color = guide_legend(order=1)) +
  new_scale_color() +
  geom_point(data = kcons_int, aes(x = "kcons_int", y = Mean, colour = Receptor), size = 7) +
  geom_point(data = kcons_int, aes(x = "kcons_int", y=rm_kcons_int_orig$b), shape = 95, size=20, colour = "yellow") +
  labs(color="kcons_int") +
  lightness(scale_color_brewer(palette="Purples"),scalefac(0.8)) +
  guides(color = guide_legend(order=2)) +
  new_scale_color() +
  geom_point(data = krec, aes(x = "krec", y = Mean, colour = Receptor), size = 7) +
  geom_point(data = krec, aes(x = "krec", y=rm_krec_orig$b), shape = 95, size=20, colour = "darkred") +
  labs(color="krec") +
  lightness(scale_color_brewer(palette="Reds"),scalefac(0.8)) +
  guides(color = guide_legend(order=3)) +
  new_scale_color() +
  geom_point(data = kdeg, aes(x = "kdeg", y = Mean, colour = Receptor), size = 7) +
  geom_point(data = kdeg, aes(x = "kdeg", y=rm_kdeg_orig$b), shape = 95, size=20, colour = "darkgreen") +
  labs(color="kdeg") +
  lightness(scale_color_brewer(palette="Greens"),scalefac(0.8)) +
  guides(color = guide_legend(order=3)) +
  new_scale_color() +
  xlab("") + ylab(TeX(" Rate constant, $(min^{-1})$")) +
  scale_y_continuous(trans= 'log10', breaks=trans_breaks('log10', function(x) 10^x,n = 4),
                     labels=trans_format('log10', math_format(10^.x)), limits = c(1e-5, 1e1),
                     sec.axis = sec_axis(trans=~.*1, name=TeX("$(min^{-1})$"),
                                         breaks=trans_breaks('log10', function(x) 10^x, n= 4),
                                         labels=trans_format('log10', math_format(10^.x)))) +
  scale_x_discrete(limits=c("kint","kcons_int", "krec", "kdeg")) +
  theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, face="bold"))

show(p)
ggsave(sprintf("%s/Intracellular Dynamics plot_geometric_mean.png", results_path), width=4500, height=3000, units="px")
dev.off()

# Whole cell scatter plot -------------------------------------------------------
p = ggplot() +
  geom_point(data = kact, aes(x = "kact", y = Mean, colour = Receptor), size = 13.5) +
  geom_point(data = kact, aes(x = "kact", y=rm_kact_orig$b), shape = 95, size=40, colour = "darkblue") +
  labs(color="kact") +
  lightness(scale_color_brewer(palette="Blues"), scalefac(0.8)) +
  guides(color = guide_legend(order=1)) +
  new_scale_color() +
  geom_point(data = kinact, aes(x = "kinact", y = Mean, colour = Receptor), size = 13.5, shape = 1, stroke = 5) +
  geom_point(data = kinact, aes(x = "kinact", y=rm_kinact_orig$b), shape = 95, size=40, colour = "darkgreen") +
  labs(color="kinact") +
  lightness(scale_color_brewer(palette="Greens"),scalefac(0.8)) +
  guides(color = guide_legend(order=2)) +
  new_scale_color() +
  geom_point(data = kp, aes(x = "kp", y = Mean, colour = Receptor), size = 13.5, shape = 15, stroke = 5) +
  geom_point(data = kp, aes(x = "kp", y=rm_kp_orig$b), shape = 95, size=40, colour = "darkred") +
  labs(color="kp") +
  lightness(scale_color_brewer(palette="Reds"),scalefac(0.8)) +
  guides(color = guide_legend(order=3)) +
  new_scale_color() +
  geom_point(data = kdp, aes(x = "kdp", y = Mean, colour = Receptor), size = 13.5, shape = 0, stroke = 5) +
  geom_point(data = kdp, aes(x = "kdp", y=rm_kdp_orig$b), shape = 95, size=40, colour = "darkgreen") +
  labs(color="kdp") +
  lightness(scale_color_brewer(palette="Greens"),scalefac(0.8)) +
  guides(color = guide_legend(order=4)) +
  new_scale_color() +
  geom_point(data = kba_1, aes(x = "kba 1", y = Mean, colour = Receptor), size = 13.5, shape = 18, stroke = 5) +
  geom_point(data = kba_1, aes(x = "kba 1", y=rm_kba_1_orig$b), shape = 95, size=40, colour = "yellow") +
  labs(color="kba 1") +
  lightness(scale_color_brewer(palette="YlOrBr"),scalefac(0.8)) +
  guides(color = guide_legend(order=5)) +
  new_scale_color() +
  geom_point(data = kba_2, aes(x = "kba 2", y = Mean, colour = Receptor), size = 13.5, shape = 9, stroke = 5) +
  geom_point(data = kba_2, aes(x = "kba 2", y=rm_kba_2_orig$b), shape = 95, size=40, colour = "darkred") +
  labs(color="kba 2") +
  lightness(scale_color_brewer(palette="Reds"),scalefac(0.8)) +
  guides(color = guide_legend(order=6)) +
  new_scale_color() +
  geom_point(data = kint, aes(x = "kint", y = Mean, colour = Receptor), size = 13.5, shape = 17) +
  geom_point(data = kint, aes(x = "kint", y=rm_kint_orig$b), shape = 95, size=40, colour = "darkblue") +
  labs(color="kint") +
  lightness(scale_color_brewer(palette="Blues"), scalefac(0.8)) +
  guides(color = guide_legend(order=7)) +
  new_scale_color() +
  geom_point(data = kcons_int, aes(x = "kcons_int", y = Mean, colour = Receptor), size = 13.5, shape = 2, stroke = 5) +
  geom_point(data = kcons_int, aes(x = "kcons_int", y=rm_kcons_int_orig$b), shape = 95, size=40, colour = "darkred") +
  labs(color="kcons_int") +
  lightness(scale_color_brewer(palette="Reds"),scalefac(0.8)) +
  guides(color = guide_legend(order=8)) +
  new_scale_color() +
  geom_point(data = krec, aes(x = "krec", y = Mean, colour = Receptor), size = 13.5, shape = 3, stroke = 5) +
  geom_point(data = krec, aes(x = "krec", y=rm_krec_orig$b), shape = 95, size=40, colour = "darkgreen") +
  labs(color="krec") +
  lightness(scale_color_brewer(palette="Greens"),scalefac(0.8)) +
  guides(color = guide_legend(order=9)) +
  new_scale_color() +
  geom_point(data = Krgt, aes(x = "Krgt", y = Mean, colour = Receptor), size = 13.5, shape = 14, stroke = 5) +
  geom_point(data = Krgt, aes(x = "Krgt", y=rm_Krgt_orig$b), shape = 95, size=40, colour = "darkblue") +
  labs(color="Krgt") +
  lightness(scale_color_brewer(palette="Blues"), scalefac(0.8)) +
  guides(color = guide_legend(order=10)) +
  new_scale_color() +
  geom_point(data = Kins, aes(x = "Kexpt", y = Mean, colour = Receptor), size = 13.5, shape = 8, stroke = 5) +
  geom_point(data = Kins, aes(x = "Kexpt", y=rm_Kins_orig$b), shape = 95, size=40, colour = "darkred") +
  labs(color="Kexpt") +
  lightness(scale_color_brewer(palette="Reds"),scalefac(0.8)) +
  guides(color = guide_legend(order=11)) +
  new_scale_color() +
  geom_point(data = kdeg, aes(x = "kdeg", y = Mean, colour = Receptor), size = 13.5, shape = 4, stroke = 5) +
  geom_point(data = kdeg, aes(x = "kdeg", y=rm_kdeg_orig$b), shape = 95, size=40, colour = "yellow") +
  labs(color="kdeg") +
  lightness(scale_color_brewer(palette="YlOrBr"),scalefac(0.8)) +
  guides(color = guide_legend(order=12)) +
  new_scale_color() +
  xlab("") + ylab(TeX(" Rate constant, $(min^{-1})$")) +
  scale_y_continuous(trans= 'log10', breaks=trans_breaks('log10', function(x) 10^x,n = 4),
                     labels=trans_format('log10', math_format(10^.x)), limits = c(1e-4, 1e3),
                     sec.axis = sec_axis(trans=~.*1, name=TeX("$(min^{-1})$"),
                                         breaks=trans_breaks('log10', function(x) 10^x, n= 4),
                                         labels=trans_format('log10', math_format(10^.x)))) +
  scale_x_discrete(limits=c("kact", "kinact","kp", "kdp", "kba 1", "kba 2", "kint", "kcons_int", "krec", "Krgt", "Kexpt", "kdeg"),
                   labels = c(expression(k[act]),expression(k[inact]),expression(k[p]),expression(k[dp]),expression(k[ba1]),expression(k[ba2]),expression(k[int]),expression(k[cons_int]),expression(k[rec]),expression(K[rgt]),expression(K[expt]), expression(k[deg]))
  ) +
  guides(color = guide_legend(order=1, ncol=1)) +  # Stack legend items vertically +
  theme(text = element_text(size = 55),
        plot.title = element_text(hjust = .5, face="bold"),
        axis.text = element_text(size = 57),  # Increase the size of axis texts
        legend.direction = "horizontal",  # Stack legend items vertically
        legend.position = "right",  # Position the legend on the right side
        legend.spacing = unit(1.5, 'cm')  # Increase space between legend labels
  )

show(p)
ggsave(sprintf("%s/Whole cell Dynamics plot_geometric_mean.png", results_path), width=14999, height=6500, units="px") # for vertical w = 11000,h=11000
dev.off()

## grant

p = ggplot() +
  geom_point(data = kact, aes(x = "kact", y = Mean, colour = Receptor), size = 13.5) +
  geom_point(data = kact, aes(x = "kact", y=rm_kact_orig$b), shape = 95, size=40, colour = "darkblue") +
  labs(color="kact") +
  lightness(scale_color_brewer(palette="Blues"), scalefac(0.8)) +
  guides(color = guide_legend(order=1)) +
  new_scale_color() +
  geom_point(data = kinact, aes(x = "kinact", y = Mean, colour = Receptor), size = 13.5, shape = 1, stroke = 5) +
  geom_point(data = kinact, aes(x = "kinact", y=rm_kinact_orig$b), shape = 95, size=40, colour = "darkgreen") +
  labs(color="kinact") +
  lightness(scale_color_brewer(palette="Greens"),scalefac(0.8)) +
  guides(color = guide_legend(order=2)) +
  new_scale_color() +
  geom_point(data = kba, aes(x = "kba", y = Mean, colour = Receptor), size = 13.5, shape = 18, stroke = 5) +
  geom_point(data = kba, aes(x = "kba", y=rm_kba_orig$b), shape = 95, size=40, colour = "yellow") +
  labs(color="kba") +
  lightness(scale_color_brewer(palette="YlOrBr"),scalefac(0.8)) +
  guides(color = guide_legend(order=5)) +
  new_scale_color() +
  geom_point(data = kint, aes(x = "kint", y = Mean, colour = Receptor), size = 13.5, shape = 17) +
  geom_point(data = kint, aes(x = "kint", y=rm_kint_orig$b), shape = 95, size=40, colour = "darkblue") +
  labs(color="kint") +
  lightness(scale_color_brewer(palette="Blues"), scalefac(0.8)) +
  guides(color = guide_legend(order=6)) +
  new_scale_color() +
  geom_point(data = krec, aes(x = "krec", y = Mean, colour = Receptor), size = 13.5, shape = 3, stroke = 5) +
  geom_point(data = krec, aes(x = "krec", y=rm_krec_orig$b), shape = 95, size=40, colour = "darkgreen") +
  labs(color="krec") +
  lightness(scale_color_brewer(palette="Greens"),scalefac(0.8)) +
  guides(color = guide_legend(order=8)) +
  new_scale_color() +
  geom_point(data = Krgt, aes(x = "krgt", y = Mean, colour = Receptor), size = 13.5, shape = 14, stroke = 5) +
  geom_point(data = Krgt, aes(x = "krgt", y=rm_Krgt_orig$b), shape = 95, size=40, colour = "darkblue") +
  labs(color="krgt") +
  lightness(scale_color_brewer(palette="Blues"), scalefac(0.8)) +
  guides(color = guide_legend(order=9)) +
  new_scale_color() +
  geom_point(data = Kins, aes(x = "kins", y = Mean, colour = Receptor), size = 13.5, shape = 8, stroke = 5) +
  geom_point(data = Kins, aes(x = "kins", y=rm_Kins_orig$b), shape = 95, size=40, colour = "darkred") +
  labs(color="kins") +
  lightness(scale_color_brewer(palette="Reds"),scalefac(0.8)) +
  guides(color = guide_legend(order=10)) +
  new_scale_color() +
  geom_point(data = kdeg, aes(x = "kdeg", y = Mean, colour = Receptor), size = 13.5, shape = 4, stroke = 5) +
  geom_point(data = kdeg, aes(x = "kdeg", y=rm_kdeg_orig$b), shape = 95, size=40, colour = "yellow") +
  labs(color="kdeg") +
  lightness(scale_color_brewer(palette="YlOrBr"),scalefac(0.8)) +
  guides(color = guide_legend(order=11)) +
  new_scale_color() +
  xlab("") + ylab(TeX(" Rate constant, $(min^{-1})$")) +
  scale_y_continuous(trans= 'log10', breaks=trans_breaks('log10', function(x) 10^x,n = 4),
                     labels=trans_format('log10', math_format(10^.x)), limits = c(1e-4, 1e3),
                     sec.axis = sec_axis(trans=~.*1, name=TeX("$(min^{-1})$"),
                                         breaks=trans_breaks('log10', function(x) 10^x, n= 4),
                                         labels=trans_format('log10', math_format(10^.x)))) +
  scale_x_discrete(limits=c("kact", "kinact", "kba", "kint", "krec", "krgt", "kins", "kdeg"),
                   labels = c(expression(k[act]),expression(k[inact]),expression(k[ba]),expression(k[int]),expression(k[rec]),expression(k[rgt]),expression(k[ins]), expression(k[deg]))
  ) +
  guides(color = guide_legend(order=1, ncol=1)) +  # Stack legend items vertically +
  theme(text = element_text(size = 55),
        plot.title = element_text(hjust = .5, face="bold"),
        axis.text = element_text(size = 65),  # Increase the size of axis texts
        axis.text.x = element_text(size = 65, face = "bold"),  # Make x-axis text bold
        legend.direction = "horizontal",  # Stack legend items vertically
        legend.position = "right",  # Position the legend on the right side
        legend.spacing = unit(1.5, 'cm'),  # Increase space between legend labels
        panel.grid.major = element_line(color = "black"),  # Set major grid lines to black
        panel.grid.minor = element_line(color = "black")  # Set minor grid lines to black
  ) +
  theme(axis.text.x = element_text(color = c("black", "black", "black", "red", "red", "black", "red", "red")))  # Set specific colors for x-axis labels

show(p)
ggsave(sprintf("%s/Whole cell Dynamics plot1_geometric_mean.png", results_path), width=14000, height=6000, units="px") # for vertical w = 11000,h=11000
dev.off()

# FIRST-ORDER PARAMETERS PLOT -------------------------------------------------
p_first_order = ggplot() +
  # kact - G protein activation
  geom_point(data = kact, aes(x = "kact", y = Mean, colour = Receptor), size = 13.5) +
  geom_point(data = kact, aes(x = "kact", y=rm_kact_orig$b), shape = 95, size=40, colour = "darkblue") +
  labs(color="kact") +
  lightness(scale_color_brewer(palette="Blues"), scalefac(0.8)) +
  guides(color = guide_legend(order=1)) +
  new_scale_color() +
  
  # kinact - G protein inactivation
  geom_point(data = kinact, aes(x = "kinact", y = Mean, colour = Receptor), size = 13.5, shape = 1, stroke = 5) +
  geom_point(data = kinact, aes(x = "kinact", y=rm_kinact_orig$b), shape = 95, size=40, colour = "darkgreen") +
  labs(color="kinact") +
  lightness(scale_color_brewer(palette="Greens"),scalefac(0.8)) +
  guides(color = guide_legend(order=2)) +
  new_scale_color() +
  
  # kdp - GRK-mediated dephosphorylation
  geom_point(data = kdp, aes(x = "kdp", y = Mean, colour = Receptor), size = 13.5, shape = 0, stroke = 5) +
  geom_point(data = kdp, aes(x = "kdp", y=rm_kdp_orig$b), shape = 95, size=40, colour = "darkgreen") +
  labs(color="kdp") +
  lightness(scale_color_brewer(palette="Greens"),scalefac(0.8)) +
  guides(color = guide_legend(order=3)) +
  new_scale_color() +
  
  # kint - Phosphorylated GPCR internalization
  geom_point(data = kint, aes(x = "kint", y = Mean, colour = Receptor), size = 13.5, shape = 17) +
  geom_point(data = kint, aes(x = "kint", y=rm_kint_orig$b), shape = 95, size=40, colour = "darkblue") +
  labs(color="kint") +
  lightness(scale_color_brewer(palette="Blues"), scalefac(0.8)) +
  guides(color = guide_legend(order=4)) +
  new_scale_color() +
  
  # kcons_int - Constitutive internalization
  geom_point(data = kcons_int, aes(x = "kcons_int", y = Mean, colour = Receptor), size = 13.5, shape = 2, stroke = 5) +
  geom_point(data = kcons_int, aes(x = "kcons_int", y=rm_kcons_int_orig$b), shape = 95, size=40, colour = "darkred") +
  labs(color="kcons_int") +
  lightness(scale_color_brewer(palette="Reds"),scalefac(0.8)) +
  guides(color = guide_legend(order=5)) +
  new_scale_color() +
  
  # krec - GPCR recycling
  geom_point(data = krec, aes(x = "krec", y = Mean, colour = Receptor), size = 13.5, shape = 3, stroke = 5) +
  geom_point(data = krec, aes(x = "krec", y=rm_krec_orig$b), shape = 95, size=40, colour = "darkgreen") +
  labs(color="krec") +
  lightness(scale_color_brewer(palette="Greens"),scalefac(0.8)) +
  guides(color = guide_legend(order=6)) +
  new_scale_color() +
  
  # Krgt - Retrograde transport
  geom_point(data = Krgt, aes(x = "Krgt", y = Mean, colour = Receptor), size = 13.5, shape = 14, stroke = 5) +
  geom_point(data = Krgt, aes(x = "Krgt", y=rm_Krgt_orig$b), shape = 95, size=40, colour = "darkblue") +
  labs(color="Krgt") +
  lightness(scale_color_brewer(palette="Blues"), scalefac(0.8)) +
  guides(color = guide_legend(order=7)) +
  new_scale_color() +
  
  # Kexpt - Golgi/ER/nucleus export
  geom_point(data = Kins, aes(x = "Kexpt", y = Mean, colour = Receptor), size = 13.5, shape = 8, stroke = 5) +
  geom_point(data = Kins, aes(x = "Kexpt", y=rm_Kins_orig$b), shape = 95, size=40, colour = "darkred") +
  labs(color="Kexpt") +
  lightness(scale_color_brewer(palette="Reds"),scalefac(0.8)) +
  guides(color = guide_legend(order=8)) +
  new_scale_color() +
  
  # kdeg - GPCR degradation
  geom_point(data = kdeg, aes(x = "kdeg", y = Mean, colour = Receptor), size = 13.5, shape = 4, stroke = 5) +
  geom_point(data = kdeg, aes(x = "kdeg", y=rm_kdeg_orig$b), shape = 95, size=40, colour = "yellow") +
  labs(color="kdeg") +
  lightness(scale_color_brewer(palette="YlOrBr"),scalefac(0.8)) +
  guides(color = guide_legend(order=9)) +
  new_scale_color() +
  
  # Axis and labels (EXACT SAME as original)
  xlab("") + ylab(TeX(" Rate constant, $(min^{-1})$")) +
  scale_y_continuous(trans= 'log10', breaks=trans_breaks('log10', function(x) 10^x,n = 4),
                     labels=trans_format('log10', math_format(10^.x)), limits = c(1e-4, 1e3),
                     sec.axis = sec_axis(trans=~.*1, name=TeX("$(min^{-1})$"),
                                         breaks=trans_breaks('log10', function(x) 10^x, n= 4),
                                         labels=trans_format('log10', math_format(10^.x)))) +
  scale_x_discrete(limits=c("kact", "kinact", "kdp", "kint", "kcons_int", "krec", "Krgt", "Kexpt", "kdeg"),
                   labels = c(expression(k[act]),expression(k[inact]),expression(k[dp]),expression(k[int]),
                             expression(k[cons_int]),expression(k[rec]),expression(K[rgt]),expression(K[expt]), 
                             expression(k[deg]))) +
  guides(color = guide_legend(order=1, ncol=1)) +
  theme(text = element_text(size = 55),
        plot.title = element_text(hjust = .5, face="bold"),
        axis.text = element_text(size = 57),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        legend.direction = "horizontal",
        legend.position = "right",
        legend.spacing = unit(1.5, 'cm'),
        panel.background = element_rect(fill = "gray90"),
        panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white")
  )

# Show and save
show(p_first_order)
ggsave(sprintf("%s/First_Order_Whole_cell_Dynamics_plot.png", results_path), 
       width=12000, height=6500, units="px", dpi=300)
dev.off()

# APPARENT FIRST-ORDER kP / kba PLOT (labelled second-order in earlier versions) ---
p_second_order = ggplot() +
  # kp - GRK-mediated phosphorylation
  geom_point(data = kp, aes(x = "kp", y = Mean, colour = Receptor), size = 13.5, shape = 15, stroke = 5) +
  geom_point(data = kp, aes(x = "kp", y=rm_kp_orig$b), shape = 95, size=40, colour = "darkred") +
  labs(color="kp") +
  lightness(scale_color_brewer(palette="Reds"),scalefac(0.8)) +
  guides(color = guide_legend(order=1)) +
  new_scale_color() +
  
  # kba_1 - β-arrestin binding
  geom_point(data = kba_1, aes(x = "kba 1", y = Mean, colour = Receptor), size = 13.5, shape = 18, stroke = 5) +
  geom_point(data = kba_1, aes(x = "kba 1", y=rm_kba_1_orig$b), shape = 95, size=40, colour = "yellow") +
  labs(color="kba 1") +
  lightness(scale_color_brewer(palette="YlOrBr"),scalefac(0.8)) +
  guides(color = guide_legend(order=2)) +
  new_scale_color() +
  
  # kba_2 - β-arrestin binding
  geom_point(data = kba_2, aes(x = "kba 2", y = Mean, colour = Receptor), size = 13.5, shape = 9, stroke = 5) +
  geom_point(data = kba_2, aes(x = "kba 2", y=rm_kba_2_orig$b), shape = 95, size=40, colour = "darkred") +
  labs(color="kba 2") +
  lightness(scale_color_brewer(palette="Reds"),scalefac(0.8)) +
  guides(color = guide_legend(order=3)) +
  new_scale_color() +
  
  # Axis and labels (EXACT SAME as original)
  xlab("") + ylab(TeX(" Rate constant, $(min^{-1})$")) +
  scale_y_continuous(trans= 'log10', breaks=trans_breaks('log10', function(x) 10^x,n = 4),
                     labels=trans_format('log10', math_format(10^.x)), limits = c(1e-4, 1e3),
                     sec.axis = sec_axis(trans=~.*1, name=TeX("$(min^{-1})$"),
                                         breaks=trans_breaks('log10', function(x) 10^x, n= 4),
                                         labels=trans_format('log10', math_format(10^.x)))) +
  scale_x_discrete(limits=c("kp", "kba 1", "kba 2"),
                   labels = c(expression(k[p]),expression(k[ba1]),expression(k[ba2]))) +
  guides(color = guide_legend(order=1, ncol=1)) +
  theme(text = element_text(size = 55),
        plot.title = element_text(hjust = .5, face="bold"),
        axis.text = element_text(size = 57),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        legend.direction = "vertical",
        legend.position = "right",
        legend.spacing = unit(1.5, 'cm'),
        panel.background = element_rect(fill = "gray90"),
        panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white")
  )

# Show and save
show(p_second_order)
ggsave(sprintf("%s/Second_Order_Whole_cell_Dynamics_plot.png", results_path), 
       width=8000, height=6500, units="px", dpi=300)
dev.off()
# ---- round-2 verification: pooled values must match Table S1 ----------------
cat("\n==== ROUND2 POOLED CHECK (geometric means, min^-1) ====\n")
chk <- function(label, rm_orig, scale = 1) {
  cat(sprintf("%-28s %.6g [%.6g, %.6g]  n=%d\n", label,
              as.numeric(rm_orig$b)/scale, as.numeric(rm_orig$ci.lb)/scale,
              as.numeric(rm_orig$ci.ub)/scale, rm_orig$k))
}
chk("kp", rm_kp_orig)
chk("krec", rm_krec_orig)
chk("krec (from x10^2)", rm_krecx10_orig, 100)
chk("krgt", rm_Krgt_orig)
chk("krgt (from x10^2)", rm_Krgtx10_orig, 100)
chk("kdeg", rm_kdeg_orig)
chk("kdeg (from x10^2)", rm_kdegx10_orig, 100)
cat("==== END ROUND2 POOLED CHECK ====\n\n")
