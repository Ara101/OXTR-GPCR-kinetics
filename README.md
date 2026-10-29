# OXTR-GPCR-kinetics

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23782402.svg)](https://doi.org/10.5281/zenodo.23782402)

Data and code for **"Kinetic timescales and measurement gaps across oxytocin receptor and GPCR signaling and trafficking"**

Keith Lionel Tukei<sup>1</sup>, Preeti Dubey<sup>1</sup>, Yingye Fang<sup>1</sup>, Yunjeong Lee<sup>1</sup>, Erin L. Reinl<sup>2</sup>, Antonina I. Frolova<sup>2</sup>, Sarah K. England<sup>2</sup>, Princess I. Imoukhuede<sup>1</sup>

<sup>1</sup> Department of Bioengineering, University of Washington, Seattle, WA, USA
<sup>2</sup> Center for Reproductive Health Sciences, Department of Obstetrics and Gynecology, WashU Medicine, St. Louis, MO, USA

---

## Overview

Synthetic oxytocin is given during roughly half of U.S. births, yet the effective dose can vary up to 45-fold between patients. Differences in how the oxytocin receptor (OXTR) signals and traffics may contribute to that variability, but the underlying rate constants have rarely been measured. This repository holds the curated dataset and the analysis code behind a meta-analysis of kinetic measurements for OXTR and related Class A G protein-coupled receptors (GPCRs), drawn from 57 studies identified through a structured literature search.

Fourteen parameters were pooled: thirteen rate constants spanning G-protein binding and activation, receptor phosphorylation, β-arrestin recruitment, internalization, recycling, retrograde transport, membrane export and degradation, plus the GPCR:G-protein binding affinity (K<sub>d</sub>). Two frameworks were used:

- **Bayesian hierarchical models (BHM)**, fit with PyMC 5.7.2: a receptor-level random intercept with partial pooling toward a global mean (ln k = μ + α<sub>receptor</sub> + ε). The BHM uses every curated value, including those without a reported uncertainty. Reported values are posterior medians with 95% highest-density intervals (NUTS, 4 chains × 2,000 draws after 2,000 tuning steps, target_accept 0.99, seed 8927; all R̂ ≤ 1.01).
- **Random-effects models (REM)**, fit with `metafor` in R (restricted maximum likelihood) on the values that report a standard error. Reported values are geometric means with 95% confidence intervals.

## Pooled estimates

| Process | Symbol | BHM | REM | Units |
|---|---|---|---|---|
| GPCR:G-protein binding affinity | K<sub>d</sub> | 1.90×10³ | — | nM |
| G-protein association | k<sub>on,G</sub> | 1.82×10⁶ | — | M⁻¹min⁻¹ |
| G-protein dissociation | k<sub>off,G</sub> | 1.35×10⁰ | — | min⁻¹ |
| G-protein activation | k<sub>act</sub> | 3.27×10⁰ | 2.48×10¹ | min⁻¹ |
| G-protein inactivation | k<sub>inact</sub> | 1.81×10⁰ | 1.40×10¹ | min⁻¹ |
| GRK-mediated phosphorylation | k<sub>P</sub> | 6.33×10⁻¹ | 4.14×10⁻¹ | min⁻¹ |
| Dephosphorylation | k<sub>dp</sub> | 5.75×10⁻² | 6.31×10⁻² | min⁻¹ |
| β-arrestin binding | k<sub>ba</sub> | 6.75×10⁻¹ | 8.94×10⁻¹ | min⁻¹ |
| Constitutive internalization | k<sub>cons,int</sub> | 1.02×10⁻² | 4.47×10⁻³ | min⁻¹ |
| Agonist-induced internalization | k<sub>int</sub> | 6.48×10⁻² | 5.88×10⁻² | min⁻¹ |
| Recycling | k<sub>rec</sub> | 2.36×10⁻² | 2.37×10⁻² | min⁻¹ |
| Golgi/ER export to membrane | k<sub>expt</sub> | 1.36×10⁻² | 3.75×10⁻² | min⁻¹ |
| Retrograde transport to Golgi | k<sub>rgt</sub> | 4.50×10⁻² | 1.31×10⁻² | min⁻¹ |
| Degradation | k<sub>deg</sub> | 1.88×10⁻³ | 4.04×10⁻³ | min⁻¹ |

K<sub>d</sub>, k<sub>on,G</sub> and k<sub>off,G</sub> are BHM-only: too few values report a standard error to support a random-effects fit. k<sub>P</sub> and k<sub>ba</sub> are apparent first-order rate constants derived as ln(2)/t<sub>½</sub> from saturating-agonist assays. Intervals and per-parameter statistics are in the manuscript's Supplementary Tables S1–S7 and in the files listed below.

## Where each result comes from

| Manuscript item | Script | Output |
|---|---|---|
| Figure 1 (pathway schematic with BHM intervals; created in BioRender, see License) | intervals from `results/bhm_refit_2026-09-24/tables/bhm_refit_summary.csv` | `results/figures/Figure1_signaling_schematic_with_values.png` |
| Figure 2 (curation workflow) | — | `results/figures/Figure2_flow_chart.png` |
| Figure 3 (receptors and ligands) | `code/python/figure3_pie_charts.py` | `results/comparative_analysis/` |
| Figures 4–5 (pooled rate constants) | `code/python/figure4_figure5_rate_constants.py` | `results/figures/python/` |
| BHM estimates (Table S1, Figures S1, S2, S5, S6, S10, S12, S14–S16) | `code/python/bhm_refit.py` | `results/bhm_refit_2026-09-24/` |
| REM estimates (Table S1, Figures S3, S4, S7–S9, S11, S15) | `code/R/rem_OXTR_parameters_log.R` | `results/figures/R/geometric_mean/` |
| REM panels of Figures S6, S11 (B), S13 and S16 (rescaled axes) | `code/R/rem_OXTR_parameters_log_round2.R` | `results/figures/R/geometric_mean_round2/` |
| Pooled summary used by Figures 4–5 and Tables S1–S3 | `code/python/build_bhm_summary.py` | `data/GPCR_parameters.xlsx` |
| Tables S2–S3 (BHM vs REM statistics) | `code/python/build_bhm_summary.py` | `results/bhm_refit_2026-09-24/table_s2_s3_statistics.csv` |
| Table S4 (heterogeneity) | `code/R/rem_OXTR_parameters_log.R` | forest-plot statistics; `data/GPCR_parameters.xlsx` (`rem_heterogeneity`) |
| Table S5 (Class A-only sensitivity) | `code/R/table_s5_scope_sensitivity.R`, then `code/python/scope_sensitivity_bhm.py` | `results/tables/Table_S5_scope_sensitivity.*`, `results/bhm_scope_sensitivity_2026-10-08/` |
| Table S6 and Figure S17 (power and precision) | `code/python/power_analysis.py` | `results/tables/Table_S6_power_analysis.xlsx`, `results/figures/python/FigureS17_power_curves.png` |
| Table S7 (degradation study clustering) | `code/R/kdeg_cluster_sensitivity.R` | `results/kdeg_cluster_sensitivity/` |
| k<sub>off,G</sub>/k<sub>on,G</sub> vs K<sub>d</sub> check (Results) | `code/python/kinetic_kd_check.py` | `results/tables/kinetic_Kd_check.csv` |

`code/python/1_comparative_analysis_bhm_rem.ipynb` is the exploratory BHM vs REM comparison notebook from which the scripts above were extracted; its outputs are cleared.

## Data

`data/OXT_All_Parameters_10162024.xlsx` (sheet `all_params`) is the curated dataset: one row per reported measurement, with the receptor, binding partner, ligand, method, cell type, value and uncertainty where reported, the source study (SR column, Supplementary References SR1–SR57) and its DOI. Receptor labels follow the source papers. Kinetic values were extracted from published figures and tables, using ImageJ and WebPlotDigitizer where values had to be read off plots. Where a study reported a half-life, the rate constant was derived as k = ln(2)/t<sub>½</sub>. All values were log-transformed before modeling and back-transformed afterwards.

| File | Contents |
|---|---|
| `OXT_All_Parameters_10162024.xlsx` | Curated dataset read by the BHM, REM, Table S5, Table S7 and power analyses |
| `OXT_All_Parameters_01282026.xlsx` | The same rows; only the export symbol differs (Kins → Kexpt, 4 rows). Read by Figures 4–5 |
| `GPCR_parameters.xlsx` | Pooled BHM and REM estimates (`*_geomean` sheets), β-arrestin isoform estimates, REM heterogeneity, and the receptor/ligand list behind Figure 3 |
| `receptor_class_map.csv` | Receptor → Class A / B / C / non-GPCR, used by the Table S5 scripts |
| `Hasbi2004*.csv`, `Namkung2018*.xlsx`, `Wright2021*.xlsx` | Values digitised from individual source figures |

Throughout, `n` counts **curated values, not distinct publications**: the random-effects fit treats each value as an independent study, and several reports contribute more than one value (Table S7 examines the effect of this for degradation).

## Reproducing the analysis

### BHM (Python)

```bash
conda env create -f code/python/environment_bhm.yml   # or: micromamba create -f ...
conda activate oxtr_bhm
python code/python/bhm_refit.py                  # every BHM model -> results/bhm_refit_2026-09-24/
python code/python/build_bhm_summary.py          # BHM columns of data/GPCR_parameters.xlsx + Table S2/S3 statistics
python code/python/figure3_pie_charts.py
python code/python/figure4_figure5_rate_constants.py
python code/python/power_analysis.py             # Table S6, Figure S17
python code/python/kinetic_kd_check.py
```

`environment_bhm.yml` pins the exact conda-forge build used for the reported values (PyMC 5.7.2, pytensor 2.14.2, arviz 0.16.1, numpy 1.25.2, Python 3.11, with the m2w64 g++ toolchain so pytensor compiles on Windows). With the same seed and environment the refit reproduces its tables exactly.

### REM (R)

R 4.4.2 with the packages in `requirements-r.txt` (missing packages are installed on first run).

```bash
Rscript code/R/rem_OXTR_parameters_log.R            # REM estimates + forest plots -> results/figures/R/geometric_mean/
Rscript code/R/rem_OXTR_parameters_log_round2.R     # rescaled forest plots -> results/figures/R/geometric_mean_round2/
Rscript code/R/table_s5_scope_sensitivity.R         # Table S5, REM half
python  code/python/scope_sensitivity_bhm.py        # Table S5, BHM half + merged table (Python environment)
Rscript code/R/kdeg_cluster_sensitivity.R           # Table S7
```

The R scripts set their working directory to `code/R`, so they can be run from any folder with `Rscript` or opened in RStudio. `code/python/power_analysis.py` reimplements the REM fit in Python (REML) and regenerates the heterogeneity and power statistics without R.

## Citation

Please cite the associated manuscript (details will be added on publication). To cite the data and code, use the Zenodo archive:

> Tukei KL, Dubey P, Fang Y, Lee Y, Reinl EL, Frolova AI, England SK, Imoukhuede PI. 2026. OXTR-GPCR-kinetics: data and code for a meta-analysis of oxytocin receptor and GPCR signaling and trafficking rate constants (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.23782403

Version 1.0.0 is the release analysed in the manuscript; https://doi.org/10.5281/zenodo.23782402 always resolves to the latest version. `CITATION.cff` holds the citation metadata.

## License

- **Code** (`code/`): MIT License, see `LICENSE`.
- **Data and results** (`data/`, `results/`): Creative Commons Attribution 4.0 International (CC BY 4.0), see `LICENSE-DATA`. The curated values come from the published studies listed in the Supplementary References (SR1–SR57); please credit those sources when reusing individual values.
- **Figure 1** (`results/figures/Figure1_signaling_schematic_with_values.png`): Created in BioRender. Tukei, K. L. (2026) https://BioRender.com/s9ppa7p. The BioRender artwork is used under a BioRender publication licence and is not covered by CC BY 4.0.
