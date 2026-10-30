"""Power and precision analysis for the OXTR / Class A GPCR pooled-kinetics
meta-analysis (REM vs BHM).

Answers four questions, all traceable to the same root issue -- most
pooled parameters rest on very few studies (k as low as 3):

  A. Is the non-significant Cochran's Q for k_inact and k_expt a genuine
     power failure (too few studies to detect real heterogeneity), or is
     it evidence that those two parameters really are more homogeneous?
  B. How many studies would be needed to pin the pooled estimate down to a
     stated precision, per parameter?
  C. Is the BHM-vs-REM divergence for k_expt "small" only because n=3 has
     no power to see it, or is the divergence for k_inact real (not a
     small-n artifact)?
  D. Does REM's narrow k_inact estimate come from three studies that
     happen to cluster, or is that just ordinary sampling noise for that
     number of studies?

Data source
-----------
data/OXT_All_Parameters_10162024.xlsx -- the file actually loaded by
code/R/rem_OXTR_parameters_log.R (the REM fitting code Table S1/S4 were
built from). Parsing logic (Symbol/column filters, the "genuine SE
reported" regex, the beta-arrestin x60 unit conversion, the natural-log +
delta-method SE transform) is reproduced line-for-line from that R script
so the SE-reporting subsets here are identical in composition to the ones
metafor's rma() was fit on.

*** Log base: confirmed NATURAL LOG (R's log(), Python's np.log), not
log10 -- see rem_OXTR_parameters_log.R line ~198 (`log(Mean)`) and the BHM
notebook (`np.log`). tau2 and all effect sizes below are on the natural-log
scale unless explicitly labelled "log10". ***

This script is a line-for-line Python re-implementation of the REM
data-loading/fitting logic (REML, as metafor::rma), so the power analysis
runs without an R environment. Its Q, tau2, I2 and p for every parameter are
printed next to the published Table S4 as a check; tau2 is the REML
between-study variance on the natural-log scale (rm$tau2), which is also
what the REM forest plots print (code/R/visualize/forest_ylee.R).

Outputs:
  - results/tables/Table_S6_power_analysis.xlsx
  - results/figures/python/FigureS17_power_curves.png
"""

import platform
import re
import sys
from itertools import combinations
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import scipy
from matplotlib import pyplot as plt
from scipy import optimize, stats

matplotlib.use("Agg")

# ---------------------------------------------------------------------------
# Reproducibility header
# ---------------------------------------------------------------------------
SEED = 20260826  # fixed seed, = date of this analysis (2026-08-26)
RNG = np.random.default_rng(SEED)

print("=" * 78)
print("power_analysis.py -- OXTR/Class A GPCR pooled-kinetics power analysis")
print("=" * 78)
print(f"random seed        : {SEED}")
print(f"python              : {sys.version.split()[0]} ({platform.platform()})")
print(f"numpy               : {np.__version__}")
print(f"pandas              : {pd.__version__}")
print(f"scipy               : {scipy.__version__}")
print(f"matplotlib          : {matplotlib.__version__}")
print("=" * 78)

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_PATH = REPO_ROOT / "data" / "OXT_All_Parameters_10162024.xlsx"
BHM_SUMMARY_PATH = REPO_ROOT / "data" / "GPCR_parameters.xlsx"
TABLE_OUT = REPO_ROOT / "results" / "tables" / "Table_S6_power_analysis.xlsx"
FIG_OUT = REPO_ROOT / "results" / "figures" / "python" / "FigureS17_power_curves.png"

N_SIM = 10_000
ALPHA = 0.05
K_SWEEP = np.arange(3, 31)
K_SWEEP_EXT = np.arange(3, 501)  # extended sweep to find 80%-power k when it exceeds 30

LN10 = np.log(10.0)
Z975 = stats.norm.ppf(0.975)
Z80 = stats.norm.ppf(0.80)

# ---------------------------------------------------------------------------
# 0. Load raw data and replicate the REM parsing pipeline
#    (mirrors code/R/rem_OXTR_parameters_log.R lines 36-222 exactly)
# ---------------------------------------------------------------------------
df_raw = pd.read_excel(DATA_PATH)

SD_PATTERN = re.compile(r"±\s*(0\.[0-9]*[1-9][0-9]*|[1-9][0-9]*)")


def se_reported(series: pd.Series) -> pd.Series:
    """True where a genuine (non-zero) SD was reported after '+/-'."""
    return series.astype(str).apply(lambda s: bool(SD_PATTERN.search(s)))


def parse_mean_sd(series: pd.Series):
    means, sds = [], []
    for s in series:
        left, right = (p.strip() for p in str(s).split("±"))
        means.append(float(left))
        sds.append(float(right))
    return np.array(means), np.array(sds)


def build_se_subset(symbol_pattern, value_col, extra_mask=None, scale=1.0):
    mask = df_raw["Symbol"].astype(str).str.contains(symbol_pattern, regex=True, na=False)
    sub = df_raw[mask].copy()
    sub = sub[sub[value_col].notna()]
    n_full = len(sub)
    if extra_mask is not None:
        sub = sub[extra_mask(sub)]
        n_full = len(sub)
    sub_se = sub[se_reported(sub[value_col])]
    mean, sd = parse_mean_sd(sub_se[value_col])
    return mean * scale, sd * scale, n_full


PARAM_LABELS = {
    "k_act": "G protein activation (kact)",
    "k_inact": "G protein inactivation (kinact)",
    "k_P": "GRK-mediated phosphorylation (kP)",
    "k_dp": "Dephosphorylation (kdp)",
    "k_ba": "β-arrestin binding (kba)",
    "k_ba1": "β-arrestin 1 binding (kba1)",
    "k_ba2": "β-arrestin 2 binding (kba2)",
    "k_cons_int": "Constitutive internalization (kcons,int)",
    "k_int": "Internalization (kint)",
    "k_rec": "Recycling (krec)",
    "k_expt": "Export (kexpt)",
    "k_rgt": "Retrograde transport (krgt)",
    "k_deg": "Degradation (kdeg)",
}

# (name, Symbol regex, value column, extra row filter, unit-scale factor)
SPECS = [
    ("k_act", "Kact, Kinact", "association rate", None, 1.0),
    ("k_inact", "Kact, Kinact", "dissociation rate", None, 1.0),
    ("k_P", "kp, kdp", "association rate", None, 1.0),
    ("k_dp", "kp, kdp", "dissociation rate", None, 1.0),
    ("k_ba", "Kba", "association rate", None, 60.0),
    ("k_ba1", "Kba", "association rate", lambda s: s["Binding Partner"] == "Beta arrestin 1", 60.0),
    ("k_ba2", "Kba", "association rate", lambda s: s["Binding Partner"] == "Beta arrestin 2", 60.0),
    ("k_cons_int", "Kcons_int", "association rate", None, 1.0),
    ("k_int", "Kint", "association rate", None, 1.0),
    ("k_rec", "Krec", "association rate", None, 1.0),
    ("k_expt", "Kins", "association rate", None, 1.0),
    ("k_rgt", "Krgt", "association rate", None, 1.0),
    ("k_deg", "Kdeg", "association rate", None, 1.0),
]

# Table S4 as published (n_se) with the full curated n of Table S1 (used ONLY
# for the verification printout below -- not as an input to any power/precision
# number).
TABLE_S4 = {
    "k_act": dict(tau2=3.11, Q=22900, I2=100, n_se=16, n_full=41),
    "k_inact": dict(tau2=1.02, Q=1230, I2=100, n_se=5, n_full=19),
    "k_P": dict(tau2=0.769, Q=310, I2=99, n_se=3, n_full=8),
    "k_dp": dict(tau2=0.157, Q=13.7, I2=85, n_se=3, n_full=5),
    "k_ba": dict(tau2=0.579, Q=1500, I2=99, n_se=13, n_full=16),
    "k_ba1": dict(tau2=0.493, Q=344, I2=99, n_se=5, n_full=7),
    "k_ba2": dict(tau2=0.614, Q=1150, I2=99, n_se=8, n_full=9),
    "k_cons_int": dict(tau2=0.364, Q=28.5, I2=89, n_se=4, n_full=6),
    "k_int": dict(tau2=0.626, Q=208, I2=97, n_se=7, n_full=12),
    "k_rec": dict(tau2=2.63, Q=502, I2=99, n_se=4, n_full=9),
    "k_expt": dict(tau2=0.0982, Q=3.92, I2=49, n_se=3, n_full=4),
    "k_rgt": dict(tau2=2.32, Q=1090, I2=100, n_se=3, n_full=7),
    "k_deg": dict(tau2=0.927, Q=223, I2=97, n_se=8, n_full=13),
}

# BHM pooled median and 95% HDI, ORIGINAL (back-transformed)
# scale. Source: data/GPCR_parameters.xlsx sheet 'GPCR_parameters_beta12_geomean'
# for all parameters except aggregate k_ba, which is not in that sheet and is
# read from sheet 'GPCR_parameters_geomean' instead. Both sheets hold the
# 2026-09-24 BHM refit (written by build_bhm_summary.py); only their use in the
# MDD formula (part C) is computed here.
bhm_sheet = pd.read_excel(BHM_SUMMARY_PATH, sheet_name="GPCR_parameters_beta12_geomean")
bhm_sheet = bhm_sheet.set_index("Parameter")

BHM_KEY_MAP = {
    "k_act": "kact", "k_inact": "kinact", "k_P": "kP", "k_dp": "kdp",
    "k_ba1": "kba 1", "k_ba2": "kba 2", "k_cons_int": "kcons_int",
    "k_int": "kint", "k_rec": "krec", "k_expt": "kexpt", "k_rgt": "krgt",
    "k_deg": "kdeg",
}
BHM_SUMMARY = {}
for pname, sheet_key in BHM_KEY_MAP.items():
    row = bhm_sheet.loc[sheet_key]
    BHM_SUMMARY[pname] = dict(
        mean=float(row["BHM (mean)"]),
        lo=float(row["BHM (lower [95% CI])"]),
        hi=float(row["BHM (upper [95% CI])"]),
    )
# aggregate k_ba -- not in the beta12 sheet
kba_row = pd.read_excel(BHM_SUMMARY_PATH, sheet_name="GPCR_parameters_geomean").set_index("Parameter").loc["kba"]
BHM_SUMMARY["k_ba"] = dict(mean=float(kba_row["BHM (mean)"]), lo=float(kba_row["BHM (lower [95% CI])"]),
                           hi=float(kba_row["BHM (upper [95% CI])"]))

# ---------------------------------------------------------------------------
# 1. REML random-effects fit (own implementation; matches metafor::rma()
#    default method="REML" on yi = log_Mean, sei = log_SD)
# ---------------------------------------------------------------------------


def reml_tau2(yi: np.ndarray, vi: np.ndarray) -> float:
    def neg_restricted_ll(tau2):
        w = 1.0 / (vi + tau2)
        ybar = np.sum(w * yi) / np.sum(w)
        ll = (
            -0.5 * np.sum(np.log(vi + tau2))
            - 0.5 * np.log(np.sum(w))
            - 0.5 * np.sum(w * (yi - ybar) ** 2)
        )
        return -ll

    res = optimize.minimize_scalar(
        neg_restricted_ll, bounds=(0, 1e4), method="bounded", options={"xatol": 1e-10}
    )
    return max(res.x, 0.0)


def fit_rem(yi: np.ndarray, vi: np.ndarray) -> dict:
    k = len(yi)
    tau2 = reml_tau2(yi, vi)
    w_fe = 1.0 / vi
    ybar_fe = np.sum(w_fe * yi) / np.sum(w_fe)
    Q = float(np.sum(w_fe * (yi - ybar_fe) ** 2))
    p = float(1 - stats.chi2.cdf(Q, df=k - 1))
    I2 = float(max(0.0, (Q - (k - 1)) / Q) * 100) if Q > 0 else 0.0
    w_re = 1.0 / (vi + tau2)
    mu = float(np.sum(w_re * yi) / np.sum(w_re))
    se = float(np.sqrt(1.0 / np.sum(w_re)))
    return dict(k=k, tau2=tau2, Q=Q, p=p, I2=I2, mu=mu, se=se, yi=yi, vi=vi)


PARAMS = {}
for name, sym, col, extra, scale in SPECS:
    mean, sd, n_full = build_se_subset(sym, col, extra, scale)
    eps = 1e-10
    yi = np.log(mean + eps)
    log_sd = sd / mean  # delta method, matches R create_log_data()
    vi = log_sd ** 2
    fit = fit_rem(yi, vi)
    fit["n_full"] = n_full
    fit["mean_bar"] = float(np.exp(fit["mu"]))
    PARAMS[name] = fit

print("\n--- Verification: recomputed vs Table S4 (natural-log scale) ---")
print(f"{'param':12s}{'k':>4s}{'Q(mine)':>11s}{'Q(tbl)':>9s}{'I2(mine)':>10s}"
      f"{'I2(tbl)':>9s}{'tau2(mine)':>12s}{'tau2(tbl)':>11s}{'ratio':>8s}")
for name in PARAM_LABELS:
    f = PARAMS[name]
    t = TABLE_S4[name]
    ratio = t["tau2"] / f["tau2"] if f["tau2"] > 0 else float("nan")
    print(f"{name:12s}{f['k']:4d}{f['Q']:11.2f}{t['Q']:9.1f}{f['I2']:10.2f}"
          f"{t['I2']:9.1f}{f['tau2']:12.4f}{t['tau2']:11.3f}{ratio:8.2f}x")

# ---------------------------------------------------------------------------
# Analysis A: power of Cochran's Q as a function of k
# ---------------------------------------------------------------------------


def simulate_q_power(vi_pool: np.ndarray, tau2: float, mu: float, k: int,
                      n_sim: int, rng: np.random.Generator) -> float:
    """Monte Carlo power of the Q test at sample size k.

    Per replicate: resample k within-study variances (with replacement)
    from the observed vi pool, draw k true effects ~ N(mu, tau2), draw k
    observed effects ~ N(true_i, vi_i), compute Q with inverse-variance
    (1/vi) weights, reject if Q > chi2_(1-alpha, k-1).
    """
    vi_draw = rng.choice(vi_pool, size=(n_sim, k), replace=True)
    true_i = rng.normal(loc=mu, scale=np.sqrt(tau2), size=(n_sim, k))
    obs_i = true_i + rng.normal(loc=0.0, scale=np.sqrt(vi_draw))
    w = 1.0 / vi_draw
    ybar = np.sum(w * obs_i, axis=1) / np.sum(w, axis=1)
    Q = np.sum(w * (obs_i - ybar[:, None]) ** 2, axis=1)
    crit = stats.chi2.ppf(1 - ALPHA, df=k - 1)
    return float(np.mean(Q > crit))


def closed_form_q_power(vi_pool: np.ndarray, tau2: float, k: int) -> float:
    """Noncentral-chi2 approximation (Hardy & Thompson 1998; Hedges &
    Pigott 2001; Jackson 2006). Uses the harmonic-mean-style 'typical'
    within-study variance from the observed vi pool, held fixed across the
    hypothetical k studies (equal-precision approximation)."""
    v_typ = 1.0 / np.mean(1.0 / vi_pool)  # harmonic mean of observed vi
    ncp = tau2 * (k - 1) / v_typ
    crit = stats.chi2.ppf(1 - ALPHA, df=k - 1)
    return float(1 - stats.ncx2.cdf(crit, df=k - 1, nc=ncp))


print("\n--- Analysis A: Cochran's Q power (Monte Carlo, N_SIM="
      f"{N_SIM}, sweep k=3..30) ---")
q_power_results = {}
for name in PARAM_LABELS:
    f = PARAMS[name]
    vi_pool, tau2, mu, n_actual = f["vi"], f["tau2"], f["mu"], f["k"]
    power_curve = np.array([
        simulate_q_power(vi_pool, tau2, mu, k, N_SIM, RNG) for k in K_SWEEP
    ])
    closed_curve = np.array([closed_form_q_power(vi_pool, tau2, k) for k in K_SWEEP])
    power_at_n = simulate_q_power(vi_pool, tau2, mu, n_actual, N_SIM, RNG) if n_actual >= 3 else np.nan

    # smallest k reaching 80% power (extend sweep if needed)
    k80 = None
    for k in K_SWEEP_EXT:
        p = simulate_q_power(vi_pool, tau2, mu, k, 2000, RNG)  # coarser N for the search
        if p >= 0.80:
            k80 = int(k)
            break
    if k80 is not None:
        k80_power = simulate_q_power(vi_pool, tau2, mu, k80, N_SIM, RNG)
    else:
        k80_power = np.nan

    q_power_results[name] = dict(
        power_curve=power_curve, closed_curve=closed_curve,
        power_at_n=power_at_n, k80=k80, k80_power=k80_power,
    )
    print(f"{name:12s} n={n_actual:2d}  power(Q, n)={power_at_n:6.3f}  "
          f"k for 80% power={k80}")

# ---------------------------------------------------------------------------
# Analysis B: studies needed for a target CI half-width, natural log and
# log10 targets (0.5 and 0.25 log10 units)
# ---------------------------------------------------------------------------

HW_TARGETS_LOG10 = {"0.5": 0.5, "0.25": 0.25}

precision_results = {}
for name in PARAM_LABELS:
    f = PARAMS[name]
    tau2 = f["tau2"]
    sigma_bar2 = float(np.mean(f["vi"]))
    row = {"sigma_bar2": sigma_bar2}
    for label, hw10 in HW_TARGETS_LOG10.items():
        hw_nat = hw10 * LN10
        k_needed = (Z975 ** 2) * (tau2 + sigma_bar2) / (hw_nat ** 2)
        # floor at k=2 (minimum size for a random-effects pooled estimate);
        # a value of 2 here means the *current* SE-subset already meets or
        # exceeds this precision target.
        row[f"k_hw_{label}"] = max(2, int(np.ceil(k_needed)))
    precision_results[name] = row

# ---------------------------------------------------------------------------
# Analysis C: minimum detectable difference (MDD) between BHM and REM,
# as a function of k (REM side only -- BHM SE held at its current,
# empirically observed value).
# ---------------------------------------------------------------------------

mdd_results = {}
for name in PARAM_LABELS:
    f = PARAMS[name]
    tau2, sigma_bar2, n_actual, mu_rem = f["tau2"], float(np.mean(f["vi"])), f["k"], f["mu"]

    bhm = BHM_SUMMARY[name]
    se_bhm_log = (np.log(bhm["hi"]) - np.log(bhm["lo"])) / (2 * Z975)
    mu_bhm_log = np.log(bhm["mean"])

    se_rem_curve = np.sqrt((tau2 + sigma_bar2) / K_SWEEP)
    mdd_log_curve = (Z975 + Z80) * np.sqrt(se_bhm_log ** 2 + se_rem_curve ** 2)
    mdd_fold_curve = np.exp(mdd_log_curve)

    se_rem_actual = f["se"]  # actual REML SE at the true current n
    mdd_log_actual = (Z975 + Z80) * np.sqrt(se_bhm_log ** 2 + se_rem_actual ** 2)
    mdd_fold_actual = float(np.exp(mdd_log_actual))

    observed_fold = float(np.exp(abs(mu_rem - mu_bhm_log)))

    # smallest k (>=2) at which the observed fold-difference clears the MDD
    k_detectable = None
    for k in K_SWEEP_EXT:
        se_rem_k = np.sqrt((tau2 + sigma_bar2) / k)
        mdd_k = np.exp((Z975 + Z80) * np.sqrt(se_bhm_log ** 2 + se_rem_k ** 2))
        if mdd_k <= observed_fold:
            k_detectable = int(k)
            break

    mdd_results[name] = dict(
        se_bhm_log=se_bhm_log, se_rem_actual=se_rem_actual,
        mdd_fold_actual=mdd_fold_actual, observed_fold=observed_fold,
        above_mdd=observed_fold > mdd_fold_actual,
        k_detectable=k_detectable, mdd_fold_curve=mdd_fold_curve,
    )
    print(f"\n{name}: observed fold-diff={observed_fold:.3f}x, "
          f"MDD at n={n_actual}: {mdd_fold_actual:.3f}x, "
          f"detectable at n={n_actual}: {observed_fold > mdd_fold_actual}, "
          f"smallest k where 28x-style gap clears MDD: {k_detectable}")

# ---------------------------------------------------------------------------
# Analysis D: bootstrap / exact-enumeration sanity check for k_inact
# ---------------------------------------------------------------------------

# Full BHM-eligible k_inact dataset: Symbol == "Kact, Kinact", any row with
# a numeric dissociation rate (SE reported or not) -- mirrors
# code/python/bhm_refit.py's k_inact selection
# exactly (Symbol == "Kact, Kinact" & dissociation rate notna()).
mask_kinact_full = df_raw["Symbol"].astype(str) == "Kact, Kinact"
kinact_full_rows = df_raw[mask_kinact_full]
kinact_full_rows = kinact_full_rows[kinact_full_rows["dissociation rate"].notna()]
kinact_full_vals = (
    kinact_full_rows["dissociation rate"].astype(str).str.split().str[0].astype(float).values
)
n_full_kinact = len(kinact_full_vals)
log_kinact_full = np.log(kinact_full_vals)

# REM point estimate from the n_rem SE-reporting studies (inverse-variance
# pooled mean, natural-log scale) -- this is PARAMS["k_inact"]["mu"].
rem_kinact_log = PARAMS["k_inact"]["mu"]
rem_kinact_geomean = float(np.exp(rem_kinact_log))
n_rem = PARAMS["k_inact"]["k"]

# Exact enumeration of all C(n_full_kinact, n_rem) subsets of the REM's size
# (n_full_kinact=19, n_rem=5 -> C(19,5)=11628, small enough to enumerate exactly --
# strictly better than Monte Carlo bootstrap for a finite population of
# this size, so we do that instead of resampling with replacement).
subset_log_means = np.array([
    np.mean(log_kinact_full[list(idx)]) for idx in combinations(range(n_full_kinact), n_rem)
])
subset_geomeans = np.exp(subset_log_means)
n_subsets = len(subset_geomeans)

percentile_rank = float(np.mean(subset_log_means <= rem_kinact_log) * 100)
two_sided_tail = float(2 * min(percentile_rank, 100 - percentile_rank) / 100)

print(f"\n--- Analysis D: k_inact {n_rem}-study subsample check ---")
print(f"full BHM-eligible k_inact dataset n = {n_full_kinact} "
      f"(Table S4 states n(full) = {TABLE_S4['k_inact']['n_full']})")
print(f"exact enumeration of all C({n_full_kinact},{n_rem}) = {n_subsets} {n_rem}-study subsets")
print(f"REM ({n_rem} SE-reporting studies) geometric mean = {rem_kinact_geomean:.3f} min^-1")
print(f"REM point estimate percentile rank among all {n_rem}-study subset means: "
      f"{percentile_rank:.1f}th percentile")
print(f"two-sided tail probability = {two_sided_tail:.3f}")

# ---------------------------------------------------------------------------
# Table S6
# ---------------------------------------------------------------------------


def sci(x, sig=3):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "—"
    return f"{x:.{sig-1}e}".replace("e+0", "×10^").replace("e-0", "×10^-").replace("e+", "×10^").replace("e-", "×10^-")


rows = []
for name in PARAM_LABELS:
    f = PARAMS[name]
    qr = q_power_results[name]
    pr = precision_results[name]
    p_str = "p < 0.001" if f["p"] < 0.001 else ("p < 0.01" if f["p"] < 0.01 else
             ("p < 0.05" if f["p"] < 0.05 else "p > 0.05"))
    rows.append({
        "Parameter": PARAM_LABELS[name],
        "n (SE-subset)": f["k"],
        "tau2 (recomputed, ln scale)": f["tau2"],
        "I2 (%)": f["I2"],
        "Q-test p": p_str,
        "Power of Q at current n": qr["power_at_n"],
        "k for 80% power (Q test)": qr["k80"],
        "k for CI half-width < 0.5 log10": pr["k_hw_0.5"],
        "k for CI half-width < 0.25 log10": pr["k_hw_0.25"],
    })

table_s6 = pd.DataFrame(rows)


def fmt_sci(x, sig=3):
    """3-sig-fig scientific notation, matching Table S4's 'X.XX x 10^n' style."""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "—"
    mant, exp = f"{x:.{sig-1}e}".split("e")
    exp = int(exp)
    return f"{mant} x 10^{exp}"


fmt_rows = []
for r in rows:
    fmt_rows.append({
        "Parameter": r["Parameter"],
        "n (SE-subset)": r["n (SE-subset)"],
        "tau2 (recomputed, ln scale)": fmt_sci(r["tau2 (recomputed, ln scale)"]),
        "I2": f"{r['I2 (%)']:.0f}%",
        "Q-test p": r["Q-test p"],
        "Power of Q at current n": f"{r['Power of Q at current n']*100:.0f}%",
        "k for 80% power (Q test)": r["k for 80% power (Q test)"],
        "k for CI half-width < 0.5 log10": r["k for CI half-width < 0.5 log10"],
        "k for CI half-width < 0.25 log10": r["k for CI half-width < 0.25 log10"],
    })
table_s6_fmt = pd.DataFrame(fmt_rows)

TABLE_OUT.parent.mkdir(parents=True, exist_ok=True)
with pd.ExcelWriter(TABLE_OUT, engine="openpyxl") as writer:
    table_s6.to_excel(writer, sheet_name="Table_S6_numeric", index=False)
    table_s6_fmt.to_excel(writer, sheet_name="Table_S6_formatted", index=False)
print(f"\nSaved {TABLE_OUT}")

# ---------------------------------------------------------------------------
# Figure S17: Q-test power curves
# ---------------------------------------------------------------------------

HIGHLIGHT = ["k_inact", "k_expt", "k_act", "k_rec"]  # the two discussed in the text, plus two high-het for contrast
COLORS = {"k_inact": "#1b7837", "k_expt": "#5aae61", "k_act": "#762a83", "k_rec": "#9970ab"}

fig, ax = plt.subplots(figsize=(7.5, 5.5))
for name in PARAM_LABELS:
    curve = q_power_results[name]["power_curve"]
    if name in HIGHLIGHT:
        ax.plot(K_SWEEP, curve, lw=2.5, label=f"{name.replace('_', '')} (τ²={PARAMS[name]['tau2']:.3f})",
                color=COLORS[name], zorder=5)
    else:
        ax.plot(K_SWEEP, curve, lw=1.0, color="0.75", zorder=1)
ax.axhline(0.80, color="black", lw=1, ls="--", alpha=0.6)
ax.text(K_SWEEP[-1], 0.815, "80% power", ha="right", va="bottom", fontsize=9)
ax.set_xlabel("Number of studies (k)")
ax.set_ylabel("Power of Cochran's Q test (α = 0.05)")
ax.set_title("Power to detect observed heterogeneity, by parameter")
ax.set_xlim(3, 30)
ax.set_ylim(0, 1.02)
ax.legend(loc="lower right", fontsize=9, frameon=False)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(FIG_OUT, dpi=300)
plt.close(fig)
print(f"Saved {FIG_OUT}")

print("\nDone.")
