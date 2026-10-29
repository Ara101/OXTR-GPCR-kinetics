"""
Bayesian hierarchical refit of every pooled parameter (2026-09-24).

Replaces the earlier specification (0_bhm_trafficking_freatures.ipynb, kept in the development repository), which entered study, method and
binding partner as integer category codes with informative slopes (so the reported "global mean" was
the prediction for whichever study happened to be coded 0), and whose fits reported R-hat > 1.01 and
hundreds of divergences.

Model, for curated value i of one parameter (natural-log scale):

    log k_i = mu + a_receptor[i] + e_i
    a_r ~ Normal(0, sigma_receptor),  e_i ~ Normal(0, sigma)
    mu ~ Normal(mean(log k), 3);  sigma_receptor, sigma ~ HalfNormal(1)

This is the design the original notebook intended ("each receptor its own rate constant, partially
pooled toward a global mean") without its integer-coded covariates. Values from one report on one
receptor (e.g. the 10 M3R k_act values of Masuho 2016) share a receptor effect, so they count as one
receptor rather than ten. A crossed study effect was tried first and dropped: with most receptors and
studies contributing a single value, receptor and study effects were not separately identifiable
(k_deg: R-hat 1.35, 36 divergences), while estimates for parameters where both models converged were
essentially unchanged. HalfNormal(2) priors left the smallest datasets (k_expt, n = 4) with hundreds
of divergences, so HalfNormal(1) is used (on the log scale this still allows ~7-fold spread per SD).

The receptor effect is non-centred and dropped when there is one receptor, or when every receptor
contributes a single value (it is then indistinguishable from the residual). The pooled ("global")
estimate is exp(mu): the rate for a typical receptor. For beta-arrestin, a second model gives each
isoform its own mu (mu_1, mu_2) with a shared receptor effect, from which P(k_ba2 > k_ba1) and the
ratio k_ba2/k_ba1 follow directly.

Reported statistic: posterior median of exp(mu) with the 95% highest-density interval of exp(mu).
The same statistic is drawn as the marker in every forest plot.

Sampler: PyMC 5.7.2 NUTS, 4 chains x 2000 draws after 2000 tuning steps, target_accept 0.99, seed 8927.

Data: data/OXT_All_Parameters_10162024.xlsx (identical to 01282026 apart from one Symbol name).
Receptor labels verified against the source papers are corrected before grouping (RELABEL below);
nothing else in the data is changed. k_ba values are converted from s^-1 to min^-1 (x60), as in every
other pipeline in this repository.

Units: k_P and k_ba are labelled as apparent first-order rate constants
(min^-1). Every curated value is ln(2)/t_1/2 from a single saturating-agonist assay, and no source
reports the GRK or beta-arrestin concentration needed for a bimolecular (M^-1 min^-1) constant.
Only the labels changed; no value or fit changed.

Run with the environment in code/python/environment_bhm.yml:
    python code/python/bhm_refit.py
"""

import re
import sys
from pathlib import Path

import arviz as az
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymc as pm

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "data" / "OXT_All_Parameters_10162024.xlsx"
OUT = REPO / "results" / "bhm_refit_2026-09-24"
SEED = 8927
SAMPLE = dict(draws=2000, tune=2000, chains=4, cores=4, target_accept=0.99, random_seed=SEED, progressbar=False)

# (Reference, workbook Receptor) -> verified receptor. Since 2026-09-24 the workbook itself carries the verified
# labels (checked against the PDFs), so these entries only matter for older copies of the workbook.
RELABEL = {
    ("Hein et al. 2006", "AT2R"): "A2AR",   # J Biol Chem 281:33345 studied A2A-adenosine and beta1-adrenergic receptors
    ("Hein et al. 2006", "β2AR"): "β1AR",
    ("Namkung et al. 2018", "AT2R"): "AT1R",  # Sci Signal 2018 profiled the angiotensin II type 1 receptor
}

# Forest-plot labels only: the abbreviations the manuscript defines. Data and result tables keep the workbook names.
DISPLAY = {"PGD2R": "DP", "GLP1": "GLP1R", "PTH": "PTHR"}

# name -> (Symbol in workbook, value column, multiplier, units, axis label)
PARAMS = {
    "Kd":        ("Kong, Koffg", "Kd", 1, "nM", "G protein equilibrium dissociation constant, $K_d$ (nM)"),
    "konG":      ("Kong, Koffg", "association rate", 1, "M$^{-1}$min$^{-1}$", "G protein association, $k_{on,G}$ (M$^{-1}$min$^{-1}$)"),
    "koffG":     ("Kong, Koffg", "dissociation rate", 1, "min$^{-1}$", "G protein dissociation, $k_{off,G}$ (min$^{-1}$)"),
    "kact":      ("Kact, Kinact", "association rate", 1, "min$^{-1}$", "G protein activation, $k_{act}$ (min$^{-1}$)"),
    "kinact":    ("Kact, Kinact", "dissociation rate", 1, "min$^{-1}$", "G protein inactivation, $k_{inact}$ (min$^{-1}$)"),
    "kP":        ("kp, kdp", "association rate", 1, "min$^{-1}$", "Phosphorylation, $k_P$ (min$^{-1}$)"),
    "kdp":       ("kp, kdp", "dissociation rate", 1, "min$^{-1}$", "Dephosphorylation, $k_{dp}$ (min$^{-1}$)"),
    "kba":       ("Kba", "association rate", 60, "min$^{-1}$", "$\\beta$-arrestin binding, $k_{ba}$ (min$^{-1}$)"),
    "kint":      ("Kint", "association rate", 1, "min$^{-1}$", "Internalization, $k_{int}$ (min$^{-1}$)"),
    "kcons_int": ("Kcons_int", "association rate", 1, "min$^{-1}$", "Constitutive internalization, $k_{cons,int}$ (min$^{-1}$)"),
    "krec":      ("Krec", "association rate", 1, "min$^{-1}$", "Recycling, $k_{rec}$ (min$^{-1}$)"),
    "kexpt":     ("Kins", "association rate", 1, "min$^{-1}$", "Export, $k_{expt}$ (min$^{-1}$)"),
    "krgt":      ("Krgt", "association rate", 1, "min$^{-1}$", "Retrograde transport, $k_{rgt}$ (min$^{-1}$)"),
    "kdeg":      ("Kdeg", "association rate", 1, "min$^{-1}$", "Degradation, $k_{deg}$ (min$^{-1}$)"),
}


def study_key(ref):
    """Same paper under different spellings ('Masuho et al., 2016' / 'Masuho et al. 2016') -> one study."""
    ref = str(ref).strip()
    year = re.search(r"(19|20)\d\d", ref)
    words = [w.strip(".,;") for w in re.split(r"[\s,&]+", ref)]
    tokens = [w for w in words if len(w) > 1 and w.lower() not in ("et", "al") and not re.match(r"\d", w)]
    return f"{tokens[0].lower()} {year.group(0) if year else ''}".strip()


def load(name):
    sym, col, mult, *_ = PARAMS[name]
    d = pd.read_excel(DATA)
    d.columns = d.columns.map(str.strip)
    x = d[(d.Symbol.astype(str) == sym) & d[col].notna()].copy()
    x["value"] = x[col].astype(str).str.split(n=1).str[0].astype(float) * mult
    x["Reference"] = x["Reference"].astype(str).str.strip()
    x["Receptor"] = [RELABEL.get((r, c), c) for r, c in zip(x.Reference, x.Receptor.astype(str).str.strip())]
    x["study"] = x.Reference.map(study_key)
    return x.reset_index(drop=True)


def build(x, isoform=None):
    y = np.log(x.value.to_numpy())
    receptors = list(dict.fromkeys(x.Receptor))
    studies = list(dict.fromkeys(x.study))
    r_idx = np.array([receptors.index(v) for v in x.Receptor])
    with pm.Model(coords={"receptor": receptors}) as model:
        if isoform is None:
            mu = pm.Normal("mu", mu=y.mean(), sigma=3)
            eta = mu
        else:
            mu = pm.Normal("mu", mu=y.mean(), sigma=3, shape=2)
            eta = mu[isoform]
        # the receptor term is kept only when some receptor holds more than one value; otherwise it is
        # indistinguishable from the residual and only makes the posterior harder to sample
        if 1 < len(receptors) < len(y):
            s_r = pm.HalfNormal("sigma_receptor", 1)
            z_r = pm.Normal("z_receptor", 0, 1, dims="receptor")
            a_r = pm.Deterministic("a_receptor", s_r * z_r, dims="receptor")
            eta = eta + a_r[r_idx]
        sigma = pm.HalfNormal("sigma", 1)
        pm.Normal("y", mu=eta, sigma=sigma, observed=y)
    return model, receptors, studies


def summarize(samples):
    """samples of log k -> (median, HDI low, HDI high) of k"""
    k = np.exp(np.asarray(samples).ravel())
    lo, hi = az.hdi(k, hdi_prob=0.95)
    return float(np.median(k)), float(lo), float(hi)


def diagnostics(name, idata):
    var = [v for v in ("mu", "sigma_receptor", "sigma") if v in idata.posterior]
    s = az.summary(idata, var_names=var, kind="diagnostics")
    return dict(parameter=name, max_rhat=float(s["r_hat"].max()), min_ess_bulk=float(s["ess_bulk"].min()),
                min_ess_tail=float(s["ess_tail"].min()), divergences=int(idata.sample_stats["diverging"].sum()))


def forest(name, rows, fname):
    """rows: list of (label, median, lo, hi, is_global)"""
    fig, ax = plt.subplots(figsize=(7.5, 0.42 * len(rows) + 1.1))
    for i, (lab, med, lo, hi, glob) in enumerate(rows[::-1]):
        c = "#1f3a93" if glob else "#2E86AB"
        ax.plot([lo, hi], [i, i], color=c, lw=2.2 if glob else 1.6)
        ax.plot(med, i, "D" if glob else "o", color=c, ms=8 if glob else 6, mec="white", mew=1)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([DISPLAY.get(r[0], r[0]) for r in rows[::-1]])
    for t, r in zip(ax.get_yticklabels(), rows[::-1]):
        if r[4]:
            t.set_fontweight("bold")
    ax.set_xlim(left=0)  # linear axis from zero, as in the original BHM forest plots
    ax.set_xlabel(PARAMS.get(name, PARAMS["kba"])[4] if name in PARAMS else name, fontsize=11)
    ax.grid(True, axis="x", alpha=0.25)
    ax.set_title("Posterior median and 95% HDI", fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(OUT / "figures" / fname, dpi=300)
    plt.close(fig)


def main():
    for sub in ("figures", "traces", "tables"):
        (OUT / sub).mkdir(parents=True, exist_ok=True)
    only = set(sys.argv[1:])
    summary, diag, receptor_rows = [], [], []

    for name in PARAMS:
        if only and name not in only:
            continue
        x = load(name)
        model, receptors, studies = build(x)
        with model:
            idata = pm.sample(**SAMPLE)
        idata.to_netcdf(OUT / "traces" / f"{name}.nc")
        post = idata.posterior
        g = summarize(post["mu"])
        summary.append(dict(Parameter=name, n=len(x), n_receptors=len(receptors), n_studies=len(studies),
                            median=g[0], hdi_lo=g[1], hdi_hi=g[2]))
        d = diagnostics(name, idata); diag.append(d)
        rows = []
        for r in receptors:
            a = post["a_receptor"].sel(receptor=r) if "a_receptor" in post else None
            if a is None:  # no receptor term: only the pooled estimate is identified
                continue
            m = summarize(post["mu"] + a)
            rows.append((r, *m, False))
            receptor_rows.append(dict(Parameter=name, receptor=r, median=m[0], hdi_lo=m[1], hdi_hi=m[2],
                                      n=int((x.Receptor == r).sum())))
        rows.append(("Global", *g, True))
        forest(name, rows, f"{name}_bhm_refit.png")
        print(f"{name}: n={len(x)} receptors={len(receptors)} studies={len(studies)} "
              f"median={g[0]:.3g} [{g[1]:.3g}, {g[2]:.3g}]  rhat={d['max_rhat']:.3f} div={d['divergences']}", flush=True)

    if not only or "kba_isoform" in only:
        x = load("kba")
        iso = np.where(x["Binding Partner"].astype(str).str.endswith("1"), 0, 1)
        model, receptors, studies = build(x, isoform=iso)
        with model:
            idata = pm.sample(**SAMPLE)
        idata.to_netcdf(OUT / "traces" / "kba_isoform.nc")
        mu = idata.posterior["mu"]
        rows = []
        for k, lab in ((0, "kba1"), (1, "kba2")):
            m = summarize(mu.isel(mu_dim_0=k))
            summary.append(dict(Parameter=lab, n=int((iso == k).sum()), n_receptors=len(receptors), n_studies=len(studies),
                                median=m[0], hdi_lo=m[1], hdi_hi=m[2]))
            rows.append(("β-arrestin 1" if k == 0 else "β-arrestin 2", *m, True))
        diff = (mu.isel(mu_dim_0=1) - mu.isel(mu_dim_0=0)).values.ravel()
        ratio = summarize(diff)
        p_gt = float((diff > 0).mean())
        summary.append(dict(Parameter="kba2/kba1", n=len(x), n_receptors=len(receptors), n_studies=len(studies),
                            median=ratio[0], hdi_lo=ratio[1], hdi_hi=ratio[2], p_kba2_gt_kba1=p_gt))
        d = diagnostics("kba_isoform", idata); diag.append(d)
        forest("kba", rows, "kba_isoform_bhm_refit.png")
        print(f"kba isoform: kba1 {rows[0][1]:.3g}, kba2 {rows[1][1]:.3g}, ratio {ratio[0]:.3g} [{ratio[1]:.3g}, {ratio[2]:.3g}], "
              f"P(kba2>kba1)={p_gt:.3f} rhat={d['max_rhat']:.3f} div={d['divergences']}", flush=True)

    tag = "_" + "_".join(sorted(only)) if only else ""
    pd.DataFrame(summary).to_csv(OUT / "tables" / f"bhm_refit_summary{tag}.csv", index=False)
    pd.DataFrame(diag).to_csv(OUT / "tables" / f"bhm_refit_diagnostics{tag}.csv", index=False)
    pd.DataFrame(receptor_rows).to_csv(OUT / "tables" / f"bhm_refit_receptor_estimates{tag}.csv", index=False)


def replot():
    """Redraw every forest plot from the saved result tables (no resampling)."""
    summ = pd.read_csv(OUT / "tables" / "bhm_refit_summary.csv").set_index("Parameter")
    rec = pd.read_csv(OUT / "tables" / "bhm_refit_receptor_estimates.csv")
    for name in PARAMS:
        rows = [(r.receptor, r.median, r.hdi_lo, r.hdi_hi, False) for r in rec[rec.Parameter == name].itertuples()]
        g = summ.loc[name]
        rows.append(("Global", g["median"], g.hdi_lo, g.hdi_hi, True))
        forest(name, rows, f"{name}_bhm_refit.png")
    rows = [(lab, summ.loc[k, "median"], summ.loc[k, "hdi_lo"], summ.loc[k, "hdi_hi"], True)
            for k, lab in (("kba1", "β-arrestin 1"), ("kba2", "β-arrestin 2"))]
    forest("kba", rows, "kba_isoform_bhm_refit.png")
    print("redrew", len(PARAMS) + 1, "forest plots")


if __name__ == "__main__":
    if sys.argv[1:] == ["--replot"]:
        replot()
    else:
        main()
