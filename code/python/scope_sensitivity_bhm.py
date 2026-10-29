"""
Table S5 -- scope sensitivity, BHM half, and the merged table.

Tests whether the pooled estimates depend on the non-Class A data: Class B and C GPCRs
(G-protein binding and activation) and non-GPCR proteins (trafficking and degradation). For every
parameter that contains such values, this script refits the receptor-level BHM of bhm_refit.py
(same model, priors, sampler settings and seed) twice: on the full set, and on Class A GPCRs only.
Receptor classes come from data/receptor_class_map.csv.

The full refit is done in the same run so the Class A-only comparison is like for like; it is also
checked against the published values in results/bhm_refit_2026-09-24/tables/bhm_refit_summary.csv.

It then merges the REM half (results/tables/Table_S5_scope_sensitivity_REM.csv, written by
code/R/table_s5_scope_sensitivity.R -- run that first) into
    results/tables/Table_S5_scope_sensitivity.csv / .xlsx
Traces and per-fit tables go to results/bhm_scope_sensitivity_2026-10-08/.

Headline values are not changed by this analysis; it is a sensitivity table only.

    Rscript code/R/table_s5_scope_sensitivity.R
    python code/python/scope_sensitivity_bhm.py
"""

import sys
from pathlib import Path

import arviz as az
import numpy as np
import pandas as pd
import pymc as pm

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bhm_refit as b  # noqa: E402  (same model, priors, SAMPLE settings and seed)

REPO = b.REPO
CLASS_MAP = REPO / "data" / "receptor_class_map.csv"
OUT = REPO / "results" / "bhm_scope_sensitivity_2026-10-08"
TABLES = REPO / "results" / "tables"
PUBLISHED = REPO / "results" / "bhm_refit_2026-09-24" / "tables" / "bhm_refit_summary.csv"
REM_CSV = TABLES / "Table_S5_scope_sensitivity_REM.csv"

TESTED = ["Kd", "konG", "koffG", "kact", "kinact", "kdeg", "krec", "kexpt", "krgt"]
UNITS = {"Kd": "nM", "konG": "M⁻¹min⁻¹"}


def classify(x, cmap):
    cls = x.Receptor.map(cmap)
    missing = sorted(set(x.Receptor[cls.isna()]))
    if missing:
        raise SystemExit(f"receptor(s) missing from {CLASS_MAP.name}: {missing}")
    return cls


def fit(x, trace):
    """Fit (or reload a finished fit from `trace`, so an interrupted run can resume)."""
    model, receptors, _ = b.build(x)
    if trace.exists():
        idata = az.from_netcdf(trace)
    else:
        with model:
            idata = pm.sample(**b.SAMPLE)
        idata.to_netcdf(trace)
    med, lo, hi = b.summarize(idata.posterior["mu"])
    d = b.diagnostics("", idata)
    return idata, dict(n=len(x), n_receptors=len(receptors), median=med, lo=lo, hi=hi,
                       max_rhat=d["max_rhat"], divergences=d["divergences"], min_ess_bulk=d["min_ess_bulk"])


SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def sci(v):
    """3 significant figures in the manuscript's style, e.g. 1.88×10⁻³"""
    e = int(np.floor(np.log10(abs(v))))
    m = round(v / 10 ** e, 2)
    if m >= 10:
        m, e = m / 10, e + 1
    return f"{m:.2f}×10{str(e).translate(SUP)}"


def fmt(m, lo, hi):
    return "—" if pd.isna(m) else f"{sci(m)} [{sci(lo)}, {sci(hi)}]"


def main():
    for sub in ("traces", "tables"):
        (OUT / sub).mkdir(parents=True, exist_ok=True)
    cmap = pd.read_csv(CLASS_MAP).set_index("receptor")["class"]
    published = pd.read_csv(PUBLISHED).set_index("Parameter")

    rows, row_classes = [], []
    for name in TESTED:
        x = b.load(name)
        x["class"] = classify(x, cmap)
        row_classes.append(x.assign(parameter=name)[["parameter", "Reference", "Receptor", "class", "value"]])
        excluded = x[x["class"] != "Class A"]
        xa = x[x["class"] == "Class A"].reset_index(drop=True)
        excl_txt = "; ".join(f"{r} ({c})" for r, c in excluded.groupby(["Receptor", "class"]).size().index) or "none"
        base = dict(parameter=name, excluded=excl_txt, n_excluded=len(excluded),
                    published=published.loc[name, "median"], published_lo=published.loc[name, "hdi_lo"],
                    published_hi=published.loc[name, "hdi_hi"])
        if excluded.empty:
            rows.append(dict(base, subset="full = Class A only", n=len(x), note="all values are Class A; nothing to exclude"))
            print(f"{name}: all {len(x)} values are Class A -- no refit needed", flush=True)
            continue
        _, full = fit(x, OUT / "traces" / f"{name}_full.nc")
        rows.append(dict(base, subset="full", **full))
        if len(xa) < 2:
            rows.append(dict(base, subset="Class A only", n=len(xa),
                             note=f"not computable: {len(xa)} Class A value{'s' if len(xa) != 1 else ''}"))
        else:
            _, a = fit(xa, OUT / "traces" / f"{name}_classA.nc")
            note = "n = 2; posterior is prior-dominated, interpret with caution" if len(xa) == 2 else ""
            rows.append(dict(base, subset="Class A only", note=note, **a))
        r = rows[-1]
        print(f"{name}: full {fmt(full['median'], full['lo'], full['hi'])} (n={full['n']}, rhat={full['max_rhat']:.3f}, "
              f"div={full['divergences']}) | Class A {fmt(r.get('median'), r.get('lo'), r.get('hi'))} (n={r['n']})", flush=True)

    bhm = pd.DataFrame(rows)
    bhm.to_csv(OUT / "tables" / "bhm_scope_sensitivity.csv", index=False)
    pd.concat(row_classes).to_csv(OUT / "tables" / "bhm_rows_by_class.csv", index=False)

    # reproduction check: full refit in this environment vs published refit
    chk = bhm[bhm.subset == "full"].copy()
    chk["rel_diff_vs_published"] = (chk["median"] / chk["published"] - 1).abs()
    chk[["parameter", "median", "published", "rel_diff_vs_published", "max_rhat", "divergences"]].to_csv(
        OUT / "tables" / "full_refit_vs_published.csv", index=False)
    print("\nfull refit vs published (relative difference of the median):")
    print(chk[["parameter", "median", "published", "rel_diff_vs_published"]].to_string(index=False))

    merge(bhm)


def merge(bhm):
    rem = pd.read_csv(REM_CSV) if REM_CSV.exists() else pd.DataFrame(columns=["parameter", "subset"])
    out = []
    for name in TESTED:
        b_rows = bhm[bhm.parameter == name]
        r_rows = rem[rem.parameter == name]
        bf = b_rows[b_rows.subset.str.startswith("full")].iloc[0]
        ba = b_rows[b_rows.subset != "full"].iloc[0] if (b_rows.subset != "full").any() else bf
        rec = dict(Parameter=name, Units=UNITS.get(name, "min⁻¹"), Values_excluded=bf.excluded)
        # BHM: reported value (published refit) and Class A-only refit
        rec["BHM full, published [95% HDI] (n)"] = f"{fmt(bf.published, bf.published_lo, bf.published_hi)} (n = {int(bf.n)})"
        if bf.n_excluded == 0:
            rec["BHM Class A only [95% HDI] (n)"] = "unchanged (all values Class A)"
            rec["BHM ratio (Class A ÷ full)"] = "1.00"
            rec["BHM Class A median inside full 95% HDI"] = "yes"
        elif pd.isna(ba.get("median")):
            rec["BHM Class A only [95% HDI] (n)"] = f"not computable (n = {int(ba.n)})"
            rec["BHM ratio (Class A ÷ full)"] = "—"
            rec["BHM Class A median inside full 95% HDI"] = "—"
        else:
            rec["BHM Class A only [95% HDI] (n)"] = f"{fmt(ba['median'], ba.lo, ba.hi)} (n = {int(ba.n)})"
            # ratio and interval check against the published (reported) full estimate; the same-run full refit
            # reproduces it to within 2% (full_refit_vs_published.csv)
            rec["BHM ratio (Class A ÷ full)"] = f"{ba['median'] / bf.published:.2f}"
            rec["BHM Class A median inside full 95% HDI"] = "yes" if bf.published_lo <= ba["median"] <= bf.published_hi else "no"
        # REM
        if r_rows.empty:
            rec["REM full [95% CI] (n)"] = "N/A (BHM only, Table S1)"
            rec["REM Class A only [95% CI] (n)"] = "N/A"
            rec["REM ratio (Class A ÷ full)"] = "—"
            rec["REM Class A mean inside full 95% CI"] = "—"
        else:
            rf = r_rows[r_rows.subset == "full"].iloc[0]
            ra = r_rows[r_rows.subset != "full"].iloc[0]
            rec["REM full [95% CI] (n)"] = f"{fmt(rf['mean'], rf.ci_lb, rf.ci_ub)} (n = {int(rf.n)})"
            if pd.isna(ra["mean"]):
                rec["REM Class A only [95% CI] (n)"] = f"not computable (n = {int(ra.n)})"
                rec["REM ratio (Class A ÷ full)"] = "—"
                rec["REM Class A mean inside full 95% CI"] = "—"
            else:
                rec["REM Class A only [95% CI] (n)"] = f"{fmt(ra['mean'], ra.ci_lb, ra.ci_ub)} (n = {int(ra.n)})"
                rec["REM ratio (Class A ÷ full)"] = f"{ra['mean'] / rf['mean']:.2f}"
                rec["REM Class A mean inside full 95% CI"] = "yes" if rf.ci_lb <= ra["mean"] <= rf.ci_ub else "no"
        notes = [n for n in (ba.get("note", ""), ) if isinstance(n, str) and n]
        rec["Note"] = "; ".join(notes)
        out.append(rec)
    table = pd.DataFrame(out)
    TABLES.mkdir(parents=True, exist_ok=True)
    table.to_csv(TABLES / "Table_S5_scope_sensitivity.csv", index=False)
    with pd.ExcelWriter(TABLES / "Table_S5_scope_sensitivity.xlsx") as xw:
        table.to_excel(xw, sheet_name="Table S5 (scope)", index=False)
        bhm.to_excel(xw, sheet_name="BHM fits", index=False)
        rem.to_excel(xw, sheet_name="REM fits", index=False)
    print("\nwrote", TABLES / "Table_S5_scope_sensitivity.xlsx")
    print(table.drop(columns=["Values_excluded"]).to_string(index=False))


if __name__ == "__main__":
    if sys.argv[1:] == ["--merge-only"]:
        merge(pd.read_csv(OUT / "tables" / "bhm_scope_sensitivity.csv"))
    else:
        main()
