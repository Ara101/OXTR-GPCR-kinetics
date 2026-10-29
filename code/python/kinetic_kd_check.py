"""
Kinetic K_d check (Results: G-protein association, dissociation and binding affinity).

Four of the five curated GPCR:G-protein K_d values are k_off/k_on ratios for the same receptor-G protein
pair (Hein et al. 2006, SR17; Ilyaskina et al. 2018, SR20); the fifth (Namkung et al. 2018, SR37) is
not. This script compares the pooled K_d with (i) those matched-pair ratios and (ii) the ratio of the
pooled k_off,G and k_on,G from the published BHM refits. The two pooled rate constants come from
independent fits over different receptor sets (k_off,G: 14 values; k_on,G: 4), so their posterior
draws are paired at random and that interval is approximate.

    python code/python/kinetic_kd_check.py   ->  results/tables/kinetic_Kd_check.csv
"""

from pathlib import Path

import arviz as az
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "data" / "OXT_All_Parameters_10162024.xlsx"
TRACES = REPO / "results" / "bhm_refit_2026-09-24" / "traces"
OUT = REPO / "results" / "tables" / "kinetic_Kd_check.csv"


def draws(name):
    return np.exp(az.from_netcdf(TRACES / f"{name}.nc").posterior["mu"].values.ravel())


def main():
    rng = np.random.default_rng(8927)
    kon, koff, kd = draws("konG"), draws("koffG"), draws("Kd")
    kinetic = koff / rng.permutation(kon) * 1e9  # min^-1 / (M^-1 min^-1) = M -> nM
    k_lo, k_hi = az.hdi(kinetic, hdi_prob=0.95)
    d_lo, d_hi = az.hdi(kd, hdi_prob=0.95)
    rows = [
        dict(quantity="K_d, pooled (BHM, n = 5)", median_nM=np.median(kd), hdi_lo_nM=d_lo, hdi_hi_nM=d_hi),
        dict(quantity="k_off,G / k_on,G, ratio of posterior medians", median_nM=np.median(koff) / np.median(kon) * 1e9,
             hdi_lo_nM=np.nan, hdi_hi_nM=np.nan),
        dict(quantity="k_off,G / k_on,G, ratio of randomly paired draws", median_nM=np.median(kinetic),
             hdi_lo_nM=k_lo, hdi_hi_nM=k_hi),
    ]
    d = pd.read_excel(DATA)
    d = d[d.Symbol.astype(str) == "Kong, Koffg"].dropna(subset=["association rate", "dissociation rate"])
    first = lambda v: float(str(v).split()[0])
    for _, r in d.iterrows():
        rows.append(dict(quantity=f"matched pair: {r.Reference.strip()}, {r.Receptor}-{r['Binding Partner']} "
                                  f"(curated K_d {first(r.Kd):.0f} nM)" if pd.notna(r.Kd) else f"matched pair: {r.Reference}",
                         median_nM=first(r["dissociation rate"]) / first(r["association rate"]) * 1e9,
                         hdi_lo_nM=np.nan, hdi_hi_nM=np.nan))
    out = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)
    print(out.round(0).to_string(index=False))


if __name__ == "__main__":
    main()
