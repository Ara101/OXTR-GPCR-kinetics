"""
Figures 4 and 5 -- BHM vs REM pooled rate constants with the curated values overlaid.

Standalone replacement for the Figure 4 / Figure 5 cells of 1_comparative_analysis_bhm_rem.ipynb.
Changes relative to those cells:

  * kP, kba1 and kba2 are apparent first-order rate constants (min^-1) and are
    plotted in Figure 4. They were derived as
    k = ln2 / t1/2 from saturating-agonist assays, and no source reports the
    GRK or beta-arrestin concentration a bimolecular constant would need.
    Figure 5 now holds konG only, the one genuinely second-order constant.
  * The y-axis limits are now computed from the data points as well as the
    intervals. The notebook set them from the intervals only (max upper x 2.5,
    min lower x 0.5), which silently dropped the M3R kact values
    (212, 342, 354 min^-1) and the lowest kdeg values (Vicentic 2000,
    ~1.0e-4; Shyng 1989, 4.3e-4 min^-1). matplotlib does not warn about this.
  * The jitter uses a fixed seed so the figure is reproducible.
  * Figure 5 shows konG (M^-1 min^-1) with a secondary axis giving the
    pseudo-first-order rate konG x [G] at an illustrative [G] = 10 uM,
    so it can be compared with Figure 4.

Pooled values are read unchanged from data/GPCR_parameters.xlsx
(sheet GPCR_parameters_beta12_geomean) -- no model is refit here.
The individual G-protein values (konG, koffG) are plotted as curated; unit
conversions (s^-1 -> min^-1) are recorded in the Note column of the dataset.

Run from anywhere:  python code/python/figure4_figure5_rate_constants.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

REPO = Path(__file__).resolve().parents[2]
DATA_POINTS = REPO / "data" / "OXT_All_Parameters_01282026.xlsx"
SUMMARY = REPO / "data" / "GPCR_parameters.xlsx"
OUT = REPO / "results" / "figures" / "python"

BHM_COLOR = "#2E86AB"
REM_COLOR = "#A23B72"
POINT_COLOR = "#656565"
GRID_COLOR = "#E0E0E0"
RNG = np.random.default_rng(20260923)

FIRST_ORDER = [
    "koffG", "kact", "kinact", "kP", "kdp", "kba 1", "kba 2",
    "kcons_int", "kint", "krec", "kexpt", "krgt", "kdeg",
]
SECOND_ORDER = ["konG"]
G_PROTEIN_M = 10e-6  # illustrative [G] for the pseudo-first-order axis of Figure 5
LABELS = {
    "konG": r"$k_{\mathrm{on,G}}$",
    "koffG": r"$k_{\mathrm{off,G}}$",
    "kact": r"$k_{\mathrm{act}}$",
    "kinact": r"$k_{\mathrm{inact}}$",
    "kP": r"$k_{\mathrm{P}}$",
    "kdp": r"$k_{\mathrm{dp}}$",
    "kba 1": r"$k_{\mathrm{ba1}}$",
    "kba 2": r"$k_{\mathrm{ba2}}$",
    "kcons_int": r"$k_{\mathrm{cons,int}}$",
    "kint": r"$k_{\mathrm{int}}$",
    "krec": r"$k_{\mathrm{rec}}$",
    "kexpt": r"$k_{\mathrm{expt}}$",
    "krgt": r"$k_{\mathrm{rgt}}$",
    "kdeg": r"$k_{\mathrm{deg}}$",
}

# workbook Symbol token -> (summary Parameter name, value column)
SYMBOL_MAP = {
    "Kong": ("konG", "association rate"),
    "Koffg": ("koffG", "dissociation rate"),
    "Kact": ("kact", "association rate"),
    "Kinact": ("kinact", "dissociation rate"),
    "kp": ("kP", "association rate"),
    "kdp": ("kdp", "dissociation rate"),
    "Kba": ("kba", "association rate"),
    "Kint": ("kint", "association rate"),
    "Kcons_int": ("kcons_int", "association rate"),
    "Krec": ("krec", "association rate"),
    "Kexpt": ("kexpt", "association rate"),
    "Krgt": ("krgt", "association rate"),
    "Kdeg": ("kdeg", "association rate"),
}


def leading_number(val):
    try:
        return float(str(val).split(" ± ")[0])
    except ValueError:
        return np.nan


def load_points():
    raw = pd.read_excel(DATA_POINTS, sheet_name="all_params")
    rows = []
    for _, r in raw.iterrows():
        for token in str(r["Symbol"]).split(","):
            token = token.strip()
            if token not in SYMBOL_MAP:
                continue
            param, col = SYMBOL_MAP[token]
            value = leading_number(r[col])
            if np.isnan(value):
                continue
            if param == "kba":
                value *= 60  # curated in s^-1; analysed in min^-1 (same x60 as the R and BHM code)
                partner = str(r["Binding Partner"])
                param = "kba 1" if partner.endswith("1") else "kba 2"
            rows.append({"Parameter": param, "Receptor": r["Receptor"], "Mean": value})
    return pd.DataFrame(rows)


def load_summary():
    s = pd.read_excel(SUMMARY, sheet_name="GPCR_parameters_beta12_geomean")
    return s.set_index("Parameter")


def draw_panel(ax, params, summary, points, jitter_width=0.04):
    x = np.arange(len(params))
    width = 0.35
    s = summary.loc[params]

    ax.errorbar(
        x - width / 2, s["BHM (mean)"],
        yerr=[s["BHM (mean)"] - s["BHM (lower [95% CI])"], s["BHM (upper [95% CI])"] - s["BHM (mean)"]],
        fmt="D", color=BHM_COLOR, ecolor=BHM_COLOR, elinewidth=1.5, capsize=5, capthick=1.5,
        markersize=9, markeredgewidth=1.5, markeredgecolor="white", alpha=0.9, zorder=3,
    )
    has_rem = s["REM (mean)"].notna().to_numpy()
    if has_rem.any():
        r = s[has_rem]
        ax.errorbar(
            x[has_rem] + width / 2, r["REM (mean)"],
            yerr=[r["REM (mean)"] - r["REM (lower [95% CI])"], r["REM (upper [95% CI])"] - r["REM (mean)"]],
            fmt="s", color=REM_COLOR, ecolor=REM_COLOR, elinewidth=1.5, capsize=5, capthick=1.5,
            markersize=8, markeredgewidth=1.5, markeredgecolor="white", alpha=0.9, zorder=3,
        )

    plotted = []
    for i, p in enumerate(params):
        vals = points.loc[points["Parameter"] == p, "Mean"].to_numpy()
        ax.scatter(RNG.normal(i, jitter_width, size=len(vals)), vals, alpha=0.45, s=40,
                   color=POINT_COLOR, edgecolors="none", zorder=2)
        plotted.append(vals)

    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[p] for p in params], fontsize=18, rotation=45, ha="right")
    ax.set_xlim(-0.7, len(params) - 0.3)

    # limits cover every interval AND every data point, padded by ~1/3 decade
    all_vals = np.concatenate(plotted + [
        s["BHM (lower [95% CI])"].to_numpy(), s["BHM (upper [95% CI])"].to_numpy(),
        s["REM (lower [95% CI])"].dropna().to_numpy(), s["REM (upper [95% CI])"].dropna().to_numpy(),
    ])
    lo, hi = all_vals.min(), all_vals.max()
    ax.set_ylim(10 ** (np.floor(np.log10(lo) * 3) / 3 - 0.2), 10 ** (np.ceil(np.log10(hi) * 3) / 3 + 0.2))

    ax.grid(True, alpha=0.4, linestyle="-", linewidth=0.7, color=GRID_COLOR, which="both", axis="y")
    ax.grid(True, alpha=0.2, linestyle="-", linewidth=0.5, color=GRID_COLOR, which="major", axis="x")
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)
        spine.set_color("gray")
        spine.set_alpha(0.8)

    n_plotted = sum(len(v) for v in plotted)
    n_in_view = sum(((v >= ax.get_ylim()[0]) & (v <= ax.get_ylim()[1])).sum() for v in plotted)
    return n_plotted, n_in_view


def legend(ax, loc="upper right", fontsize=16, include_rem=True):
    handles = [
        Line2D([0], [0], marker="D", color="w", markerfacecolor=BHM_COLOR, markeredgecolor="white",
               markersize=10, markeredgewidth=1.5, label="BHM (95% CrI)"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor=REM_COLOR, markeredgecolor="white",
               markersize=9, markeredgewidth=1.5, label="REM (95% CI)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=POINT_COLOR, markeredgecolor="none",
               markersize=6, label="Individual data points"),
    ]
    if not include_rem:
        handles = [h for h in handles if not h.get_label().startswith("REM")]
    ax.legend(handles=handles, loc=loc, fontsize=fontsize, framealpha=0.95, edgecolor="gray",
              fancybox=True, shadow=True)


def main():
    plt.style.use("seaborn-v0_8-whitegrid")
    OUT.mkdir(parents=True, exist_ok=True)
    points = load_points()
    summary = load_summary()
    report = []

    # ---- Figure 4: first-order and apparent first-order rate constants -------
    fig, ax = plt.subplots(figsize=(16, 8))
    n, n_in = draw_panel(ax, FIRST_ORDER, summary, points)
    ax.set_ylabel("Rate constant (min$^{-1}$), log scale", weight="semibold", fontsize=18, labelpad=10)
    ax.set_xlabel("First-order rate constants", weight="semibold", fontsize=18, labelpad=10)
    legend(ax, loc="lower left")
    fig.tight_layout()
    fig.savefig(OUT / "Figure4_first_order.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    report.append(("Figure4", n, n_in))

    # ---- Figure 5: the second-order constant konG, with a pseudo-first-order axis --
    fig, ax = plt.subplots(figsize=(5.5, 8))
    n, n_in = draw_panel(ax, SECOND_ORDER, summary, points, jitter_width=0.06)
    ax.set_xlim(-0.8, 0.8)
    ax.set_ylabel("Rate constant (M$^{-1}$min$^{-1}$), log scale", weight="semibold", fontsize=16, labelpad=8)
    sec = ax.secondary_yaxis("right", functions=(lambda k: k * G_PROTEIN_M, lambda r: r / G_PROTEIN_M))
    sec.set_ylabel("Pseudo-first-order rate at 10 µM G protein (min$^{-1}$)", fontsize=13, labelpad=8)
    legend(ax, loc="lower right", fontsize=11, include_rem=False)  # konG is BHM-only
    fig.savefig(OUT / "Figure5_second_order.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    report.append(("Figure5", n, n_in))

    print("data points plotted / inside axis limits:")
    for name, n, n_in in report:
        print(f"  {name:18s} {n:3d} / {n_in:3d}")
    print("points per parameter:", points.groupby("Parameter").size().to_dict())


if __name__ == "__main__":
    main()
