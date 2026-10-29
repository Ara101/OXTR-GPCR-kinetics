"""
Figure 3: ligand and receptor composition of the curated kinetic values (pie charts).

Reads the 'ligands' sheet of data/GPCR_parameters.xlsx (one row per curated kinetic value used in the pooled
parameters: Parameter, receptor, ligand) and draws one pie per column. Categories below 2% of values are grouped
into "Other". Same layout as the original notebook helper (1_comparative_analysis_bhm_rem.ipynb,
create_dual_chart_with_arrows, pie only).

Outputs
  results/comparative_analysis/Pie chart quantification of Ligands within the Stud.png
  results/comparative_analysis/Pie chart quantification of Receptors within the Stud.png
  results/comparative_analysis/Figure3_ligands_receptors.jpg   (both pies side by side, Figure 3)
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "data" / "GPCR_parameters.xlsx"
OUT = REPO / "results" / "comparative_analysis"
REPORT = OUT / "Figure3_ligands_receptors.jpg"
THRESHOLD = 2  # percent


def pie(series, title, path):
    counts = series.astype(str).str.strip().value_counts()
    total = counts.sum()
    pct = counts / total * 100
    main, other = counts[pct >= THRESHOLD], counts[pct < THRESHOLD]
    data = pd.concat([main, pd.Series({"Other": other.sum()})]) if len(other) else counts
    labels = []
    for name, n in data.items():
        p = n / total * 100
        if name == "Other" and len(other):
            labels.append(f"Other\n({len(other)} items, {p:.1f}%)")
        else:
            shown = name if len(name) <= 15 else name[:12] + "..."
            labels.append(f"{shown}\n{n} ({p:.1f}%)")
    fig, ax = plt.subplots(figsize=(9, 7))
    explode = [0.05] * len(data)
    if len(other):
        explode[-1] = 0.1
    wedges, _ = ax.pie(data.values, colors=plt.cm.Set3(np.arange(len(data))), explode=explode, startangle=90,
                       wedgeprops=dict(edgecolor="w", linewidth=1))
    box = dict(boxstyle="round,pad=0.3", fc="w", ec="black", lw=0.72)
    kw = dict(arrowprops=dict(arrowstyle="-"), bbox=box, zorder=0, va="center", fontsize=10, weight="bold")
    for w, lab in zip(wedges, labels):
        ang = (w.theta2 + w.theta1) / 2.0
        y, x = np.sin(np.deg2rad(ang)), np.cos(np.deg2rad(ang))
        kw["arrowprops"].update({"connectionstyle": f"angle,angleA=0,angleB={ang}"})
        ax.annotate(lab, xy=(x, y), xytext=(1.3 * x, 1.3 * y), horizontalalignment={-1: "right", 1: "left"}[int(np.sign(x))], **kw)
    ax.set_title(f"{title}\n(Threshold: {THRESHOLD}%)", fontsize=14)
    ax.set_aspect("equal")
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return total


def main():
    d = pd.read_excel(DATA, sheet_name="ligands")
    OUT.mkdir(parents=True, exist_ok=True)
    lig = OUT / "Pie chart quantification of Ligands within the Stud.png"
    rec = OUT / "Pie chart quantification of Receptors within the Stud.png"
    n = pie(d["ligand"], "Ligand Distribution", lig)
    pie(d["receptor"], "Receptor Distribution", rec)
    # Figure 3: ligands (A) left, receptors (B) right, on a white 16:9 canvas
    a, b = Image.open(lig).convert("RGB"), Image.open(rec).convert("RGB")
    h = max(a.height, b.height)
    W = int(round((a.width + b.width) * 1.04))
    H = max(h, int(W * 9 / 16))
    canvas = Image.new("RGB", (W, H), "white")
    gap = (W - a.width - b.width) // 3
    canvas.paste(a, (gap, (H - a.height) // 2))
    canvas.paste(b, (2 * gap + a.width, (H - b.height) // 2))
    canvas.save(REPORT, "JPEG", quality=92, optimize=True)
    print(f"Figure 3 written from {n} curated values")


if __name__ == "__main__":
    main()
