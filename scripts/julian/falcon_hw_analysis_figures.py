"""Figuras de análisis de las corridas en hardware real IBM (chunked-QAOA).

Lee `results/julian/hw_analysis_<instance>.json` (de `falcon_hw_retrieve.py`) + `runs_summary.csv`.
Genera 3 PNGs en results/julian/figures/:
  figHW1 - ruido / rendimiento de factibilidad (fracción one-hot válida en HW vs 100% sim + #factibles)
  figHW2 - distribución de SRS de muestras factibles por bloque (marca mejor post-selección + óptimo)
  figHW3 - comparación de SRS: HW-chunked vs DP/histórico/sim-chunked por instancia (infactible marcado)

Correr: .venv-quantum/bin/python scripts/julian/falcon_hw_analysis_figures.py
"""
from __future__ import annotations

import csv
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
JULIAN = os.path.join(ROOT, "results", "julian")
FIGDIR = os.path.join(JULIAN, "figures")
os.makedirs(FIGDIR, exist_ok=True)


def _load_analyses():
    out = {}
    for f in sorted(glob.glob(os.path.join(JULIAN, "hw_analysis_*.json"))):
        d = json.load(open(f))
        out[d["instance"]] = d
    return out


def _save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, name), dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("escrito", os.path.join("results/julian/figures", name))


# --------------------------------------------- figHW1: ruido / feasibility yield
def fig_noise(analyses):
    labels, blocks = [], []
    for inst, d in analyses.items():
        for b in d["per_block"]:
            labels.append(f"{inst}\nb{b['block_idx']} ({b['n_qubits']}q)")
            blocks.append(b)
    x = np.arange(len(blocks))
    oh = [b["one_hot_frac_hw"] for b in blocks]
    nfeas = [b["n_feasible_samples"] for b in blocks]

    fig, ax1 = plt.subplots(figsize=(max(7, 1.3 * len(blocks)), 5))
    ax1.bar(x, oh, 0.6, color="tab:blue", label="fracción one-hot válida en HW")
    ax1.axhline(1.0, color="tab:green", ls="--", lw=2, label="sim ideal (100%)")
    ax1.set_ylabel("Fracción de shots one-hot válida", color="tab:blue")
    ax1.set_ylim(0, 1.05)
    ax1.set_xticks(x); ax1.set_xticklabels(labels, fontsize=8)
    for xi, v in zip(x, oh):
        ax1.annotate(f"{v:.0%}", (xi, v), ha="center", va="bottom", fontsize=8, color="tab:blue")
    ax2 = ax1.twinx()
    ax2.plot(x, nfeas, "o-", color="tab:red", label="# muestras factibles (post-selección)")
    ax2.set_ylabel("# muestras factibles distintas", color="tab:red")
    ax1.set_title("Ruido de hardware vs rendimiento de factibilidad (chunked-QAOA XY, ibm_kingston)\n"
                  "el ruido tira ~⅔ de los shots fuera del subespacio, pero XY+post-selección igual halla factibles")
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper right")
    _save(fig, "figHW1_noise_feasibility.png")


# --------------------------------------------- figHW2: distribución de SRS factibles
def fig_distribution(analyses):
    # elegir la instancia con más muestras factibles totales
    inst = max(analyses, key=lambda k: sum(b["n_feasible_samples"] for b in analyses[k]["per_block"]))
    d = analyses[inst]
    nb = len(d["per_block"])
    fig, axes = plt.subplots(1, nb, figsize=(4.2 * nb, 4), sharey=True)
    if nb == 1:
        axes = [axes]
    for ax, b in zip(axes, d["per_block"]):
        dist = b["feasible_srs_dist"]                       # lista de (SRS, count)
        if not dist:
            ax.set_title(f"bloque {b['block_idx']}: sin factibles"); continue
        srs_vals = np.array([s for s, _ in dist])
        counts = np.array([c for _, c in dist])
        ax.hist(srs_vals, bins=20, weights=counts, color="tab:blue", alpha=0.75,
                label="muestras factibles (por shots)")
        ax.axvline(b["best_SRS"], color="tab:red", lw=2, label=f"mejor post-sel {b['best_SRS']:.3f}")
        ax.set_title(f"{inst} bloque {b['block_idx']} ({b['n_qubits']}q)\none-hot HW={b['one_hot_frac_hw']:.0%}")
        ax.set_xlabel("SRS del bloque (más alto = mejor)")
        ax.legend(fontsize=7)
    axes[0].set_ylabel("shots (peso)")
    fig.suptitle("Distribución de las soluciones FACTIBLES muestreadas en hardware (por bloque)\n"
                 "el ruido dispersa; la post-selección se queda con la de mejor energía factible", y=1.04)
    _save(fig, "figHW2_sampled_distribution.png")


# --------------------------------------------- figHW3: SRS HW vs baselines
def fig_srs_compare(analyses):
    rows = list(csv.DictReader(open(os.path.join(ROOT, "results", "runs_summary.csv"))))

    def pick(method, inst):
        cands = [r for r in rows if r["method"] == method and r["instance"] == inst
                 and r["SRS"] not in ("", None) and not r.get("window_label")]
        if not cands:
            return None
        best = max(cands, key=lambda r: float(r["SRS"]))
        return {"SRS": float(best["SRS"]), "feasible": best["feasible"] == "True"}

    insts = list(analyses.keys()) or ["small"]
    series = [("historical", "histórico u=0"), ("dp", "DP exacto"),
              ("qaoa_chunked", "QAOA-chunked sim"), ("qaoa_chunked_hardware", "QAOA-chunked HW")]
    fig, ax = plt.subplots(figsize=(max(7, 2.2 * len(insts)), 5))
    nM = len(series); width = 0.8 / nM
    x = np.arange(len(insts))
    colors = ["tab:gray", "tab:green", "tab:cyan", "tab:red"]
    for j, (m, lbl) in enumerate(series):
        xs = x + (j - nM / 2) * width + width / 2
        first = True
        for xi, inst in zip(xs, insts):
            r = pick(m, inst)
            if not r:
                continue
            ax.bar(xi, r["SRS"], width, color=colors[j], edgecolor="black",
                   hatch="" if r["feasible"] else "xx", label=lbl if first else None)
            first = False
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x); ax.set_xticklabels([f"{i}" for i in insts])
    ax.set_ylabel("SRS (más alto = mejor; negativo)")
    ax.set_title("SRS: QAOA-chunked en HARDWARE real vs baselines (hatch ✗ = infactible)\n"
                 "small T12/L3: HW alcanza el óptimo DP y es factible (XY + post-selección)")
    ax.legend(fontsize=8, loc="lower right")
    _save(fig, "figHW3_srs_hardware_vs_baselines.png")


def main():
    analyses = _load_analyses()
    if not analyses:
        print("No hay hw_analysis_*.json; correr falcon_hw_retrieve.py primero.")
        return
    print("instancias con análisis HW:", list(analyses))
    fig_noise(analyses)
    fig_distribution(analyses)
    fig_srs_compare(analyses)
    print("Listo:", FIGDIR)


if __name__ == "__main__":
    main()
