"""Figuras de escalamiento / ventaja eventual / feasibilidad del enfoque cuántico (Falcon).

Genera PNGs en results/julian/figures/ para el pptx. Marco HONESTO (spec §7): a 1 embalse el DP es
exacto y polinomial -> no hay ventaja cuántica; la ventaja eventual es en el régimen multi-embalse
(estado DP exponencial en #embalses R, qubits lineal). Las figuras de escalamiento comparan CADA
enfoque en tiempo/qubits/memoria; las de feasibilidad usan datos reales de results/runs_summary.csv.

Correr: .venv-quantum/bin/python scripts/julian/falcon_advantage_figures.py
"""
from __future__ import annotations

import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FIGDIR = os.path.join(ROOT, "results", "julian", "figures")
os.makedirs(FIGDIR, exist_ok=True)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

L = 5                     # niveles oficiales
BLK = 5                   # weeks/chunk (medium/large) -> qubits/chunk = BLK*L = 25
SV_WALL = 30              # ~ límite statevector exacto (2^30*16B ≈ 16 GB)
RAM = 16 * 2**30          # 16 GB en bytes


def _save(fig, name):
    path = os.path.join(FIGDIR, name)
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("escrito", os.path.relpath(path, ROOT))


# ---------------------------------------------------------------- Fig A: tiempo vs T
def fig_time():
    T = np.arange(6, 105)
    brute = L ** T.astype(float)                       # O(L^T)
    dp = 2.0 * T**2 * L**2                             # O(T^2 L^2) transiciones
    qaoa_full = (2.0 ** (T * L)) * 100                 # statevector: 2^{TL} * iters
    qaoa_chunk = np.ceil(T / BLK) * (2.0 ** (BLK * L)) * 100   # (T/b) 2^{bL} iters -> lineal en T
    qaoa_hw = (T * L) ** 2 * 4096                      # HW: ~2q-gates * shots (poli, sin 2^n)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.semilogy(T, brute, "k--", label=r"Brute force  $O(L^T)$")
    ax.semilogy(T, qaoa_full, color="tab:red", label=r"QAOA sim full  $O(2^{TL}\cdot it)$")
    ax.semilogy(T, dp, color="tab:green", lw=2.5, label=r"DP exacto  $O(T^2L^2)$")
    ax.semilogy(T, qaoa_chunk, color="tab:blue", lw=2.5,
                label=r"QAOA sim chunked  $O(\frac{T}{b}2^{bL})$ (lineal en T)")
    ax.semilogy(T, qaoa_hw, color="tab:purple", lw=2.5,
                label=r"QAOA hardware  $O(\mathrm{depth}\cdot shots)$ (poli)")
    # marcadores DP medidos (docs/ANALISIS §3): T26/L5=30024, T52/L5=116676 transiciones
    ax.scatter([26, 52], [30024, 116676], color="tab:green", zorder=5, s=40,
               label="DP transiciones (medido)")
    for x in (12, 26, 52):
        ax.axvline(x, color="gray", ls=":", alpha=0.4)
    ax.set_xlabel("Horizonte T (semanas)")
    ax.set_ylabel("Operaciones (escala relativa, log)")
    ax.set_title("Complejidad temporal por enfoque (L=5)\nchunked y hardware evitan el muro exponencial de la simulación")
    ax.legend(fontsize=8, loc="upper left")
    ax.set_ylim(1, 1e60)
    _save(fig, "figA_time_scaling.png")


# ---------------------------------------------------------------- Fig B: memoria vs T
def fig_memory():
    T = np.arange(6, 105)
    sv_full = (2.0 ** (T * L)) * 16                    # statevector full
    sv_chunk = np.full_like(T, 2.0 ** (BLK * L) * 16, dtype=float)  # constante
    dp_mem = 2.0 * T**2 * L * 32                       # estados * ~32 B

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.semilogy(T, sv_full, color="tab:red", label=r"Statevector full  $2^{TL}\cdot16$B")
    ax.semilogy(T, sv_chunk, color="tab:blue", lw=2.5,
                label=r"Statevector chunked  $2^{bL}\cdot16$B (constante)")
    ax.semilogy(T, dp_mem, color="tab:green", lw=2.5, label=r"DP estados  $\sim T^2L\cdot32$B")
    ax.axhline(RAM, color="black", ls="--", label="16 GB (laptop M4 / T4 GPU)")
    ax.set_xlabel("Horizonte T (semanas)")
    ax.set_ylabel("Memoria (bytes, log)")
    ax.set_title("Memoria por enfoque (L=5)\nla sim full explota ~30 qubits; chunked queda plano; DP minúsculo")
    ax.legend(fontsize=8, loc="upper left")
    ax.set_ylim(1e3, 1e80)
    _save(fig, "figB_memory_scaling.png")


# ---------------------------------------------------------------- Fig 1: qubits vs T
def fig_qubits():
    T = np.arange(6, 105)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(T, T * L, color="tab:red", label=f"one-hot  T·L (L={L})")
    ax.plot(T, T * int(np.ceil(np.log2(L))), color="tab:orange", label="binary  T·⌈log₂L⌉")
    ax.plot(T, T * (L - 1), color="tab:brown", label="domain-wall  T·(L−1)")
    ax.axhline(BLK * L, color="tab:blue", lw=2.5, label=f"chunked = b·L = {BLK*L} (constante)")
    ax.axhspan(SV_WALL, ax.get_ylim()[1] if False else 600, color="red", alpha=0.06)
    ax.axhline(SV_WALL, color="black", ls="--", label="muro statevector ~30 qubits")
    ax.set_xlabel("Horizonte T (semanas)")
    ax.set_ylabel("Qubits (circuito full-run)")
    ax.set_title("Qubits vs horizonte: full-run crece lineal con T; chunking lo acota\n"
                 "→ QPU de hoy corre 'en batches' cualquier T; full-run largo necesita FT-QPU")
    ax.set_ylim(0, 560)
    ax.legend(fontsize=8, loc="upper left")
    _save(fig, "fig1_qubits_vs_T.png")


# ---------------------------------------------------------------- Fig 2: DP vs brute (1 embalse)
def fig_dp_vs_brute():
    T = np.arange(4, 60)
    brute = L ** T.astype(float)
    dp_states = 2.0 * T**2 * L
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.semilogy(T, brute, "k--", label=r"candidatos $L^T$ (fuerza bruta)")
    ax.semilogy(T, dp_states, color="tab:green", lw=2.5, label=r"estados DP $\sim T^2L$ (polinomial)")
    ax.scatter([12, 26, 52], [433, 6505, 25086], color="tab:green", zorder=5, s=45,
               label="DP estados (medido, docs §3)")
    ax.set_xlabel("Horizonte T (semanas)")
    ax.set_ylabel("Estados / candidatos (log)")
    ax.set_title("Un solo embalse: el DP colapsa la explosión $L^T$ → exacto y polinomial\n"
                 "(no hay ventaja cuántica en el problema oficial)")
    ax.legend(fontsize=9, loc="center right")
    ax.set_ylim(1, 1e40)
    _save(fig, "fig2_dp_vs_brute.png")


# ---------------------------------------------------------------- Fig 3: crossover multi-embalse
def fig_multireservoir(Tref=26):
    R = np.arange(1, 13)
    dp_states = ((4 * Tref + 1) * L) ** R.astype(float)     # estado conjunto de R embalses
    qubits_full = R * Tref * L
    qubits_chunk = R * BLK * L
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.semilogy(R, dp_states, color="tab:green", lw=2.5,
                label=r"estados DP $\sim ((4T{+}1)L)^R$ (exp. en R)")
    ax.semilogy(R, qubits_full, color="tab:red", lw=2.5, label="qubits QUBO full  R·T·L")
    ax.semilogy(R, qubits_chunk, color="tab:blue", lw=2.5, label="qubits chunked  R·b·L")
    ax.axhline(1e12, color="black", ls="--", label="muro clásico tratable (~1e12 estados)")
    ax.axhline(1000, color="gray", ls=":", label="≈ qubits FT-QPU (proyección)")
    # crossover R*: primer R donde dp_states supera 1e12
    Rstar = int(R[np.argmax(dp_states > 1e12)])
    ax.axvline(Rstar, color="purple", ls="-.", alpha=0.7)
    ax.text(Rstar + 0.1, 1e14, f"R*≈{Rstar}\nDP intratable,\nqubits chicos", color="purple", fontsize=8)
    ax.set_xlabel("Nº de embalses acoplados R")
    ax.set_ylabel("Estados (DP) / qubits (QUBO), log")
    ax.set_title(f"VENTAJA EVENTUAL (proyección, T={Tref}, L={L}): multi-embalse\n"
                 "DP exponencial en R; qubits lineal → una (FT-)QPU corre lo que el DP no puede")
    ax.legend(fontsize=8, loc="upper left")
    ax.set_ylim(1, 1e40)
    fig.text(0.5, -0.02, "Extensión fuera del benchmark oficial (1 embalse). Supuesto: DP conjunto "
             "con curse-of-dimensionality (refs [1-3] del spec).", ha="center", fontsize=7, style="italic")
    _save(fig, "fig3_multireservoir_crossover.png")


# ---------------------------------------------------------------- Fig C: depth vs qubits (medido)
def fig_depth():
    import falcon_config as fcfg
    import falcon_constants as fc
    import falcon_data as fd
    import falcon_qaoa as qa
    import falcon_qaoa_hardware as hw
    import falcon_qubo as fq
    from qiskit import transpile
    from qiskit_ibm_runtime.fake_provider import FakeBrisbane

    df = fd.build_weekly_benchmark(write=False)
    c = fc.load_official_constants()
    Rall = df["R_obs_m3_week"].to_numpy()
    dSall = df["DeltaS_obs_m3"].to_numpy()
    S0 = df.attrs["S0_m3"]
    fake = FakeBrisbane()

    ns, depths, twoq = [], [], []
    for Tb in (2, 3, 4, 5):                             # bloques L=5 -> n = 10,15,20,25
        R, dS = Rall[:Tb], dSall[:Tb]
        pr = fc.instance_params(Rall, Tb, L)
        w = fc.compute_weights(max(Tb, 2), c["S_min_m3"], pr["u_max"])
        cfg = fcfg.FalconConfig(T=Tb, L=L, balance="soft", onehot="xy_mixer")
        Q, const, vi, meta = fq.build_qubo(cfg, R, dS, S0, S_min=c["S_min_m3"],
                                           delta_u=pr["delta_u"], levels=pr["levels"], weights=w,
                                           B=c["eta"] * float(R.sum()))
        s, p, _ = qa._ising_terms(Q)
        groups = [[vi.idx(t, l) for l in range(L)] for t in range(Tb)]
        qc = hw._measured_circuit_xy(meta["n_qubits"], s, p, [0.5], [0.5], groups)
        tqc = transpile(qc, fake, optimization_level=3)
        ns.append(meta["n_qubits"]); depths.append(tqc.depth()); twoq.append(tqc.num_nonlocal_gates())

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(ns, depths, "o-", color="tab:blue", label="nuestra QAOA-XY chunked (medido)")
    ax.scatter([28], [1900], color="tab:red", s=70, marker="s",
               label="Ivan (T26/L5 domain-wall, ~1900)")
    for x, d, g in zip(ns, depths, twoq):
        ax.annotate(f"{d}\n({g} 2q)", (x, d), fontsize=7, ha="center", va="bottom")
    ax.set_xlabel("Qubits del bloque")
    ax.set_ylabel("Profundidad transpilada (FakeBrisbane, opt=3)")
    ax.set_title("Recurso en hardware: profundidad transpilada por bloque\n"
                 "bloques chicos = circuitos poco profundos; el costo del XY-mixer crece con el tamaño")
    ax.legend(fontsize=9)
    ax.set_ylim(0, 2100)
    _save(fig, "figC_depth_vs_qubits.png")


# ---------------------------------------------------------------- Fig 4: feasibilidad (datos reales)
def fig_feasibility():
    rows = list(csv.DictReader(open(os.path.join(ROOT, "results", "runs_summary.csv"))))
    insts = ["debug", "small", "medium", "large"]
    # método -> (label, selector de variant preferido)
    methods = [
        ("historical", "histórico u=0", None),
        ("dp", "DP exacto", None),
        ("qaoa", "QAOA sim", None),
        ("qaoa_chunked", "QAOA-chunked sim", None),
        ("qaoa_hardware", "QAOA hardware", None),
        ("qaoa_chunked_hardware", "QAOA-chunked HW", None),
    ]

    def pick(method, inst):
        cands = [r for r in rows if r["method"] == method and r["instance"] == inst
                 and r["SRS"] not in ("", None)]
        if not cands:
            return None
        # SOLO corridas canónicas (ventana base start=0, sin E1 first/middle/stress)
        canon = [r for r in cands if not r.get("window_label")]
        use = canon if canon else cands
        # entre canónicas, la mejor SRS (mejor encoding/p para métodos cuánticos)
        best = max(use, key=lambda r: float(r["SRS"]))
        return {"SRS": float(best["SRS"]), "feasible": best["feasible"] == "True"}

    fig, ax = plt.subplots(figsize=(9, 5))
    nM = len(methods)
    width = 0.8 / nM
    x = np.arange(len(insts))
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, nM))
    for j, (m, lbl, _) in enumerate(methods):
        vals, hatches, present = [], [], []
        for inst in insts:
            r = pick(m, inst)
            vals.append(r["SRS"] if r else np.nan)
            hatches.append("" if (r and r["feasible"]) else "xx")
            present.append(r is not None)
        xs = x + (j - nM / 2) * width + width / 2
        for xi, v, h, pr in zip(xs, vals, hatches, present):
            if not pr or np.isnan(v):
                continue
            ax.bar(xi, v, width, color=colors[j], edgecolor="black", hatch=h,
                   label=lbl if xi == xs[[i for i, p in enumerate(present) if p][0]] else None)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x); ax.set_xticklabels([f"{i}\n(T,L)" for i in insts])
    ax.set_xticklabels(["debug T5/L3", "small T12/L3", "medium T26/L5", "large T52/L5"])
    ax.set_ylabel("SRS (más alto = mejor; negativo)")
    ax.set_title("Feasibilidad del enfoque cuántico (datos reales de results/)\n"
                 "hatch ✗ = infactible. QAOA-chunked+relajaciones da soluciones factibles y cerca del DP")
    ax.legend(fontsize=8, ncol=3, loc="lower right")
    _save(fig, "fig4_feasibility.png")


if __name__ == "__main__":
    fig_time()
    fig_memory()
    fig_qubits()
    fig_dp_vs_brute()
    fig_multireservoir()
    fig_depth()
    fig_feasibility()
    print("\nListo:", FIGDIR)
