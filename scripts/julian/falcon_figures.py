"""Genera las figuras (y sus notas) para la presentacion del Falcon Challenge.

Mapeadas a la rubrica: Problem Formulation (25%), Quantum Impl. (25%), Baseline (20%),
Benchmarking (20%), Presentation (10%). Enfasis: JUSTIFICAR las relajaciones de
restricciones con datos (no afirmaciones).

- Figuras clasicas (DP/baseline/exhaustive) se recomputan EN VIVO (rapido, .venv).
- Resultados cuanticos (qaoa / qaoa_chunked / qaoa_hardware) se leen de
  results/runs_summary.csv (ya corridos, posiblemente en .venv-quantum / IBM HW).

Salidas: results/figures/*.png  +  results/figures/FIGURE_NOTES.md
Uso: .venv/bin/python scripts/julian/falcon_figures.py
"""
from __future__ import annotations

import math
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import falcon_baselines as bl
import falcon_config as fcfg
import falcon_constants as fc
import falcon_data as fd
import falcon_qubo as fq
import falcon_solvers as sv
import falcon_srs as srs
import falcon_storage as st

REPO = Path(__file__).resolve().parent.parent.parent
FIGDIR = REPO / "results" / "figures"
CSV = REPO / "results" / "runs_summary.csv"
INSTANCES = [(5, 3, "debug"), (12, 3, "small"), (26, 5, "medium"), (52, 5, "large"), (52, 7, "large7")]
NOTES: list[dict] = []


def note(fid, title, proves, rubric, talking, procon):
    NOTES.append({"id": fid, "title": title, "proves": proves, "rubric": rubric,
                  "talking": talking, "procon": procon})


def save(fig, fid):
    FIGDIR.mkdir(parents=True, exist_ok=True)
    p = FIGDIR / f"{fid}.png"
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {p.relative_to(REPO)}")


# --------------------------------------------------------------------------- #
def _load():
    df = pd.read_csv(CSV)
    wk = fd.build_weekly_benchmark(write=False)
    c = fc.load_official_constants()
    return df, wk, c


def latest(df, instance, method, variant=None, window=False, L=None):
    """Ultima fila (por timestamp) para (instance, method[, variant][, L]); no-ventana por defecto.

    `L` desambigua large (T52/L5 vs T52/L7), que comparten instance='large'.
    """
    m = (df["instance"] == instance) & (df["method"] == method)
    if variant is not None:
        m &= (df["variant"].astype(str) == variant)
    if L is not None:
        m &= (df["L"] == L)
    if not window:
        m &= ~(df["window_label"].astype(str).str.len() > 0)
    sub = df[m].sort_values("timestamp")
    return sub.iloc[-1] if len(sub) else None


def pick(*rows):
    """Primer row no-None (evita `Series or Series`, que rompe por truthiness ambigua)."""
    for r in rows:
        if r is not None:
            return r
    return None


def _instance_arrays(wk, c, T, L):
    Rall = wk["R_obs_m3_week"].to_numpy()
    R, dS = Rall[:T], wk["DeltaS_obs_m3"].to_numpy()[:T]
    S0 = wk.attrs["S0_m3"]
    p = fc.instance_params(Rall, T, L)
    w = fc.compute_weights(T, c["S_min_m3"], p["u_max"])
    B = c["eta"] * float(R.sum())
    return R, dS, S0, p, w, B


# ---------- A. Justificacion de relajaciones -------------------------------- #
def fig_A2_storage_band(wk, c):
    """S(t) observado vs cotas 0 y S_max: nunca se acercan -> drop es lossless."""
    S = np.concatenate([[wk.attrs["S0_m3"]], wk["S_obs_m3"].to_numpy()]) / 1e6
    Smax, Smin = c["S_max_m3"] / 1e6, c["S_min_m3"] / 1e6
    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.axhspan(0, Smax, color="0.93", label="feasible band [0, S_max]")
    ax.plot(range(len(S)), S, "-o", ms=3, color="#1f77b4", label="observed storage S(t)")
    ax.axhline(Smax, color="#d62728", ls="--", lw=1.2, label=f"S_max = {Smax:,.0f} Mm³")
    ax.axhline(Smin, color="#ff7f0e", ls="--", lw=1.2, label=f"S_min (25%) = {Smin:,.0f} Mm³")
    ax.axhline(0, color="k", lw=1)
    ax.set(xlabel="week", ylabel="storage (Mm³)",
           title="Storage never approaches 0 or S_max  →  bounds 0≤S≤S_max are droppable")
    ax.legend(fontsize=7, loc="upper left")
    save(fig, "A2_storage_band")
    note("A2", "Storage vs physical bounds",
         f"Storage stays ~10-20% of S_max all year; min≈{S.min():.0f}Mm³ vs 0, "
         f"max≈{S.max():.0f}Mm³ vs {Smax:.0f}Mm³.",
         "Quantum Impl. (25%) — justifies storage_bounds='drop' (0 qubits)",
         "Dropping 0≤S≤S_max costs zero optimality here.",
         "Pro: saves O(T) slacks. Con: must re-enable if a wetter regime approaches the bounds.")


def fig_A1_A3_constraints(wk, c):
    """Tabla-figura 'que restriccion ata' + utilizacion de balance |Σu|/B (DP)."""
    rows, utils, labels = [], [], []
    for T, L, name in INSTANCES:
        R, dS, S0, p, w, B = _instance_arrays(wk, c, T, L)
        dp = bl.dp_optimal(R, dS, S0, S_min=c["S_min_m3"], S_max=c["S_max_m3"],
                           delta_u=p["delta_u"], L=L, weights=w, B=B)
        u = dp["u_star"]; S = st.simulate_storage(S0, dS, u)
        util = abs(u.sum()) / B if B else 0.0
        utils.append(util); labels.append(f"{name}\nT{T}/L{L}")
        n_forbidden = int(sum(1 for t in range(T) for a in p["levels"] if R[t] + a < -1e-9))
        rows.append((name, S.min() / 1e6, S.max() / 1e6, util, n_forbidden))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.3), gridspec_kw={"width_ratios": [1.3, 1]})
    # left: balance utilization
    ax1.bar(labels, [u * 100 for u in utils], color="#2ca02c")
    ax1.axhline(100, color="#d62728", ls="--", lw=1, label="budget B (η·ΣR_obs)")
    ax1.set(ylabel="|Σu| / B  (%)", title="Balance constraint is (near-)binding → keep it EXACT",
            ylim=(0, 115))
    for i, u in enumerate(utils):
        ax1.text(i, u * 100 + 2, f"{u*100:.0f}%", ha="center", fontsize=8)
    ax1.legend(fontsize=8)
    # right: decision table
    ax2.axis("off")
    tbl = [["constraint", "binds?", "relaxation", "qubits"],
           ["0 ≤ S ≤ S_max", "never", "drop", "0"],
           ["|Σu| ≤ B", "~100%*", "exact 1-slack", "~log₂(0.8T)"],
           ["R(t) ≥ 0", "sparse", "prohibit lvl", "0"],
           ["one-hot Σ=1", "always", "penalty / XY", "0"]]
    t = ax2.table(cellText=tbl, loc="center", cellLoc="center",
                  colWidths=[0.30, 0.16, 0.30, 0.24])
    t.auto_set_font_size(False); t.set_fontsize(8.5); t.scale(1, 1.6)
    for j in range(4):
        t[(0, j)].set_facecolor("#dddddd"); t[(0, j)].set_text_props(weight="bold")
    ax2.set_title("Per-constraint relaxation decision", fontsize=10)
    save(fig, "A1_A3_constraints")
    util_act = [u for (T, L, name), u in zip(INSTANCES, utils) if not (name in ("debug", "small"))]
    note("A1/A3", "Which constraints bind + how relaxed",
         "Balance binds at %.0f–%.0f%% where the optimizer acts (medium/large); debug/small are u*=0 "
         "(drought) so balance is slack there. Storage bounds never bind; R≥0 forbids only %d–%d levels. "
         "(*table '~100%%' = the binding regime.)"
         % (min(util_act) * 100, max(util_act) * 100,
            min(r[4] for r in rows), max(r[4] for r in rows)),
         "Quantum Impl. (25%) + Benchmarking (20%) — the relaxation rationale",
         "Each relaxation is chosen by whether the constraint actually binds in our data (measured, not assumed).",
         "Pro: minimal qubits, provably lossless for the dropped ones. Con: data-regime-specific.")


def fig_A4_encoding_qubits():
    """Qubits por encoding (one-hot / binary / domain-wall) vs instancia, con techo statevector."""
    names, oh, bina, dw = [], [], [], []
    for T, L, name in INSTANCES:
        names.append(f"{name}\nT{T}/L{L}")
        oh.append(T * L)
        bina.append(T * math.ceil(math.log2(L)))
        dw.append(T * (L - 1))                       # domain-wall: L-1 bits/semana
    x = np.arange(len(names)); wdt = 0.27
    fig, ax = plt.subplots(figsize=(9, 4.3))
    ax.bar(x - wdt, oh, wdt, label="one-hot (T·L)", color="#1f77b4")
    ax.bar(x, dw, wdt, label="domain-wall (T·(L−1)) — used on IBM HW", color="#2ca02c")
    ax.bar(x + wdt, bina, wdt, label="binary (T·⌈log₂L⌉)", color="#ff7f0e")
    ax.axhline(30, color="#d62728", ls="--", lw=1.3, label="statevector / HW block ceiling (~30 qubits)")
    ax.set(xticks=x, ylabel="qubits (full instance)",
           title="Encoding sizing vs ~30-qubit limit → compact encodings + chunking are necessary")
    ax.set_xticklabels(names, fontsize=8)
    ax.legend(fontsize=7.5)
    save(fig, "A4_encoding_qubits")
    note("A4", "Qubit count by encoding (one-hot / domain-wall / binary)",
         "medium T26/L5: one-hot=130q, domain-wall=104q, binary=78q — all >30 → chunking. Ivan's IBM "
         "hardware run used domain-wall (4 bits/week) in blocks 7+7+7+5 (≤30q).",
         "Quantum Impl. (25%) — justifies encoding + chunking choices",
         "Compact encodings + per-block chunking are what make QAOA fit a real device / statevector.",
         "Pro: domain-wall enabled the actual IBM-hardware run; binary fits small in statevector. "
         "Con: binary exact only for L≤4; every full instance still needs chunking.")


def fig_A5_penalty_sweep(wk, c):
    """SRS del optimo del QUBO (exhaustivo) vs multiplicador de penalties: plano => maxabs lo absorbe."""
    T, L = 12, 3
    R, dS, S0, p, w, B = _instance_arrays(wk, c, T, L)
    half = (L - 1) // 2
    mults = [0.1, 0.3, 1.0, 3.0, 10.0, 100.0]
    srs_vals, feas = [], []
    Js = fq._historical_Jscale(S0, dS, c["S_min_m3"], w, T)
    for m in mults:
        pen = {"P_onehot": 10 * Js * m, "P_R": 10 * Js * m, "P_bal_slack": 10 * Js * m}
        cfg = fcfg.FalconConfig(T=T, L=L, penalties=pen)
        Q, const, vi, meta = fq.build_qubo(cfg, R, dS, S0, S_min=c["S_min_m3"],
                                           delta_u=p["delta_u"], levels=p["levels"], weights=w, B=B)
        out = sv.exhaustive_qubo(Q, const, vi, meta, p["levels"], R, B, half)
        srs_vals.append(out["SRS_star"])
        dv = sv.decode_and_verify(out["lv"], vi, p["levels"], R, dS, S0, S_min=c["S_min_m3"],
                                  S_max=c["S_max_m3"], u_max=p["u_max"], B=B, weights=w)
        feas.append(dv["feasible"])
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    ax.plot(mults, srs_vals, "-o", color="#9467bd")
    ax.set_xscale("log")
    ax.set(xlabel="penalty multiplier × base", ylabel="QUBO optimum SRS (exhaustive)",
           title="QUBO optimum is flat across penalty magnitude (maxabs normalization absorbs it)")
    for xi, yi, fe in zip(mults, srs_vals, feas):
        ax.annotate("✓" if fe else "✗", (xi, yi), fontsize=9, ha="center", va="bottom")
    save(fig, "A5_penalty_sweep")
    note("A5", "Penalty magnitude is not the lever",
         "Exhaustive-QUBO SRS is constant (=%.5f) across 0.1–100× penalties (feasible)." % srs_vals[-1],
         "Quantum Impl. (25%) — tuning insight",
         "Penalty multiplier doesn't change the optimum; depth and the XY-mixer are the real levers.",
         "Pro: no brittle penalty tuning. Con: penalties still must be 'large enough' (very small breaks feasibility).")


# ---------- B. Validez de datos --------------------------------------------- #
def fig_B1_deltaS(wk, c):
    """ΔS derivado vs −(oficial) semanal: coinciden exacto (validacion de datos + signo)."""
    if "DeltaS_official_m3" not in wk.columns:
        print("  [B1] sin DeltaS_official_m3, salteo"); return
    der = wk["DeltaS_obs_m3"].to_numpy() / 1e6
    off = wk["DeltaS_official_m3"].to_numpy() / 1e6
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    ax.plot(der, "-", lw=2, color="#1f77b4", label="derived  S(t+1)−S(t)")
    ax.plot(off, "--", lw=1.3, color="#ff7f0e", label="−(official Change-in-Storage), weekly")
    corr = float(np.corrcoef(der, off)[0, 1])
    ax.set(xlabel="week", ylabel="ΔS (Mm³/week)",
           title=f"Derived ΔS == −(official IBWC series), weekly  (corr={corr:+.3f})")
    ax.legend(fontsize=8)
    save(fig, "B1_deltaS_validation")
    note("B1", "ΔS cross-validated vs official IBWC",
         "Weekly derived == −official: corr=%+.3f, err≈0; official used raw drives storage negative." % corr,
         "Baseline (20%) + Problem Formulation (25%) — data validity",
         "Our ΔS matches the official series once the sign convention is fixed; results aren't preliminary.",
         "Pro: data is validated + sign trap documented. Con: official daily series has gaps/opposite sign.")


def fig_B2_drought(wk, c):
    """Storage como % de capacidad (1 año) con S_min(25%): regimen de sequia."""
    S = np.concatenate([[wk.attrs["S0_m3"]], wk["S_obs_m3"].to_numpy()]) / c["S_max_m3"] * 100
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    ax.fill_between(range(len(S)), 0, S, color="#8c564b", alpha=0.35)
    ax.plot(range(len(S)), S, "-o", ms=3, color="#8c564b", label="storage (% capacity)")
    ax.axhline(25, color="#ff7f0e", ls="--", lw=1.5, label="S_min = 25% (conservation threshold)")
    ax.set(xlabel="week", ylabel="% of S_max", ylim=(0, 30),
           title="Multi-year drought: storage ~10–20% of capacity — 0/2182 days above S_min (6-yr audit)")
    ax.legend(fontsize=8)
    save(fig, "B2_drought_regime")
    note("B2", "Drought regime frames the problem",
         "Storage 10–20% of S_max; 0 of 2182 days (2020–2026) above S_min → C_crit dominates the SRS.",
         "Problem Formulation (25%) — societal/hydrological context",
         "The reservoir is critically low all year, so avoiding deeper shortfall (C_crit) is the priority.",
         "Pro: explains why u=0 optima appear at coarse L. Con: window is a drought; wetter years may differ.")


# ---------- C. Baseline & correctness --------------------------------------- #
def fig_C1_srs_bars(df, wk, c):
    """SRS por instancia: baselines + DP + mejor cuantico, con factibilidad anotada."""
    fig, axes = plt.subplots(1, len(INSTANCES), figsize=(15, 4.2), sharey=False)
    for ax, (T, L, name) in zip(axes, INSTANCES):
        inst = "large" if name.startswith("large") else name
        bars = []
        h = latest(df, inst, "historical", L=L); bars.append(("hist", h))
        tp = latest(df, inst, "threshold", "pure", L=L); bars.append(("thr-pure", tp))
        tb = latest(df, inst, "threshold", "balanced", L=L); bars.append(("thr-bal", tb))
        d = latest(df, inst, "dp", L=L); bars.append(("DP", d))
        q = pick(latest(df, inst, "qaoa", L=L), latest(df, inst, "qaoa_chunked", L=L),
                 latest(df, inst, "qubo_exhaustive", L=L))
        if q is not None:
            bars.append(("quantum", q))
        bars = [(lab, r) for lab, r in bars if r is not None]
        if not bars:
            ax.set_visible(False); continue
        labs = [lab for lab, _ in bars]
        vals = [float(r["SRS"]) for _, r in bars]
        feas = [bool(r["feasible"]) for _, r in bars]
        # feasible = verde solido; infeasible = rojo con hachurado (no comparar su altura)
        best_feas = max((v for v, fe in zip(vals, feas) if fe), default=None)
        cont = ax.bar(labs, vals, color=["#2ca02c" if fe else "#d62728" for fe in feas])
        for bar, fe, v in zip(cont, feas, vals):
            if not fe:
                bar.set_hatch("////"); bar.set_edgecolor("white")
                ax.text(bar.get_x() + bar.get_width() / 2, v, "INFEASIBLE\n(over budget)",
                        ha="center", va="top", fontsize=6, color="#7a0000")
            elif best_feas is not None and abs(v - best_feas) < 1e-9:
                ax.text(bar.get_x() + bar.get_width() / 2, v, "★ best\nfeasible",
                        ha="center", va="bottom", fontsize=6.5, color="#0a5a0a")
        ax.set_title(f"{name} T{T}/L{L}", fontsize=9)
        ax.tick_params(axis="x", labelsize=7, rotation=35)
        ax.axhline(0, color="k", lw=0.6)
    axes[0].set_ylabel("SRS (higher = better)")
    fig.suptitle("SRS by method — compare heights ONLY among feasible (solid green) bars. "
                 "Hatched red = INFEASIBLE: 'better' SRS only by overspending the balance budget |Σu|≤B.",
                 y=1.03, fontsize=9)
    save(fig, "C1_srs_by_method")
    note("C1", "SRS comparison per instance (feasibility-aware)",
         "Among FEASIBLE methods, DP ≥ historical ≥ threshold-balanced; threshold-pure/clamped look "
         "higher (e.g. medium −0.289) but are INFEASIBLE (|Σu|−B>0). Feasible-but-worse-than-hist "
         "(debug/small threshold-balanced) is legitimate — u=0 is optimal in the drought.",
         "Baseline (20%) + Benchmarking (20%)",
         "Compare only feasible bars; exact DP is the honest strong baseline for the quantum comparison.",
         "Pro: valid metric with feasibility made explicit. Con: naive threshold 'wins' only if you "
         "ignore the balance constraint it violates.")


# ---------- D. Quantum & benchmarking --------------------------------------- #
def fig_D1_quantum_vs_baseline(df, wk, c):
    """Medium T26/L5: full DP >= DP-chunked >= QAOA-chunked (descomposicion del gap)."""
    T, L = 26, 5
    R, dS, S0, p, w, B = _instance_arrays(wk, c, T, L)
    h = latest(df, "medium", "historical")
    dp_full = bl.dp_optimal(R, dS, S0, S_min=c["S_min_m3"], S_max=c["S_max_m3"],
                            delta_u=p["delta_u"], L=L, weights=w, B=B)["SRS_star"]
    import falcon_chunked as ch
    dpc = ch.staged_solve(wk["R_obs_m3_week"].to_numpy(), wk["DeltaS_obs_m3"].to_numpy(), S0, T,
                          S_min=c["S_min_m3"], S_max=c["S_max_m3"], L=L, block_size=5, solver="dp")["SRS"]
    # QAOA-chunked: preferir el factible (variant *_xy); sino el que haya
    qc = pick(latest(df, "medium", "qaoa_chunked", "blk5_p1_xy"),
              latest(df, "medium", "qaoa_chunked"))
    labels = ["historical\n(u=0)", "QAOA-chunked", "DP-chunked", "full DP\n(optimum)"]
    vals = [float(h["SRS"]), float(qc["SRS"]) if qc is not None else np.nan, dpc, dp_full]
    feas = [True, bool(qc["feasible"]) if qc is not None else False, True, True]
    cols = ["#7f7f7f", "#2ca02c" if feas[1] else "#d62728", "#1f77b4", "#9467bd"]
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    ax.bar(labels, vals, color=cols)
    ax.axhline(dp_full, color="#9467bd", ls=":", lw=1)
    ax.set(ylabel="SRS", title="Medium T26/L5: quantum vs baseline + gap decomposition")
    for i, (v, fe) in enumerate(zip(vals, feas)):
        if not np.isnan(v):
            ax.text(i, v, f"{v:.4f}\n{'feasible' if fe else 'INFEASIBLE'}", ha="center",
                    va="bottom" if v > -0.35 else "top", fontsize=8)
    qc_txt = f"{float(qc['SRS']):.4f} ({'feasible' if feas[1] else 'infeasible'})" if qc is not None else "pending"
    save(fig, "D1_quantum_vs_baseline_medium")
    note("D1", "Quantum vs baseline on the official benchmark",
         "full DP=%.4f ≥ DP-chunked=%.4f ≥ QAOA-chunked=%s; chunking gap vs QAOA gap separated."
         % (dp_full, dpc, qc_txt),
         "Benchmarking (20%) + Quantum Impl. (25%) — the headline comparison",
         "QAOA solves each block feasibly; the gap to full-DP is the price of chunking (a scaling limit), per spec §7.",
         "Pro: honest, feasible quantum result at 25q/block. Con: eta_local chunking caps ΔSRS; no quantum advantage (expected).")


def fig_D2_approx_ratio(df):
    """AR por instancia/encoding, incluyendo hardware IBM si existe."""
    pts = []
    for inst, method, variant, lab in [
        ("debug", "qaoa", "onehot_p1", "debug one-hot (sim)"),
        ("debug", "qaoa_hardware", "onehot_p1", "debug one-hot (IBM HW)"),
        ("small", "qaoa", "binary_p1", "small binary (sim)"),
        ("medium", "qaoa_chunked", "blk5_p1_xy", "medium XY-chunked (sim)")]:
        r = latest(df, inst, method, variant)
        if r is not None and pd.notna(r["approximation_ratio"]):
            pts.append((lab, float(r["approximation_ratio"]), bool(r["feasible"])))
    if not pts:
        print("  [D2] sin datos AR, salteo"); return
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    cols = ["#2ca02c" if fe else "#d62728" for _, _, fe in pts]
    ax.bar([l for l, _, _ in pts], [a for _, a, _ in pts], color=cols)
    ax.axhline(1.0, color="k", ls="--", lw=1, label="AR=1.0 (optimum)")
    ax.set(ylabel="approximation ratio (SRS/DP*)", title="QAOA approximation ratio by instance/encoding")
    ax.tick_params(axis="x", labelsize=8, rotation=15)
    ax.legend(fontsize=8)
    save(fig, "D2_approx_ratio")
    note("D2", "QAOA quality + XY-mixer + hardware",
         "debug AR=1.0 on sim AND IBM hardware; XY-mixer keeps 100% probability in the one-hot subspace.",
         "Quantum Impl. (25%)",
         "QAOA reaches the optimum on debug (incl. real hardware); XY-mixer removes the one-hot penalty cleanly.",
         "Pro: validated on hardware; XY-mixer principled. Con: AR degrades at p=1 as size grows.")


def fig_D4_windows(df):
    """Robustez de ventana: DP vs historical por ventana (first/middle/stress), con ΔSRS.

    ΔSRS = SRS_DP − SRS_hist por ventana: 0 => el optimo es u=0 (sin headroom, sequia);
    >0 => el DP encuentra una politica no trivial. Muestra si el optimo es no trivial en
    cada ventana, no solo el valor de SRS.
    """
    def _srs(inst, method, wl):
        sub = df[(df["instance"] == inst) & (df["method"] == method) & (df["window_label"] == wl)]
        return float(sub["SRS"].iloc[-1]) if len(sub) else np.nan

    order = ["first", "middle", "stress"]
    insts = ["debug", "small", "medium"]
    fig, axes = plt.subplots(1, len(insts), figsize=(13, 4.4), sharey=False)
    for ax, inst in zip(axes, insts):
        x = np.arange(len(order)); wdt = 0.38
        hist = [_srs(inst, "historical", wl) for wl in order]
        dp = [_srs(inst, "dp", wl) for wl in order]
        ax.bar(x - wdt / 2, hist, wdt, label="historical (u=0)", color="#7f7f7f")
        ax.bar(x + wdt / 2, dp, wdt, label="DP optimum", color="#1f77b4")
        for xi, (hv, dv) in enumerate(zip(hist, dp)):
            if not (np.isnan(hv) or np.isnan(dv)):
                d = dv - hv
                ax.annotate(f"ΔSRS\n{d:+.4f}", (xi, max(hv, dv)), ha="center", va="bottom",
                            fontsize=7, color=("#0a5a0a" if d > 1e-9 else "#555555"))
        ax.set_title(f"{inst}", fontsize=10)
        ax.set_xticks(x); ax.set_xticklabels(order, fontsize=8)
        if ax is axes[0]:
            ax.set_ylabel("SRS (higher = better)")
        ax.legend(fontsize=7)
    fig.suptitle("Window robustness: DP vs historical per window (ΔSRS>0 ⇒ non-trivial optimum found). "
                 "first=weeks[0,T)  ·  middle=start (52−T)//2  ·  stress=lowest-mean-storage window",
                 y=1.03, fontsize=9)
    save(fig, "D4_window_robustness")
    note("D4", "Window robustness (E1) — DP vs historical per window",
         "ΔSRS(DP−hist): debug/small ≈ 0 in ALL windows (u=0 optimal, drought+coarse L); medium >0 in "
         "every window → the optimum is genuinely non-trivial there, not a start-of-year artifact.",
         "Benchmarking (20%)",
         "Windows: first=weeks[0,T); middle=start (52−T)//2; stress=T-week window of lowest mean storage "
         "(deepest deficit). Comparing DP to historical per window shows WHERE optimization actually helps.",
         "Pro: separates 'found optimum' from 'optimum beats baseline' across windows. Con: all windows "
         "are within one drought year.")


# ---------- E. Scaling ------------------------------------------------------- #
# Estados DP explorados (medido, snapshot 2026-06-30, docs/ANALISIS_DP_Y_RESULTADOS.md §3).
DP_STATES = {("debug", 5, 3): 100, ("small", 12, 3): 433, ("medium", 26, 5): 6505,
             ("large", 52, 5): 25086, ("large7", 52, 7): 49739}


def fig_E1_scaling(df):
    """Colapso exponencial->polinomial: espacio L^T vs estados DP realmente explorados."""
    names, spaces, states, rts = [], [], [], []
    for T, L, name in INSTANCES:
        inst = "large" if name.startswith("large") else name
        r = latest(df, inst, "dp", L=L)
        if r is None:
            continue
        st_cnt = DP_STATES.get((name, T, L))
        if st_cnt is None:
            continue
        names.append(f"{name}\nT{T}/L{L}")
        spaces.append(L ** T); states.append(st_cnt); rts.append(float(r["runtime_seconds"]))
    x = np.arange(len(names))
    fig, ax = plt.subplots(figsize=(9, 4.6))
    ax.plot(x, [math.log10(s) for s in spaces], "-o", color="#d62728", lw=2,
            label="log₁₀(L^T) candidate schedules (brute)")
    ax.plot(x, [math.log10(s) for s in states], "-s", color="#1f77b4", lw=2,
            label="log₁₀(DP states actually explored)")
    for xi, sp, stt in zip(x, spaces, states):
        ax.annotate(f"{sp:.0e}", (xi, math.log10(sp)), fontsize=7, color="#d62728",
                    ha="center", va="bottom")
        ax.annotate(f"{stt:,}", (xi, math.log10(stt)), fontsize=7, color="#1f77b4",
                    ha="center", va="top")
    ax.set(xticks=x, ylabel="log₁₀(count)",
           title="DP collapses the search: L^T→10⁴³ candidates, but only ~10⁴ states explored (O(T²L²))")
    ax.set_xticklabels(names, fontsize=8)
    ax.legend(fontsize=8, loc="center left")
    ax2 = ax.twinx()
    ax2.plot(x, rts, ":^", color="#2ca02c", label="DP runtime (s)")
    ax2.set_ylabel("DP runtime (s)", color="#2ca02c")
    ax2.tick_params(axis="y", labelcolor="#2ca02c")
    save(fig, "E1_scaling")
    note("E1", "Scaling: L^T search space vs DP states explored (the collapse)",
         "medium L^T=1.49e18 → only 6,505 DP states (0.02s); large L7 8.81e43 → 49,739 states (0.20s). "
         "DP state=(t, C_t=Σk_j, k_prev): storage depends only on the integer cumulative sum C_t, so "
         "exponentially many schedules fuse into O(T²·L²) states — exact, lossless (brute==dp).",
         "Benchmarking (20%) — scaling analysis",
         "Exploiting optimal substructure turns a 10⁴³ search into ~10⁴ states: exact ground truth at every size.",
         "Pro: DP is exact + sub-second at all scales. Con: it's a classical baseline; QAOA can't match its scaling.")


def fig_E2_dsrs(df):
    """ΔSRS del DP vs instancia: el headroom crece con T/L."""
    names, ds = [], []
    for T, L, name in INSTANCES:
        inst = "large" if name.startswith("large") else name
        r = latest(df, inst, "dp", L=L)
        if r is None or pd.isna(r["dSRS_vs_historical"]):
            continue
        names.append(f"{name}\nT{T}/L{L}"); ds.append(float(r["dSRS_vs_historical"]))
    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.bar(names, ds, color="#2ca02c")
    ax.set(ylabel="ΔSRS vs historical (DP)",
           title="Optimization headroom grows with horizon T and levels L")
    for i, v in enumerate(ds):
        ax.text(i, v, f"{v:+.4f}", ha="center", va="bottom", fontsize=8)
    ax.tick_params(axis="x", labelsize=8)
    save(fig, "E2_dsrs_by_instance")
    note("E2", "ΔSRS headroom by instance",
         "ΔSRS_vs_hist: debug/small≈0 (drought+coarse L), medium=+0.021, large=+0.051/+0.075.",
         "Benchmarking (20%) + Problem Formulation (25%)",
         "The debug instances are trivially u=0; the real optimization payoff is at L=5, longer horizons.",
         "Pro: shows where the method matters. Con: small ΔSRS overall (drought caps achievable gains).")


# ---------- F. Policy -------------------------------------------------------- #
def fig_F1_policy(wk, c):
    """Politica optima DP en medium: u(t) + trayectoria de storage."""
    T, L = 26, 5
    R, dS, S0, p, w, B = _instance_arrays(wk, c, T, L)
    dp = bl.dp_optimal(R, dS, S0, S_min=c["S_min_m3"], S_max=c["S_max_m3"],
                       delta_u=p["delta_u"], L=L, weights=w, B=B)
    u = dp["u_star"] / 1e6
    S = st.simulate_storage(S0, dS, dp["u_star"]) / 1e6
    Sh = st.simulate_storage(S0, dS, np.zeros(T)) / 1e6
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
    ax1.bar(range(1, T + 1), u, color="#1f77b4")
    ax1.set(ylabel="u(t) (Mm³/week)", title=f"DP optimal release-adjustment policy (medium T26/L5, ΔSRS=+{dp['SRS_star']-srs.compute_srs(srs.compute_costs(Sh*1e6,np.zeros(T),c['S_min_m3']),w):.4f})")
    ax1.axhline(0, color="k", lw=0.6)
    ax2.plot(range(T + 1), Sh, "--", color="#7f7f7f", label="historical (u=0)")
    ax2.plot(range(T + 1), S, "-o", ms=3, color="#2ca02c", label="optimized S_opt")
    ax2.axhline(c["S_min_m3"] / 1e6, color="#ff7f0e", ls="--", lw=1, label="S_min")
    ax2.set(xlabel="week", ylabel="storage (Mm³)"); ax2.legend(fontsize=8)
    save(fig, "F1_dp_policy_medium")
    note("F1", "Optimal policy (medium)",
         "DP conserves in ~11/26 weeks (u=−2Δu), lifting storage while respecting the balance budget.",
         "Problem Formulation (25%) + Presentation (10%)",
         "The optimizer redistributes releases in time (not just withholding) to reduce critical shortfall.",
         "Pro: interpretable, operationally meaningful. Con: gains are modest given the drought.")


def write_notes():
    lines = ["# Figure notes — Falcon presentation\n",
             "> Generado por `scripts/julian/falcon_figures.py`. Cada figura: qué prueba, rúbrica, "
             "talking point y pro/con. PNGs en `results/figures/`.\n",
             "## Feasibility explainer (para C1 / D1 / D4)\n",
             "Una política `u(t)` es **factible** solo si cumple las **4 restricciones oficiales**: "
             "`R(t)=R_obs+u≥0`, `|u(t)|≤u_max`, `0≤S(t)≤S_max`, y **`|Σu(t)|≤B=η·ΣR_obs`** (presupuesto "
             "de balance). En nuestros datos las primeras tres casi nunca atan: **la que decide "
             "factibilidad es el balance**. *Factibilidad y SRS son ejes independientes*: una política "
             "puede ser factible y peor que el histórico (p.ej. threshold-balanced en debug/small — "
             "legítimo, en sequía `u=0` ya es óptimo), y una **infactible puede tener un SRS más alto "
             "solo por gastar de más el presupuesto** (p.ej. threshold-pure en medium, −0.289, pero "
             "`|Σu|−B>0`). **Regla de lectura: comparar SRS solo entre barras factibles.**\n",
             "## Definición de ventanas (para D4)\n",
             "*first* = semanas `[0,T)`; *middle* = start `(52−T)//2`; *stress* = la ventana de T "
             "semanas con **menor storage medio** (déficit más profundo, la más exigente).\n"]
    for n in NOTES:
        lines += [f"## {n['id']} — {n['title']}  (`{n['id'].replace('/','_')}...png`)",
                  f"- **Proves:** {n['proves']}",
                  f"- **Rubric:** {n['rubric']}",
                  f"- **Talking point:** {n['talking']}",
                  f"- **Pro/con:** {n['procon']}\n"]
    (FIGDIR / "FIGURE_NOTES.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"  wrote {(FIGDIR/'FIGURE_NOTES.md').relative_to(REPO)}")


def main():
    df, wk, c = _load()
    print("Generando figuras en results/figures/ ...")
    for f, args in [(fig_A2_storage_band, (wk, c)), (fig_A1_A3_constraints, (wk, c)),
                    (fig_A4_encoding_qubits, ()), (fig_A5_penalty_sweep, (wk, c)),
                    (fig_B1_deltaS, (wk, c)), (fig_B2_drought, (wk, c)),
                    (fig_C1_srs_bars, (df, wk, c)), (fig_D1_quantum_vs_baseline, (df, wk, c)),
                    (fig_D2_approx_ratio, (df,)), (fig_D4_windows, (df,)),
                    (fig_E1_scaling, (df,)), (fig_E2_dsrs, (df,)), (fig_F1_policy, (wk, c))]:
        try:
            f(*args)
        except Exception as e:
            print(f"  [WARN] {f.__name__} fallo: {type(e).__name__}: {e}")
    write_notes()
    print("Listo.")


if __name__ == "__main__":
    main()
