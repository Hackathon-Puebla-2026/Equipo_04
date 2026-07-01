"""Smoke test de Fase 0: datos -> constantes oficiales -> SRS de baselines.

Demuestra que con las constantes oficiales (S_min ~ 822M m^3) el problema es no
trivial (SRS != 0) y que el DP exacto domina a los baselines.

Uso: python scripts/falcon_smoke_fase0.py
(Al ejecutarse desde scripts/, los modulos hermanos se importan directo.)
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # permite import desde cualquier CWD

import falcon_baselines as bl
import falcon_constants as fc
import falcon_data as fd
import falcon_srs as srs
import falcon_storage as st


def _srs_of(u, S0, dS, S_min, weights):
    S = st.simulate_storage(S0, dS, u)
    return srs.compute_srs(srs.compute_costs(S, u, S_min), weights), S


def main(T: int = 12, L: int = 5) -> None:
    df = fd.build_weekly_benchmark(write=False)
    c = fc.load_official_constants()
    R = df["R_obs_m3_week"].to_numpy()
    dS = df["DeltaS_obs_m3"].to_numpy()
    S0 = df.attrs["S0_m3"]

    p = fc.instance_params(R, T, L)
    w = fc.compute_weights(T, c["S_min_m3"], p["u_max"])
    B = c["eta"] * float(R[:T].sum())

    print(f"Instancia T={T}, L={L}  (semanas disponibles: {df.attrs['n_weeks']})")
    print(f"S_max = {c['S_max_m3']:,.0f} m^3 | S_min = {c['S_min_m3']:,.0f} m^3 (oficial)")
    print(f"delta_u = {p['delta_u']:,.0f} m^3/sem | u_max = {p['u_max']:,.0f} | B = {B:,.0f}")
    print(f"weights w1={w['w1']:.3e} w2={w['w2']:.3e} w3={w['w3']:.3e}\n")

    uh = bl.historical(T)
    srh, Sh = _srs_of(uh, S0, dS[:T], c["S_min_m3"], w)

    ut = bl.threshold_rule(R[:T], dS[:T], S0, c["S_min_m3"], p["delta_u"])
    srt, _ = _srs_of(ut, S0, dS[:T], c["S_min_m3"], w)

    dp = bl.dp_optimal(R[:T], dS[:T], S0, S_min=c["S_min_m3"], S_max=c["S_max_m3"],
                       delta_u=p["delta_u"], L=L, weights=w, B=B)

    print(f"{'baseline':12s} {'SRS':>16s}  {'dSRS vs hist':>16s}")
    print(f"{'historical':12s} {srh:16.6e} {0.0:16.6e}")
    print(f"{'threshold':12s} {srt:16.6e} {srt - srh:16.6e}")
    print(f"{'dp_optimal':12s} {dp['SRS_star']:16.6e} {dp['SRS_star'] - srh:16.6e}  feasible={dp['feasible']}")
    print(f"\nmin storage historico = {Sh.min():,.0f} m^3  (< S_min => problema no trivial: "
          f"{Sh.min() < c['S_min_m3']})")
    assert srh != 0.0, "SRS historico = 0: revisar constantes (deberia ser no trivial)"
    print("\nSMOKE FASE 0 OK")


if __name__ == "__main__":
    main()
