# Auditoría de datos y fuerza bruta (hallazgos)

> **Tipo: guía / hallazgos** (no es spec). Documenta una auditoría escéptica del pipeline: ¿están bien
> el código de fuerza bruta, las restricciones y los datos? ¿Hay soluciones no triviales (`u≠0`)?
> Script: `scripts/julian/falcon_brute_audit.py`. Datos: `FalconChallenge/data/` + `falcon_data.py`.

## TL;DR

- **Fuerza bruta: CORRECTA.** Dos implementaciones independientes coinciden exactamente en las 41
  ventanas T12/L3, y ambas hallan `u≠0` en casos plantados donde el óptimo no es cero.
- **Restricciones: CORRECTAS.** El filtro de factibilidad del brute coincide con
  `falcon_storage.check_constraints` (0 discrepancias en 20 000 muestras) y con el spec. La **restricción
  que ata es el balance** `|Σu| ≤ B = η·ΣR_obs`.
- **Datos: CONSISTENTES.** Unidades correctas; `S_obs/S_max` coincide **exacto** con la serie
  independiente de `%Conservación`. El `ΔS` oficial coincide con el derivado (signo-corregido).
- **`u=0` en T12/L3 NO es un bug:** es propiedad del **régimen de sequía** (storage siempre `<< S_min`).
  Con **L=5 / T mayor** el óptimo pasa a ser **no trivial**.

## 1. Fuerza bruta (correctitud)

- **Dos implementaciones independientes** coinciden: `brute_vectorized` (numpy, contador base-L
  vectorizado, en el script de auditoría) vs `falcon_baselines.brute_force_optimal` (itertools). En las
  **41 ventanas T12/L3**: `max|ΔSRS| = 0.00e+00`, ambas enumeran los `3¹² = 531 441` combos.
- **Caso plantado micro (T=3, a mano):** con `Δu=10`, `S0=100`, `ΔS=[-50,0,0]`, `S_min=100`,
  `w=(1, 1e-3, 1e-3)`, el óptimo es `u*=[-10,-10,-10]` (SRS −2900 vs −7500 de `u=0`). El brute lo
  encuentra. ⇒ el brute **sí** halla `u≠0` cuando conviene.
- **Caso plantado datos reales (T12/L5, objetivo solo-`C_crit`, `S_min` sobre el storage):** óptimo
  `u≠0` en 12/12 semanas. Confirma el pipeline completo (arrays reales, storage, costos, factibilidad,
  enumeración).

## 2. Restricciones

Para una ventana T12/L3 representativa: `Δu=2.745M`, `u_max=Δu`, `B=η·ΣR_obs≈12.99M` (`B/Δu≈4.73`).
Bajo `u=0`: `min S≈349M`, `max S≈402M` → las cotas `0≤S≤S_max` **no atan**; `R≥0` rara vez ata; **ata el
balance** (presupuesto ~4.7·Δu de retiro total). El filtro inline del brute (`R≥0`, `0≤S≤S_max`,
`|Σu|≤B`) da el **mismo veredicto** que `check_constraints` en 20 000 schedules aleatorios (0
discrepancias).

## 3. Datos

- **Unidades:** storage `Value (m^3)` (~3.75e8), release `Value (m^3/s)` (~12). Sin bug de escala
  (TCM/m³). `S0` coincide con `df.attrs['S0_m3']`.
- **Cross-check independiente:** `S_obs/S_max` vs la serie `Percentage.Conservation-Web-Telemetry`
  coincide **exacto** (`max Δ = 0.000 pp`, media 12.6%).
- **Régimen (6 años):** integrando el `ΔS` oficial (2020-2026) el storage queda en **~6-15% de
  capacidad, 0/2182 días sobre `S_min`** (25%). Sequía multianual real. (La reconstrucción a 6 años
  deriva y llega a valores negativos por acumulación de ruido diario: **no** sirve para storage absoluto
  fuera del año con anclaje observado.)

### 3.1 El `ΔS` oficial y la "trampa de signo"

El dataset oficial `Discharge.Total.Change-in-Storage@08461200` (diario, m³):

- 43% de valores diarios negativos — **normal** (storage baja cuando release > inflow).
- ⚠️ **Convención de signo OPUESTA** a la ec. de balance del spec `S(t+1)=S(t)+ΔS−u`. Usado crudo:
  suma anual **−272.0M** mientras el storage observado **subió +273.2M**; `S_opt(u=0)` **deriva −545.7M**
  y cae a ~3% de capacidad. Esta es la "trampa de valores negativos" que motivó derivar `ΔS` del storage.
- ✅ **Negado y agregado por semana coincide EXACTO con el derivado** `S(t+1)−S(t)`:
  `corr=+1.000`, `median|err|=0.000M`, misma suma anual. (El ruido diario de 0.78M se telescopea al
  sumar por semana; además hay 3 huecos multi-día en el diario.)

| resolución | comparación | median\|err\| | corr | suma |
|---|---|---:|---:|---|
| diaria | `+oficial` vs derivado | 2.27M | −0.77 | −272 vs +273 |
| diaria | `−oficial` vs derivado | 0.78M | +0.77 | +272 vs +273 |
| **semanal** | **`−oficial` vs derivado** | **0.000M** | **+1.000** | **+273.2 vs +273.2** |

**Decisión:** el **driver sigue siendo el derivado** (robusto al signo y a los huecos), y el oficial se
usa como **invariante de cross-validación** (`falcon_data.validate_deltaS_vs_official`, columna
`DeltaS_official_m3`). Con esto **`ΔS` deja de ser "preliminar"**: está validado contra la serie oficial
del IBWC. Los números **no cambian** (mismo `ΔS`), así que no se re-corrió nada.

## 4. ¿Cuándo el óptimo es `u≠0`? (búsqueda por fuerza bruta)

- **L=3, T=12 — todas las ventanas:** `0/41` con óptimo `u≠0`. El mejor schedule `u≠0` **pierde** vs
  `u=0` por ≥ 5.09e-3 SRS en toda ventana. (Ojo: schedules **factibles** con `u≠0` existen siempre
  -p.ej. la regla de umbral-; lo que no existe es un **óptimo** con `u≠0`.)
- **L=5, T=12 — first/middle/stress (brute vectorizado, 5¹²≈244M):** óptimo **`u≠0`**:
  first 4/12 (ΔSRS +2.74e-3), middle 4/12 (+2.80e-3), stress 5/12 (+4.75e-3).
- **Sensibilidad (barrido de `S_min` a T12/L3, pesos oficiales):** ningún valor de `S_min` vuelve `u≠0`
  el óptimo a **L=3** — es la **granularidad de nivel** (a L=3 el paso máximo es ±Δu) sumada al régimen
  de sequía lo que lo mantiene trivial. A **L=5** (pasos ±2Δu) el óptimo ya se inclina a `u≠0`.

**Conclusión:** `T12/L3` es una instancia de *debug* cuyo óptimo `u=0` es una propiedad real
(sequía + niveles gruesos), **no** un bug. El `ΔSRS` con headroom vive en **L=5+** (medium `T26/L5`
tiene óptimo no trivial en las 27 ventanas). Enfocar QAOA/benchmarking ahí.

## 5. Reproducir

```bash
.venv/bin/python scripts/falcon_data.py            # imprime la validación ΔS vs oficial (invariante OK)
.venv/bin/python scripts/julian/falcon_brute_audit.py   # auditoría completa (~5-6 min por el L=5)
```
