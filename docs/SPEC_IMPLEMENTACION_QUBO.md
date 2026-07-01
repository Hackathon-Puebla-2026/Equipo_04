# SPEC de implementación: pipeline QUBO Falcon

Plan activo del equipo: un **builder de QUBO configurable, reutilizado por todos los solvers**,
construido por fases (MVP primero) y atacado tarea por tarea con el checklist de abajo. Cada tarea es
suficientemente específica para empezar a codear. Digest de hallazgos: `docs/HALLAZGOS_CLAVE.md`;
auditoría de datos/brute: `docs/AUDITORIA_DATOS_Y_BRUTE.md`.

---

## 0. Objetivo y alcance

Encontrar una política de ajuste de liberaciones `u(t)` para Falcon que maximice el SRS, formulada como
QUBO y resuelta con baselines clásicos, DP exacto y QAOA. Entregar comparación cuántico vs baseline,
escalamiento (T12→T26→T52) y reportes reproducibles. El reservorio modelado es simplificado (no es
política operativa real).

---

## 1. Cómo nos evalúan (rúbrica) y prioridades

`FalconChallenge/HackathonRubric.md`:

| Criterio | Peso | Qué cubre | Tareas de este spec que lo alimentan |
|---|---:|---|---|
| Problem Formulation | 25% | problema claro, literatura, SDG, impacto social | §11 writeup (SDG 6.4/6.5/13.1, refs [1-3] del spec), `docs/FalconChallenge_V6.md` |
| Baseline | 20% | métricas válidas, baseline state-of-the-art | SRS + histórico + umbral + **DP exacto** + baseline informado por literatura (Fase 0) |
| Quantum Implementation | 25% | solución cuántica correcta, código limpio, documentado, enfoque justificado | Fase 3 QAOA + justificación de encoding (one-hot/domain-wall/XY-mixer), `scripts/falcon_*.py` con docstrings |
| Benchmarking | 20% | escalamiento + comparación cuántico vs baseline + resultados | Fase 3 tablas ΔSRS/runtime/escala, QAOA vs DP/baselines |
| Quality of Presentation | 10% | entendimiento, resultados, pros/cons, comunicación | §11 writeup + slides (en inglés) |

Requisitos extra: **explicar el uso de IA** en el proyecto; **slides en inglés**; presentación oral en
inglés o español. Implicación de prioridad: Quantum (25%) y Formulation (25%) pesan más; Baseline y
Benchmarking 20% c/u dependen fuertemente del **DP exacto** (ground truth) y de las tablas de
escalamiento. No dejar el writeup para el final.

---

## 2. Decisiones fijadas (convenciones)

- **Constantes oficiales** de `FalconChallenge/data/falcon_reservoir_constants.json`:
  `S_max = 3,288,726,000 m³`, `S_min = 0.25·S_max ≈ 822,181,500 m³`. (Las preliminares del equipo
  usaban `S_max=máx observado`; se descartan.)
- **Unidades: m³ (SI) en todo el pipeline y reportes.** Conversión de release: integración nativa
  15-min `Σ(caudal·900 s)`.
- **`Δu` por instancia**: `Δu = 0.25 · mediana del release semanal sobre la ventana de T semanas`.
  Además **registrar** el `Δu` de ventana completa (52 semanas) como referencia en la metadata de cada
  corrida. `u_max = 2Δu`, niveles `L=5: {−2,−1,0,1,2}·Δu`, `L=3: {−1,0,1}·Δu`.
- **`S_obs(t)` = valor de cierre de semana** (último valor diario de la semana). `S0` = frontera
  inicial. `ΔS_obs(t) = S(t+1) − S(t)` es el **driver**, ahora **cross-validado** contra el dataset
  oficial `Discharge.Total.Change-in-Storage@08461200` (ya presente): a nivel semanal, derivado ==
  **−(oficial)** EXACTO (corr +1.000, err≈0). ⚠️ El oficial usa **signo opuesto** (usado crudo, el
  storage se va a negativo); por eso el driver sigue siendo el derivado. `ΔS` **deja de ser
  preliminar**. Invariante: `falcon_data.validate_deltaS_vs_official`; detalle:
  `docs/AUDITORIA_DATOS_Y_BRUTE.md`.
- **Semanas: bins fijos de 7 días anclados al primer timestamp; se descarta el remanente parcial →
  52 semanas limpias.** (ivan tiene 53 por usar semanas calendario con bordes parciales; la diferencia
  es esperada y se explica al validar.)
- **Código compartido: módulos planos `scripts/falcon_*.py`** (GUIDELINES §3). Experimentos
  personales en `scripts/<name>/`, `notebooks/<name>/`.
- **Dataset semanal: reconstruido canónicamente** por `falcon_data.py` desde `FalconChallenge/data/`;
  se valida contra `notebooks/ivan/Act2/falcon_weekly_benchmark_preliminary.csv`. Se guarda en
  `data/processed/` (gitignored, reproducible).
- **Abstracción de encoding de primera clase** (one-hot / domain-wall / binary intercambiables).

---

## 3. Unidades y normalización numérica

Todo en m³. Pero en m³ los pesos quedan diminutos (`w₁ ~ 1e-18`) y arruinan el condicionamiento del
QUBO/QAOA. Por eso el builder trabaja con **normalización interna**: dividir `Q` (y `const`) por una
escala de referencia (p.ej. `Δu²` o `max|Q|`) para que los coeficientes sean O(1) al resolver, y
**deshacer la escala al reportar** energías/SRS. SI para la ciencia, normalizado para el solver.

---

## 4. `FalconConfig` (qué hace configurable al builder)

```python
instance:    {T, L, levels}                 # levels desde Δu
encoding:    "onehot" | "domainwall" | "binary"
terms:       {c_dev: bool, c_smooth: bool, c_crit: "soft" | "deficit_slack"}
constraints: {onehot: "penalty" | "xy_mixer",
              balance: "soft" | "slack",
              release_nonneg: "prohibit" | "slack",
              storage_bounds: "drop" | "postselect" | "slack"}  # default "drop" (ver abajo)
penalties:   {P_onehot, P_bal, P_R, P_crit, ...}   # auto-escala desde J_scale (§8)
normalize:   "delta_u2" | "maxabs" | None
```

Los términos de costo consumen solo "una expresión lineal de `u(t)`" que provee `falcon_encodings.py`,
así cambiar el encoding toca solo 3 funciones, no los términos.

**Defaults recomendados (justificados por datos, ver `docs/ANALISIS_DP_Y_RESULTADOS.md` §8):**
- `storage_bounds = "drop"`: `0 ≤ S ≤ S_max` **nunca ata** -> no se codifica (0 qubits). Reactivar
  solo si cambian datos/constantes y el storage se acerca a las cotas.
- `balance = "slack"` (**preferir el exacto de 1 slack**, `M+s=M_cap`): cuesta solo `≈log₂(0.8·T)`
  qubits (2-6 según instancia), es correcto y casi gratis. El `"soft"` `P_bal(Σu)²` (0 qubits) es
  aproximado (sesga `Σu→0`); usarlo solo como atajo de MVP rápido.
- `c_crit = "soft"` para el MVP (0 qubits, aproxima); `"deficit_slack"` (Opción B) solo en la versión
  fiel (agrega `~2·T·b` qubits, el más caro).
- `release_nonneg = "prohibit"` (nivel prohibido / penalización lineal, 0 slacks).
- `onehot = "penalty"` en MVP; `"xy_mixer"` como variante (Fase 2/3).

---

## 5. Arquitectura (módulos planos en `scripts/` + firmas)

```
scripts/
  falcon_config.py    # dataclass FalconConfig
  falcon_constants.py # load_official_constants(); instance_params(...); compute_weights(...)
  falcon_data.py      # load_ibwc_csv(path); build_weekly_benchmark(data_dir) -> DataFrame (m³, 52 sem)
  falcon_storage.py   # simulate_storage(S0, dS, u); compute_release(R_obs, u); check_constraints(...)
  falcon_srs.py       # compute_costs(S, u, Smin); compute_srs(costs, weights)
  falcon_encodings.py # build_var_index(T,L); linear_expr_u(t); linear_expr_storage(t)
  falcon_qubo.py      # add_square_of_linear_expression(...); build_qubo(cfg,data); qubo_energy(...); to_quadratic_program(Q)
  falcon_solvers.py   # exhaustive(...); simulated_annealing(...); decode_and_verify(samples,...)
  falcon_baselines.py # historical(...); threshold_rule(..., clamp_release=False); dp_optimal(...)  # DP exacto (§7)
  falcon_results.py   # [HECHO] record_run(...): esquema fijo -> results/runs/{run_id}.json + results/runs_summary.csv
  falcon_run_baselines.py # [HECHO] corre+registra los 3 baselines en small/medium/large; deja SRS_hist (ref) guardada
  falcon_qaoa.py      # (Fase 3) build/run/decode QAOA; flag de device CPU/GPU/MPS
```

---

## 6. Checklist por fases

### Fase 0: fundamentos

- [x] **Entorno (mínimo)**: `.venv` con `pandas numpy pytest` (Python 3.14). *Done.* El stack cuántico
  (`qiskit, dimod`) + pyright LSP se instalan en Fase 1/3 (posible Python 3.11/3.12 por wheels).
- [x] **`falcon_constants.py`** :: `load_official_constants() -> dict` (lee el JSON, m³);
  `instance_params(weekly_release, T, L) -> {delta_u, u_max, levels, delta_u_full_ref}`
  (Δu = 0.25·mediana sobre T semanas); `compute_weights(T, S_min, u_max) -> {w1,w2,w3}`. *Done:*
  `S_min==0.25·S_max` (822,181,500 m³); pesos y niveles validados; registra Δu de ventana completa.
- [x] **`falcon_data.py`** :: `load_ibwc_csv(path)` (unidad desde header `Value (...)`);
  `build_weekly_benchmark(data_dir) -> DataFrame[week, week_start, week_end, S_obs_m3, R_obs_m3_week,
  DeltaS_obs_m3, DeltaS_source]`. Reglas: S close-of-week; R por integración nativa 15-min; bins de
  7 días anclados, descartar parcial → **52 semanas**; todo en m³. *Done:* 52 filas; validación contra
  CSV de ivan (merge por semana, |Δ| < tol; documentar diferencias por bordes parciales).
- [x] **`falcon_storage.py`** :: `simulate_storage(S0, deltaS, u) -> S[len T+1]`;
  `compute_release(R_obs, u)`; `check_constraints(R_obs, u, S, Smax, umax, eta_B) -> {feasible,
  violations}`. *Done:* recursión `S(t+1)=S(t)+ΔS−u` coincide con ejemplo a mano.
- [x] **`falcon_srs.py`** :: `compute_costs(S, u, Smin) -> {Ccrit, Cdev, Csmooth}`
  (Ccrit=Σmax(0,Smin−S)², Cdev=Σu², Csmooth=Σ(u−u₋₁)²); `compute_srs(costs, weights) -> float`.
  *Done:* coincide con cálculo a mano en ejemplo de 3 semanas. Reusa patrón `evaluate_srs` de ivan.
- [x] **`falcon_baselines.py`** :: `historical(T)`; `threshold_rule(R_obs, deltaS, S0, Smin, Δu)`;
  `dp_optimal(data, params) -> {u*, SRS*}` (DP exacto §7). *Done:* `dp_optimal` ≥ SRS de histórico y
  umbral; `dp_optimal == brute_force` en T=3 (L=3 y L=5); T=52 en ~80 ms.
- [x] **Smoke Fase 0** `scripts/falcon_smoke_fase0.py`: con `S_min` oficial el SRS es no trivial
  (hist −0.296, thr −0.300, dp −0.294, feasible) y DP domina a los baselines. `pyrightconfig.json`
  agregado (extraPaths=scripts).

### Fase 1: el builder QUBO (centro)

- [x] **`falcon_config.py`** :: dataclass `FalconConfig` con defaults decididos (§4). *Done.*
- [x] **`falcon_encodings.py`** (one-hot) :: `build_var_index(T,L)`; `linear_expr_u(t)` (dict
  `{var: coef}` + const); `linear_expr_storage(t)`; `add_exprs`/`scale_expr`. *Done:* roundtrip exacto.
- [x] **`falcon_qubo.py`** :: `add_square_of_linear_expression(Q, const_ref, expr, weight)`
  (`(c+Σaᵢxᵢ)² = c² + Σ(2caᵢ+aᵢ²)xᵢ + 2Σ_{i<j}aᵢaⱼxᵢxⱼ`); `build_qubo(cfg, data) -> (Q, const,
  var_index)`. **Set MVP con los defaults decididos (§4)**: `C_dev + C_smooth + soft-storage C_crit +
  one-hot penalty + release-prohibition P_R`; **`storage_bounds="drop"`** (no se codifica, nunca ata);
  **balance**: primer corte `soft P_bal(Σu)²` para validar el pipeline, luego cambiar a
  **`balance="slack"` exacto de 1 slack** (`M+s=M_cap`, ~log₂ qubits, preferido). Normalización interna
  (§3); `qubo_energy(Q, x, const)`; `to_quadratic_program(Q)`. *Done:* construye T12/L3 sin error;
  nº vars = `T·L` (+ ~log₂ del slack de balance si exacto); sin slacks de storage.
- [x] **Gate de energía** (`falcon_energy_gate.py`): `xᵀQx·scale + const == J == −SRS` en 9 bitstrings
  factibles, T12/L3 y T26/L5. *Done:* max |Δ| ~1e-13.
- [x] **`falcon_solvers.py`** :: `exhaustive_qubo` (vectorizado `X@Q`) para T12/L3; `simulated_annealing_qubo`;
  `decode_and_verify` (post-selección). *Done:* `exhaustive_qubo == dp` en T12/L3; SA factible en T26/L5
  (gap ~0.03, heurístico a afinar). Registro en `results/` via `falcon_run_qubo.py`
  (`qubo_exhaustive`, `qubo_sa`).

### Fase 2: crecer el modelo (todos los trucos)

- [ ] `C_crit` Opción B (slacks déficit/superávit) en `falcon_qubo.py`.
- [ ] Restricciones exactas con slacks (balance, cotas de storage).
- [ ] Encoding **domain-wall** (`T·(L−1)` bits) en `falcon_encodings.py`.
- [x] Encoding **binary/log** (`T·⌈log₂L⌉` bits; necesario para que la instancia chica entre en
  statevector, §10). *Done (2026-07-01):* `falcon_encodings.py::BinaryVarIndex` + `build_binary_var_index`;
  `build_qubo` despacha por encoding (validez=penalizar codewords inválidos `l≥L`, release-prohibit y
  balance agnósticos vía `balance_M_expr`/`add_codeword_penalty`); exacto para `b≤2` (L≤4). Validado:
  roundtrip, gate energía==−SRS (`3e-14`), y **mínimo global == DP\*** en T5/L3. T12/L3 = 24 qubits (vs
  36 one-hot). Sin regresión en one-hot (gate + `exhaustive==dp` OK). Ver `docs/ANALISIS_DP_Y_RESULTADOS.md` §9.
- [ ] **XY-mixer** para eliminar la penalización one-hot (Fase 3 lo usa).
- [ ] *Done de cada variante:* reproduce el mismo SRS para el mismo `u` decodificado.
- [x] **E2 (clásico) - chunking temporal**: `scripts/falcon_chunked.py::staged_solve(..., solver="dp")`.
  Parte T en bloques de `block_size`; secuencial: `S0_blk`=storage final del bloque previo,
  `k_prev_init`=último nivel del previo (linkea `C_smooth`); balance por bloque `balance_split ∈
  {"eta_local","global_greedy"}`. Concatena `u`, evalúa SRS/factibilidad global con los módulos
  canónicos, y `gap_vs_full` vs `dp_optimal` global. *Done (2026-07-01):* `dp_optimal` extendido con
  `k_prev_init`; `build_qubo` con `u_prev` (frontera C_smooth, validado Δenergía==w3·(u₀−u_prev)²).
  **Hallazgo T26/L5 bloques de 5:** DP-chunked (eta_local) = −0.311534 = histórico, con **gap_vs_full
  = 0.021**: el troceo con balance por-bloque **NO alcanza el óptimo global** (full DP −0.290423 SÍ es
  no trivial) porque la redistribución óptima cruza fronteras de bloque dentro del η global. → probar
  `global_greedy` / bloques más grandes. Ver `docs/ANALISIS_DP_Y_RESULTADOS.md` §10.

### Fase 3: cuántico + escalamiento

- [x] **`falcon_qaoa.py`** :: `Q → to_quadratic_program → Ising → QAOA` en simulador, instancia chica,
  CPU local → GPU/MPS en WCentroid. Reusar `docs/georgia_qubo_snippets.md`
  (`build_quadratic_program`, `precompute_diagonal`, `build_qaoa_circuit`). Convenciones: seed 42,
  ≥5 restarts si `Q` mal condicionada, schema JSON unificado. *Done (2026-07-01):* **Qiskit + Aer
  statevector** (`.venv-quantum`, py3.12); circuito manual RZ/RZZ/RX; `<H>`+muestreo desde statevector
  + diagonal precomputada; COBYLA, seed 42, ≥5 restarts, transpile-once. Gate energía==−SRS `~1e-15`;
  `H_ising`==`diag−const` `1e-14`. Decodifica + verifica factibilidad (top-256). Ver `docs/…RESULTADOS.md` §9.
- [~] **Benchmark/escala**: tablas ΔSRS vs baselines (DP, histórico, umbral), runtime y escala
  T12→T26→T52 → `results/`. *Parcial (2026-07-01):* `falcon_run_qaoa.py` corre y registra debug T5/L3
  (one-hot, 17q, **AR 1.000**) y small T12/L3 (binary, 24q, AR 1.31, factible). Tabla en §9. **Falta**
  T26/T52 (no entran en statevector → MPS/sampling o chunking E2). Calibración de penalties (§8):
  hecha, no es la palanca (maxabs absorbe el multiplicador); la palanca es profundidad/XY-mixer.
- [ ] **E1 - robustez de ventana**: `falcon_run_windows.py` (o extensión del runner) corre cada método
  (historical, threshold pure/clamped/balanced, dp, brute si enumerable) en **3 ventanas por instancia**
  (excepto T52 que ocupa todo el año): `first` (start=0), `middle` (start=(52-T)//2), `stress` (auto:
  ventana de T semanas con menor storage medio o mayor varianza). Agregar campo `window_start` /
  `window_label` al esquema de `falcon_results.py` y al `run_id` (evita colisiones). *Done:* tabla que
  muestra la variación de SRS/ΔSRS/factibilidad entre ventanas; concluir si `first` es representativa.
- [x] **E2 (cuántico) - chunking en QUBO/QAOA**: `staged_solve(..., solver="qaoa")` cambia el solver
  por-bloque (DP → nuestro `falcon_qaoa.run_qaoa`); cada bloque = QUBO chico que entra en statevector.
  *Done (2026-07-01):* runner `falcon_run_chunked.py`, `method="qaoa_chunked"`. **T26/L5 bloques de 5 =
  25 qubits/bloque** corre **QAOA REAL** (no el fallback silencioso de Ivan): 4/6 bloques factibles, 2
  reportados infactibles explícitamente. QAOA-chunked = −0.356259 (**infactible**, gap_vs_full 0.066):
  peor que DP-chunked por el límite p=1 + balance soft por-bloque que deriva el balance global.
  Descomposición de gap (chunking 0.021 + QAOA 0.045). Honesto y consistente con spec §7 (sin ventaja
  cuántica). Ver §10. Mejorar: `global_greedy`, mayor p (init INTERP), XY-mixer, balance duro por bloque.
- [ ] **Writeup (rúbrica 1 y 5)**: formulación del problema + SDG 6.4/6.5/13.1 + impacto + literatura;
  pros/cons; nota de uso de IA; slides en inglés. *Done:* cubre los 5 criterios de la rúbrica.

### Opcionales / baselines extra (no bloquean el flujo QUBO)

- [ ] **`falcon_milp.py` :: `milp_optimal(...)` (OPCIONAL)**: óptimo exacto **independiente** para
  cross-check del DP, sobre todo en **medium (T26/L5)** donde brute force es imposible (`5^26≈1.5e18`).
  Modelo: `x_{t,ℓ}` binario one-hot (`Σ_ℓ x_{t,ℓ}=1`), `u_t=Σ_ℓ aₗ x_{t,ℓ}`, storage lineal en bits;
  `C_crit` linealizado con déficit `d_t≥0, d_t≥S_min−S_t` (min Σ d_t²  →  MIQP, o aprox. lineal Σ d_t
  si el solver es lineal); `C_dev/C_smooth` cuadráticos (MIQP) o linealizados; restricciones R≥0,
  `0≤S≤S_max`, balance `|Σu|≤B`. Registrar con `record_run(method="milp")`. Dependencia:
  `python-mip`/`PuLP+CBC` (lineal) o solver MIQP (posible Python 3.11/3.12). *Done:* `milp == dp` en
  medium (o explica la diferencia si se usa aproximación lineal de C_crit). **Sensibilidad, reportar
  aparte del score oficial.**
- [ ] **Regla de umbral balanceada** (de ivan) y **annealing quantum-inspired** como baselines extra
  comparables (registrar con `record_run`). Opcionales.

---

## 7. DP exacto sobre lattice (ground truth para todas las instancias)

Como `u` es discreto, el storage vive en un **lattice entero exacto**: `S_t = H_t − Δu·c_t`, con
`c_t = Σ_{k<t} kₖ` (suma entera acumulada, `k∈{−2..2}`). Estado DP = `(t, c_t, k_{t−1})`
(`k_{t−1}` solo para `C_smooth`); `c_t∈[−2t,2t]`. Estados `O(T²·L)`, transiciones `O(T²·L²)` →
trivial incluso en T=52. Todos los costos y restricciones (`C_crit`, `C_dev`, `C_smooth`, balance
`|Δu·c_T|≤B`, no-negatividad, cotas de storage) son funciones **exactas** de los enteros: **sin
discretización, sin error**. Por eso el DP da el **óptimo global exacto del problema discretizado en
TODAS las instancias** → es a la vez baseline fuerte y ground truth para validar QUBO/QAOA a cualquier
tamaño (mejor que brute force, que muere en T>~6).

---

## 8. Calibración de penalizaciones

Por capas: (1) **auto-escala desde `J_scale`** (spec §14:
`J_scale ≈ |w₁C_crit^hist|+|w₂C_dev^rule|+|w₃C_smooth^rule|`; `P_onehot=10·J_scale`,
`P_bal=J_scale/max(B²,eps)`, `P_R=10·J_scale`, `P_crit=10·J_scale/storage_scale²`); (2) **barrido de
sensibilidad** `P∈{0.1,1,10,100}·base` **una vez en la instancia chica**, fijar multiplicadores y
reusar en instancias grandes. Reportar SRS, factibilidad, violaciones one-hot/balance, storage mínimo,
runtime. Combina con la normalización interna (§3).

---

## 9. Reuso (referencia + regresión, no dependencia)

- **ivan** (`notebooks/ivan/`): primitiva `add_square_linear`, `official_parameters`, `evaluate_srs`,
  decoder tolerante a one-hot, estrategia de slacks por expansión binaria. Reimplementar limpio en
  `scripts/`, no importar del notebook.
- **emilio** (`notebooks/emilio/`): términos `Q_crit/Q_dev/Q_smooth/Q_onehot`, validadores brute-force,
  SA con filtrado de factibilidad.
- **Georgia** (`docs/georgia_qubo_snippets.md`): QUBO→Qiskit, `precompute_diagonal`, QAOA manual,
  verificación de energía, convenciones (seed 42, ≥5 restarts, schema JSON).

Coordinar con ivan/emilio antes de que el builder compartido reemplace sus piezas (GUIDELINES §3).

---

## 10. Notas de hardware (ver `CLAUDE.md` "Compute environments")

Statevector tope ~30 qubits (T4 16 GB). T12/L3: one-hot=36 (no entra) pero **binary/domain-wall=24
(entra)** → los encodings compactos son necesarios para QAOA exacta en simulador. Medio/grande →
MPS/sampling. Device tras un flag: CPU (M4) local; GPU/MPS (T4, WCentroid) para lo pesado.

---

## 11. Preguntas abiertas (iterar)

- ¿Convención exacta del borde de semana sobre nuestros datos (anclaje al primer timestamp vs ISO)?
  Decidido: anclado + descartar parcial → 52; reconfirmar al validar contra ivan.
- ¿Mantener además un `Δu` de referencia fijo para sanity cross-instancia? Decidido: registrar, no
  optimizar con él.
- Baseline "informado por literatura" (rúbrica 20%): ¿DP exacto cuenta como state-of-the-art o sumamos
  MILP/SA citando refs [1-3]?
- Granularidad/ábaco del writeup de impacto social (rúbrica 25%): definir antes de la última semana.
