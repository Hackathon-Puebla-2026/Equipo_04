# CLAUDE.md — Falcon Challenge: EDA, Benchmark Clásico y Formulación QUBO/QAOA

Input de trabajo para implementar el pipeline QUBO/QAOA del reto Falcon. **Guía/ejemplo, no fuente de
verdad** (la spec `docs/FalconChallenge_V6.md` manda). Objetivo: (1) analizar datos del embalse,
(2) construir el dataset semanal, (3) baselines clásicos, (4) formular el QUBO, (5) resolver clásico
y luego QAOA. Basado en los notebooks de referencia `QUBO_Mathematical_Definition`,
`QUBO_PenaltyMethod`, `QUBO_Examples_MaximumCut`.

> Constantes ya resueltas en `data/falcon_reservoir_constants.json`: `S_max=3,288,726 TCM`,
> `S_min=0.25·S_max≈822,181.5 TCM`. Código QUBO/QAOA reutilizable buffereado en
> `docs/georgia_qubo_snippets.md`.

---

## 1. Reglas de seguridad sobre datos

1. **No modificar/borrar/renombrar/mover** archivos originales en `data/` o `data/source_csv/`.
2. No descargar datos sin autorización. No usar `.htm/.gif/.json/.txt/.pdf` como sustituto de las
   series numéricas.
3. Procesados → solo `data/processed/` y `results/`. Código nuevo → `src/`. (Ojo: ubicación real de
   archivos la define `docs/GUIDELINES.md`.)
4. Si falta un dataset oficial, reportarlo y marcar la aproximación como **preliminar / no oficial**.
5. No `git add/commit/push` ni comandos destructivos. No instalar paquetes sin permiso.
6. No presentar resultados sintéticos/derivados como benchmark oficial.

---

## 2. Objetivo del reto

Optimizar el ajuste semanal de liberación `u(t)` del International Falcon Reservoir. Release optimizado
`R(t)=R_obs(t)+u(t)` (con `u=0` → operación histórica). Dinámica
`S_opt(t+1)=S_opt(t)+ΔS_obs(t)−u(t)`. Maximizar `SRS = −(w₁C_crit + w₂C_dev + w₃C_smooth)`;
para QUBO minimizar `J = w₁C_crit + w₂C_dev + w₃C_smooth + penalizaciones`.
(Detalle de variables y dinámica: ver `CLAUDE.md` raíz y `docs/FalconChallenge_V6.md`.)

---

## 3. Términos de costo del benchmark

### 3.1 Almacenamiento crítico
`C_crit = Σ_{t=0}^{T} [max(0, S_min − S_opt(t))]²` — penaliza solo semanas bajo `S_min`.

### 3.2 Desviación vs histórico
`C_dev = Σ_{t=0}^{T-1} u(t)²` — evita cambios excesivos vs descarga histórica.

### 3.3 Cambios bruscos entre semanas
`C_smooth = Σ_{t=1}^{T-1} [u(t) − u(t-1)]²` — favorece ajustes suaves.

---

## 4. Restricciones oficiales

`R(t) ≥ 0` · `|u(t)| ≤ u_max` · `0 ≤ S_opt(t) ≤ S_max` ·
`|Σ_{t} u(t)| ≤ η·Σ_{t} R_obs(t)`.

La última evita "mejorar" el almacenamiento reduciendo sistemáticamente liberaciones: obliga a
redistribuir en el tiempo, no a retener agua.

---

## 5. Parámetros oficiales del benchmark

`L=5`, `u(t) ∈ {−2Δu, −Δu, 0, Δu, 2Δu}`, `Δu = 0.25·R̃ᵒᵇˢ_week` (mediana semanal), `u_max = 2Δu`,
`η = 0.10`, `S_min = 0.25·S_max`, `S_scale = S_min`. Pesos:
`w₁ = 1/((T+1)·S_scale²)`, `w₂ = 0.1/(T·u_max²)`, `w₃ = 0.1/((T-1)·(2u_max)²)`.
(También en `CLAUDE.md` raíz.)

---

## 6. Instancias de escalamiento

| Instancia | T | L | Uso |
|---|---:|---:|---|
| Small | 12 | 3 | validación / debugging |
| Medium | 26 | 5 | **benchmark oficial** |
| Large | 52 | 5 o 7 | escalamiento |

Agendas candidatas: `N_schedules = Lᵀ`.

---

## 7. Datasets requeridos

### 7.1 Oficiales mínimos

| Variable | Dataset | Estación | Uso |
|---|---|---|---|
| `S_obs(t)` | `Total Storage.Web-Daily-tcm@08461200` | 08461200 | almacenamiento |
| `ΔS_obs(t)` | `Discharge.Total.Change-in-Storage@08461200` | 08461200 | cambio de almacenamiento |
| `R_obs(t)` | `Discharge.Best Available@08461300` | 08461300 | release histórico |
| `S_max` | Falcon total conservation capacity | overview | define `S_min` (ya resuelto, ver constants.json) |

### 7.2 CSV anexos disponibles

En `data/` hay exports de `Discharge@08461300` (m³/s, convertir a volumen: `R_week = R_{m³/s}×604800`),
`Total Storage@08461200` (**ojo: header dice `Value (m^3)` aunque el dataset se llame `tcm`** — leer
unidad real del header), `Reservoir Elevation` y `Lake Area` (solo EDA, no indispensables al QUBO).

### 7.3 Faltante importante

No está el CSV oficial `Discharge.Total.Change-in-Storage@08461200`. Si no aparece, derivar
`ΔS_obs(t) ≈ S_obs(t+1) − S_obs(t)` y etiquetar `derived_from_total_storage_not_official` (valida el
pipeline pero **no es benchmark oficial**).

---

## 8. Dataset semanal esperado

Construir `data/processed/falcon_weekly_benchmark.csv` con columnas:
`week, week_start, week_end, S_obs_m3, DeltaS_obs_m3, R_obs_m3_week, DeltaS_source`.

Agregación: `S_obs_m3` = último valor semanal (cierre de semana); `R_obs_m3_week` = **volumen acumulado
de la semana** = `Σ(caudal × Δt)` dentro de la semana (equivalente a `promedio_caudal × 604800` si el
muestreo es uniforme); `DeltaS_obs_m3` = oficial si existe (a volumen semanal), si no diferencia semanal
de `S_obs_m3`. Alinear por semana y descartar semanas con faltantes críticos.

> **Acumulado vs mediana — no confundir.** `R_obs(t)` (la *serie* semanal que entra al balance de
> almacenamiento, a `R(t)≥0` y a `|Σu|≤η·ΣR_obs`) es siempre el **volumen acumulado** por semana (suma,
> nunca mediana). La **mediana** aparece solo en el *escalar* `R̃ᵒᵇˢ_week` de §5 (`Δu = 0.25·R̃ᵒᵇˢ_week`):
> es la mediana **entre semanas** de esos volúmenes ya acumulados (mediana, no media, por robustez ante
> semanas de crecida). Pipeline: caudal 15-min → suma por semana → `R_obs(t)`; luego mediana entre
> semanas → `R̃ᵒᵇˢ_week` → `Δu`.

---

## 9. Baselines clásicos

### 9.1 Histórico
`u_hist(t)=0`, `R_hist(t)=R_obs(t)`.

### 9.2 Regla de umbral
`u_rule(t) = −Δu` si `S_rule(t) < S_min`, si no `0`; `R_rule(t)=R_obs(t)+u_rule(t)`.

### 9.3 Clásico fuerte (≥1)
Búsqueda exhaustiva `T=12,L=3` (`3¹²=531,441`, factible), simulated annealing, dynamic programming,
MILP/MIQP, o evolutionary. Para `T=26,L=5` la exhaustiva no es práctica (`5²⁶≈1.49×10¹⁸`).
Código reutilizable (brute force, SA, GA) en `docs/georgia_qubo_snippets.md`.

---

# 10. Principios QUBO (de los notebooks de referencia)

**10.1 Definición.** Minimizar `xᵀQx = Σᵢ Q_ii x_i + Σ_{i<j} Q_ij x_i x_j`, `x ∈ {0,1}ⁿ`. Lineales en
la diagonal, cuadráticos fuera; usar triangular superior. Como `xᵢ²=xᵢ`, un lineal va en la diagonal.
Maximizar → `max f = min −f`.

**10.2 Método de penalización.** Restricción → `f(x) + Σ Pᵢ gᵢ(x)`, con `gᵢ=0` si se cumple, `>0` si se
viola. Igualdad `Σ aᵢyᵢ = b` → `P(Σ aᵢyᵢ − b)²`. Desigualdad `Σ aᵢyᵢ ≤ b` → slack `≥0`:
`P(Σ aᵢyᵢ + s − b)²` (slack a binario).

**10.3 Entero → binario.** Para `ȳ_i ≤ y_i ≤ Ȳ_i`: `y_i = ȳ_i + Σ_{j=0}^{N-2} 2ʲ x_jⁱ +
(Ȳ_i − ȳ_i − Σ 2ʲ) x_{N-1}ⁱ`, con `N = ⌈log₂(Ȳ_i − ȳ_i + 1)⌉`.

**10.4 Analogía Max-Cut.** Arista `(i,j)`: corte `= x_i + x_j − 2x_i x_j`; maximizar → minimizar
`−x_i − x_j + 2x_i x_j`. Patrón para Falcon: variables binarias por decisión → escribir `u(t)`, `S(t)`
y restricciones como expresiones lineales en binarias → elevar al cuadrado → cargar en `Q`.

---

# 11. Formulación QUBO recomendada para Falcon

## 11.1 Variable de decisión one-hot
`x_{t,ℓ} ∈ {0,1}`, `x_{t,ℓ}=1 ⟺ u(t)=aₗ`. Niveles: oficial `L=5` `aₗ∈{−2Δu,−Δu,0,Δu,2Δu}`; small
`L=3` `aₗ∈{−Δu,0,Δu}`. One-hot: `Σ_ℓ x_{t,ℓ}=1`, entonces `u_t = Σ_ℓ aₗ x_{t,ℓ}`. Bits de decisión
`N = T×L` (T12/L3 → 36; T26/L5 → 130).

## 11.2 Penalización one-hot
`P_onehot Σ_t (Σ_ℓ x_{t,ℓ} − 1)²`. Expansión por semana:
`(Σ_ℓ x_{t,ℓ} − 1)² = 1 − Σ_ℓ x_{t,ℓ} + 2 Σ_{ℓ<m} x_{t,ℓ} x_{t,m}` (la constante se ignora en `Q`).
→ diagonal `−P_onehot` por bit; off-diagonal `+2P_onehot` por par `(ℓ,m)`, `ℓ<m`.

## 11.3 Expresión lineal de almacenamiento
`H_t = S_0 + Σ_{k=0}^{t-1} ΔS_obs(k)` (histórico acumulado sin ajustes). Entonces
`S_t = H_t − Σ_{k<t} u_k = H_t − Σ_{k<t} Σ_ℓ aₗ x_{k,ℓ}`. **Lineal en las binarias** → cualquier costo
cuadrático en `S_t` es QUBO-able.

## 11.4 `C_dev`
`w₂ Σ_t (Σ_ℓ aₗ x_{t,ℓ})²`. Con one-hot se simplifica a `w₂ Σ_t Σ_ℓ aₗ² x_{t,ℓ}` (diagonal). Para
robustez usar el helper general (§12).

## 11.5 `C_smooth`
`w₃ Σ_{t=1}^{T-1} (Σ_ℓ aₗ x_{t,ℓ} − Σ_ℓ aₗ x_{t-1,ℓ})²`. Genera **acoplamientos entre semanas
consecutivas**.

## 11.6 `C_crit`: dos opciones
El oficial usa hinge `max(0, S_min − S_t)²`, no directamente cuadrático.

**Opción A — MVP (no oficial).** Costo suave `C_storage^soft = Σ_t (S_t − S_target)²`, con
`S_target = S_min` o superior. Simple, sin auxiliares, útil para validar `Q`/QAOA/pipeline; pero
penaliza también estar por arriba de `S_min` → reportar como aproximación.

**Opción B — fiel (déficit/superávit).** Auxiliares `D_t,E_t ≥ 0` con `S_t + D_t − E_t = S_min`
(`D_t`=déficit, `E_t`=superávit). Minimizando `D_t²` con la igualdad, `D_t = max(0, S_min − S_t)`.
Agregar `w₁ Σ_t D_t² + P_crit Σ_t (S_t + D_t − E_t − S_min)²`. Codificar
`D_t = q_S Σ_r 2ʳ d_{t,r}`, `E_t = q_S Σ_r 2ʳ e_{t,r}` (`q_S`=resolución en m³; bits cubren el rango;
más qubits pero más fiel). Recomendación: A para validar, B para la versión final, comparar.

## 11.7 Descarga no negativa (`R_t = R_obs,t + u_t ≥ 0`)
**Simple:** si `R_obs,t + aₗ < 0`, prohibir ese nivel con `P_R x_{t,ℓ}`. **Slack:**
`R_obs,t + u_t − r_t = 0`, `r_t ≥ 0` binario (más qubits, no recomendado en MVP).

## 11.8 Balance acumulado (`|Σ_t u_t| ≤ B`, `B = η Σ_t R_obs,t`)
**Aproximación QUBO:** `P_bal (Σ_t u_t)²` (no es la desigualdad exacta pero favorece balance ~0).
**Exacta con slacks:** dos desigualdades `±Σ_t u_t ≤ B` con `η₊,η₋ ≥ 0`:
`P_bal (Σ u_t + η₊ − B)² + P_bal (−Σ u_t + η₋ − B)²` (slacks binarios). Recomendación: MVP usa
`P_bal(Σu)²` + filtrar factibilidad al decodificar; final usa slacks si los qubits alcanzan.

## 11.9 Almacenamiento físico (`0 ≤ S_t ≤ S_max`)
**MVP:** no incluir en QUBO; simular `S_t` tras muestrear, rechazar infeasibles, reportar tasa de
factibilidad. **Penalizado:** `S_t + s_t^up = S_max` y `−S_t + s_t^low = 0` (slacks `≥0` binarios,
penalizaciones cuadráticas; muchos qubits auxiliares).

---

# 12. Helper recomendado para construir QUBO

```python
def add_square_of_linear_expression(Q, offset, expr, weight):
    """
    Adds weight * (c + sum_i a_i x_i)^2 to Q.
    expr: {'constant': c, 'linear': {var_index: coefficient}}
    For binary x_i, x_i^2 = x_i. Expansion:
        weight * [ c^2 + sum_i (2*c*a_i + a_i^2) x_i + sum_{i<j} 2*a_i*a_j x_i*x_j ]
    """
```

`(c + Σᵢ aᵢxᵢ)² = c² + Σᵢ (2c·aᵢ + aᵢ²) xᵢ + 2 Σ_{i<j} aᵢaⱼ xᵢxⱼ`. Usarlo para: one-hot, `C_dev`,
`C_smooth`, soft storage, igualdad déficit/superávit, balance, y constraints con slack.

---

# 13. QUBO objective recomendado

## 13.1 MVP
`J_MVP = w₂C_dev + w₃C_smooth + λ_S Σ_t (S_t − S_target)² + P_onehot Σ_t (Σ_ℓ x_{t,ℓ} − 1)²
+ P_bal (Σ_t u_t)² + P_R Σ_{(t,ℓ)∈I_R} x_{t,ℓ}`, donde `I_R` = niveles inválidos con
`R_obs,t + aₗ < 0`. Ruta más rápida a un prototipo QUBO/QAOA.

## 13.2 Official-compatible
`J = w₁ Σ_t D_t² + w₂C_dev + w₃C_smooth + P_crit Σ_t (S_t + D_t − E_t − S_min)²
+ P_onehot Σ_t (Σ_ℓ x_{t,ℓ} − 1)² + P_bal·(balance) + P_R·(release) + P_S·(storage bounds)`.
Etiquetar **official-compatible** solo si: `ΔS_obs` viene del dataset oficial; `S_max` es la capacidad
oficial; pesos oficiales; constraints aplicadas/rechazadas consistentemente; unidades compatibles.

---

# 14. Elección de coeficientes de penalización

`P` debe ser grande para evitar violaciones, pero no tanto que opaque diferencias entre soluciones
factibles. Estrategia: (1) escala del objetivo
`J_scale ≈ |w₁C_crit^hist| + |w₂C_dev^rule| + |w₃C_smooth^rule|`; (2) probar
`P_onehot=10·J_scale`, `P_bal=J_scale/max(B²,eps)`, `P_R=10·J_scale`,
`P_crit=10·J_scale/storage_scale²`; (3) barrido `P ∈ {0.1,1,10,100}·base`; (4) reportar SRS,
factibilidad, violaciones one-hot/balance, almacenamiento mínimo, runtime.

---

# 15. Conversión QUBO → Ising / QAOA

`x_i = (1 − z_i)/2`, `z_i ∈ {−1,1}`. Pasos: construir `Q` (triangular superior) → `QuadraticProgram`
de Qiskit Optimization → `to_ising()` → QAOA en simulador para `T=12,L=3` → comparar contra exhaustiva
y SA. **No empezar con hardware real.** Conversión `Q→QuadraticProgram` y circuito QAOA manual (evita
cuelgues de `PauliEvolutionGate`) buffereados en `docs/georgia_qubo_snippets.md`.

---

# 16. Estructura de archivos esperada

```text
FalconChallenge/
├── data/{source_csv/<exports>, processed/falcon_weekly_benchmark.csv}
├── src/{data_loader, eda_falcon, storage_model, srs_score, classical_baselines,
│        qubo_utils, qubo_falcon, solve_classical, solve_qaoa}.py
├── results/{data_quality_report.md, benchmark_summary.csv, benchmark_timeseries.csv,
│            qubo_summary.md, qaoa_results.csv, figures/}
└── CLAUDE.md
```

(La ubicación real en este repo la define `docs/GUIDELINES.md`: código compartido en `scripts/` raíz,
carpetas por integrante. Esta estructura es la sugerida por la guía, adaptar.)

---

# 17. Implementación mínima esperada (firmas)

**17.1 `data_loader.py`:** `load_ibwc_export_csv(path)`, `identify_dataset(path)`,
`build_weekly_benchmark(input_dir, output_path)`. Manejar header especial (`#Data Set Export…` +
`Timestamp, Value (...)`); renombrar columnas a `timestamp, value`.

**17.2 `storage_model.py`:** `simulate_storage(S0, deltaS_obs, u)`, `compute_release(R_obs, u)`,
`check_constraints(R_obs, u, S, Smax, umax, eta)`.

**17.3 `srs_score.py`:** `compute_costs(S, u, Smin)`, `compute_weights(T, Smin, umax)`,
`compute_srs(costs, weights)`.

**17.4 `classical_baselines.py`:** `historical_policy(T)`,
`threshold_policy(R_obs, deltaS_obs, S0, Smin, delta_u)`,
`exhaustive_search_small_instance(levels, R_obs, deltaS_obs, S0, params)`.

**17.5 `qubo_utils.py`:** `add_linear(Q,i,coeff)`, `add_quadratic(Q,i,j,coeff)`,
`add_square_of_linear_expression(Q, offset, expr, weight)`, `qubo_energy(Q, x, offset=0.0)`,
`make_upper_triangular(Q)`.

**17.6 `qubo_falcon.py`:** `build_variable_index(T, L, aux_config=None)`,
`linear_expr_u(t, levels, var_index)`, `linear_expr_storage(t, S0, deltaS_obs, levels, var_index)`,
`build_falcon_qubo_mvp(data, params, penalties)`,
`build_falcon_qubo_official_like(data, params, penalties, aux_config)`,
`decode_solution(x, var_index, levels, data, params)`.

**17.7 `solve_qaoa.py`:** QUBO → Qiskit `QuadraticProgram`; QAOA small instance; guardar
muestras/soluciones; comparar con baselines. Si Qiskit no está instalado, fallar con mensaje claro.

---

# 18. Salidas esperadas

- `results/data_quality_report.md`: archivos usados, columnas, unidades, rango de fechas, faltantes,
  si existe `ΔS_obs` oficial, si el benchmark es oficial o preliminar.
- `results/benchmark_summary.csv`: `mode, method, T, L, SRS, DeltaSRS_vs_historical,
  DeltaSRS_vs_threshold, Ccrit, Cdev, Csmooth, min_storage, weeks_below_Smin,
  release_balance_error, feasible, runtime_seconds`.
- `results/qubo_summary.md`: nº de binarias, auxiliares, términos lineales/cuadráticos, penalties,
  tipo (`MVP_SOFT_STORAGE` u `OFFICIAL_LIKE_WITH_DEFICIT_SLACKS`), aviso si usa `ΔS_obs` derivado.
- `results/figures/`: `storage_comparison.png`, `release_adjustments.png`, `release_comparison.png`,
  `srs_comparison.png`, `qubo_energy_distribution.png` (si hay muestras QAOA/annealing).

---

# 19. Orden recomendado de trabajo

1. Inspeccionar CSV anexos. 2. Construir `falcon_weekly_benchmark.csv`. 3. EDA + reporte de calidad.
4. Simulación clásica. 5. SRS. 6. Baseline histórico. 7. Baseline regla de umbral.
8. Búsqueda exhaustiva `T=12,L=3`. 9. QUBO MVP. **10. Validar QUBO comparando energía QUBO vs costo
clásico decodificado.** 11. QUBO official-like con déficit/superávit si el tamaño lo permite.
12. QAOA en simulador. 13. Comparar contra baselines. 14. Guardar reportes reproducibles.

---

# 20. Criterios de aceptación

No modifica datos originales; produce dataset semanal procesado; reporta si es oficial o preliminar;
calcula SRS histórico y de regla de umbral; ≥1 optimización clásica para `T=12,L=3`; QUBO reproducible
(matriz o dict); decodifica binarios → `u(t)`; verifica restricciones tras decodificar; compara vs
baselines; reporta runtime y escalabilidad; distingue resultados oficiales / preliminares / sintéticos
/ aproximados por QUBO MVP.

---

# 21. Arranque sugerido

Primero etapas 1-8 (inspección → dataset → EDA → simulación clásica → SRS → baselines → exhaustiva
`T=12,L=3`), sin QAOA todavía. Antes de QAOA, validar que costo clásico y energía QUBO coincidan para
varias soluciones decodificadas. Entregables: `data_quality_report.md`, `benchmark_summary.csv`,
`benchmark_timeseries.csv`, figuras, y si el resultado es oficial o preliminar.

---

# 22. Notas críticas para evitar errores

1. `Discharge.Best Available@08461300` suele estar en m³/s: convertir a volumen semanal antes de restar
   de almacenamiento en m³.
2. El `Total Storage` puede llamarse `tcm` pero el CSV indica `Value (m^3)`: **leer la unidad real del
   header.**
3. Si `ΔS_obs` se deriva de almacenamiento, etiquetar como preliminar.
4. Si se usa `S_max = max(S_obs)` por falta de capacidad oficial, etiquetar preliminar (ya tenemos la
   oficial: ver `data/falcon_reservoir_constants.json`).
5. El QUBO MVP con almacenamiento suave no reproduce exactamente `C_crit`: reportar como aproximación.
6. El QUBO con slacks aumenta mucho los qubits: usarlo primero en `T=12`.
7. **Siempre decodificar y verificar factibilidad después de resolver QUBO/QAOA.**
8. No comparar QAOA solo contra histórico: también contra regla de umbral y exhaustiva/SA cuando se
   pueda.

---

# 23. Resumen ejecutivo

Ruta: `CSV anexos → limpieza semanal → benchmark clásico → QUBO one-hot L=3,T=12 → validación vs
exhaustiva → QUBO L=5,T=26 → QAOA/SA → comparación ΔSRS`. Formulación inicial one-hot: `x_{t,ℓ}=1 ⟺
u_t=aₗ`, `u_t=Σ_ℓ aₗ x_{t,ℓ}`, `S_t = H_t − Σ_{k<t} Σ_ℓ aₗ x_{k,ℓ}`. Minimizar
`J = w₁C_crit + w₂C_dev + w₃C_smooth + penalties`. Primero `C_crit` con aproximación suave para validar;
después déficit/superávit con slacks para acercarse al benchmark oficial.
