# Resumen de hallazgos del equipo (EDA + QUBO preliminar)

Detalle del trabajo de **ivan** y **emilio** hasta ahora, para alinear al equipo antes de construir el
QUBO builder compartido. El digest decisivo está en `docs/HALLAZGOS_CLAVE.md` (siempre cargado); este
archivo es el detalle on-demand.

> Fuentes analizadas: `notebooks/ivan/Falcon_QUBO_Documentacion_Colab(1).ipynb` + `notebooks/ivan/Act2/`
> y `notebooks/emilio/Pruebas.ipynb`. Todos los números son los reportados por esos notebooks/CSV.

---

## 1. Estado general

- **EDA y dataset semanal: hechos** (ivan, carpeta `Act2/`). Datos de 1 año listos y validados.
- **Formulación QUBO one-hot: hecha por ambos** (ivan con slacks exactos; emilio con filtrado de
  factibilidad). Validada con búsqueda exhaustiva en instancias chicas.
- **Solvers: solo clásicos** (simulated annealing). **QAOA: pendiente por todos.**
- **Todo es PRELIMINAR / NO OFICIAL**: ΔS derivado y `S_max` = máx observado (ver §4).

---

## 2. Fundamento de datos (consensuado entre ambos)

5 datasets IBWC (estaciones `08461200`/`08461300`), **2025-06-29 → 2026-06-29 (1 año, 53 semanas)**,
**0 faltantes y 0 duplicados** tras limpieza (`eda_series_summary.csv`):

| Serie | Cadencia | Filas válidas | Unidad | Resumen |
|---|---|---|---|---|
| `R_obs` (Discharge@08461300) | 15 min | 35,023 | m³/s | mín 2.47, mediana **18.87**, máx 71.34, media 21.92 |
| `S_obs` (Total Storage@08461200) | diaria | 366 | m³ | mín **315.0M**, mediana 401.4M, máx **648.84M**, std 69.7M |
| Elevación@08461200 | diaria | 366 | m | 76.91 - 80.31 |
| Lake Area@08461200 | 15 min | 35,027 | m² | 74.4M - 123.8M |
| % Conservation@08461200 | 15 min | 35,027 | % | 9.58 - 19.73 |

16 figuras en `notebooks/ivan/Act2/figures/` (series, histogramas, scatter S-elevación/área, métricas
semanales, baselines).

---

## 3. Dataset semanal (ivan: `Act2/falcon_weekly_benchmark_preliminary.csv`)

**53 semanas**. Columnas:
`week, week_end_date, S_obs_m3, R_obs_m3_week, elevation_m_week_mean, lake_area_m2_week_mean,
DeltaS_obs_m3, DeltaS_source`.

Agregación: `S_obs_m3` = valor diario al cierre de semana; `R_obs_m3_week` = `mean(caudal diario)·604800`
(volumen semanal); `DeltaS_obs_m3` = `S(t+1)−S(t)` **derivado**, etiquetado en todas las filas como
`derived_from_total_storage_not_official`. Semana 1 termina 2025-06-29; semana 53 termina 2026-06-28.

---

## 4. Constantes calculadas y DISCREPANCIAS (lo más accionable)

| Constante | ivan (T=12, L=3) | emilio (T=26, año completo) | Oficial (`falcon_reservoir_constants.json`) |
|---|---|---|---|
| `S_max` | 648.84M m³ (máx obs, **NO oficial**) | 648.84M m³ (máx obs, **NO oficial**) | **3,288.726M m³** (capacidad de conservación) |
| `S_min` = 0.25·S_max | 162.21M m³ | 162.21M m³ | **822.18M m³** |
| `Δu` = 0.25·mediana semanal | 2,688,917 m³ (mediana 10.76M) | 2,953,742 m³ (mediana ~11.8M) | depende de ventana - estandarizar |
| `u_max` = 2·Δu | 5.38M m³ | 5.91M m³ | = 2·Δu |
| `η` | 0.10 | 0.10 | 0.10 |
| `S0` | ~375.6M m³ | 378.7M m³ | (semana inicial) |

**Consecuencia crítica:** el storage máximo observado (648.84M m³) es **menor** que el `S_min` oficial
(822.18M m³). Es decir, con la capacidad oficial el embalse está **siempre por debajo del umbral
crítico** durante todo el año (el embalse está a ~20% de su capacidad: sequía). Por lo tanto:

- Hoy ambos obtienen **SRS ≈ 0** porque usaron `S_min` preliminar = 162.2M (el storage nunca baja de
  ~349M), lo que vuelve el problema **trivial** (nada que optimizar; baselines = óptimo).
- Con `S_min` oficial = 822.18M, `C_crit ≠ 0` en todas las semanas → los baselines dejan de ser cero y
  **el problema se vuelve no trivial** (hay margen real de optimización reduciendo liberaciones dentro
  del balance permitido).

**Recomendación:** adoptar las constantes oficiales ya resueltas (`S_max`, `S_min`) como verdad común,
y **fijar una convención única para `Δu`** (qué ventana y si la mediana se toma del año completo o por
instancia T). Sin esto, los SRS no son comparables entre integrantes.

Faltantes confirmados (`pdf_vs_data_contrast.csv`, `benchmark_readiness.csv`): `S_max` oficial (se usó
proxy), `Discharge.Total.Change-in-Storage@08461200` (ΔS derivado), evaporación (opcional, ausente).

---

## 5. Formulaciones QUBO

Ambos usan **one-hot** (`x_{t,ℓ}=1 ⟺ u(t)=aₗ`), `u_t=Σ aₗ x_{t,ℓ}`, storage lineal en bits.

- **ivan:** clase `QUBOBuilder` con `add_square_linear(coeffs, const, weight)` (expande
  `(c+Σaᵢxᵢ)²`). Restricciones **exactas con slacks por expansión binaria** (déficit/superávit `d,r`,
  cotas de storage `g,h`, balance `b±`). T=12/L=5 → **278 variables, 7466 términos** (60 de decisión +
  ~218 auxiliares). Funciones reusables: `official_parameters`, `evaluate_srs`, decoder con tolerancia
  one-hot.
- **emilio:** funciones por término `Q_crit / Q_dev / Q_smooth / Q_onehot` (numpy, monolíticas),
  one-hot por penalización con λ=1e10·max(w). T=26/L=5 → **130 variables** (solo decisión, sin slacks;
  la factibilidad se maneja en el solver, no en Q).

Esto mapea directo a la discusión de encoding/constraints: **slacks exactos (ivan)** vs
**post-selección/filtrado (emilio)** para las desigualdades.

---

## 6. Solvers y resultados

- **ivan:** `neal.SimulatedAnnealingSampler` (dimod), T=12, 200 reads, 0.878 s. Resultado SRS −0.0797;
  pero **one-hot imperfecto** (P_one=25 muy débil; alguna semana con suma>1) y **balance violado**
  (|Σu|=71.9M vs permitido 42.7M; P_balance=50 insuficiente). El decoder recupera por argmax.
- **emilio:** brute force T≤3 (valida correctitud del QUBO; T=3 energía −1.108e-02). SA numpy v1
  **infeasible**; **SA v2 con filtrado de factibilidad → feasible** (energía −0.5237, 9.5 s, rechaza
  20.8% de los movimientos). Sin Qiskit/PennyLane/D-Wave-QPU.
- **Nadie corrió QAOA todavía** (oportunidad clara para el track cuántico).

---

## 7. Baselines

Histórico (`u=0`) y regla de umbral (`u=−Δu` si `S<S_min`): ambos dan **SRS ≈ 0** con el `S_min`
preliminar (storage nunca baja de 162M; mínimo simulado ~349M, `weeks_below_Smin=0`, `feasible=True`).
Esto **cambiará** al usar `S_min` oficial (§4).

---

## 8. Issues abiertos (accionable)

a. **Adoptar `S_max`/`S_min` oficiales** (el cambio más importante; vuelve el problema no trivial).
b. **Fijar convención de `Δu`** (ventana + año completo vs por instancia) para que los SRS comparen.
c. Falta el dataset oficial `Discharge.Total.Change-in-Storage@08461200` → ΔS sigue **derivado**
   (etiquetar todo como preliminar hasta conseguirlo).
d. Evaporación ausente (extensión opcional).
e. **Pesos de penalización sin calibrar**: ivan débiles (one-hot/balance se violan); emilio λ=1e10 con
   `w~1e-18` → riesgo numérico. Hace falta barrido de sensibilidad (ver `CLAUDE_Falcon_QUBO_Input.md`
   §14).
f. Rutas Windows hardcodeadas en el notebook de emilio (no portable; usar rutas relativas/env).
g. **QAOA pendiente.**

---

## 9. Activos reutilizables para el builder compartido

- **De ivan:** `QUBOBuilder.add_square_linear` (helper general §12), `official_parameters`,
  `evaluate_srs`, decoder con tolerancia one-hot, estrategia de slacks por expansión binaria.
- **De emilio:** funciones por término `Q_crit/Q_dev/Q_smooth/Q_onehot`, validadores brute-force
  (T≤3), SA con filtrado de factibilidad (post-selección).

Encajan en el builder configurable planificado (código compartido en la raíz de `scripts/`, ver
`docs/GUIDELINES.md`). Conviene unificar: un solo builder por términos componibles + un índice de
variables que abstraiga el encoding + constantes oficiales en un único lugar.
