# Análisis del DP (complejidad y tiempos) y resultados actuales

> **Snapshot: 2026-06-30.** Números y tiempos medidos en esta fecha (corrida en
> `results/runs_summary.csv`, timestamps `20260630-*`). Se desactualizarán al cambiar datos, código o
> convenciones; re-medir y actualizar la fecha si se regeneran los resultados.

Consolida el estado tras Fase 0. Datos crudos: `results/runs_summary.csv` (+ JSON por corrida en
`results/runs/`). Comparación con ivan: `results/comparacion_vs_ivan.md`. Todo en m^3, constantes
oficiales, con `ΔS_obs` derivado (resultados marcados `preliminary`).

## 1. Estado del pipeline (Fase 0 completa)

Módulos compartidos en `scripts/`:

- `falcon_constants.py`: constantes oficiales (`S_max=3,288,726,000 m^3`, `S_min=822,181,500 m^3`),
  `instance_params` (Δu por instancia, `u_max=half·Δu`, niveles), pesos oficiales `w1,w2,w3`.
- `falcon_data.py`: reconstrucción canónica del dataset semanal (52 semanas, m^3) desde
  `FalconChallenge/data/`; `S_obs` cierre de semana, `R_obs` integración nativa, `ΔS` derivado.
- `falcon_storage.py`: dinámica `S(t+1)=S(t)+ΔS(t)-u(t)`, release y chequeo de restricciones.
- `falcon_srs.py`: costos `Ccrit/Cdev/Csmooth` y `SRS`.
- `falcon_baselines.py`: `historical`, `threshold_rule` (pura + clamped), `dp_optimal` (exacto),
  `brute_force_optimal` (árbitro, instancias chicas).
- `falcon_results.py`: writer estandarizado (`results/runs_summary.csv` + `results/runs/*.json`).
- `falcon_run_baselines.py`: corre y registra todos los baselines en las instancias estándar.

## 2. Resultados actuales

Referencia primaria del ΔSRS = `SRS_hist` (histórico `u=0`). Números de `results/runs_summary.csv`:

| Instancia | T | L | SRS_hist | threshold (pura) | SRS_DP | SRS_brute | ΔSRS_DP vs hist | DP factible |
|---|---:|---:|---:|---:|---:|---:|---:|:--:|
| debug | 5 | 3 | -0.276985 | -0.370566 (infactible) | -0.276985 | -0.276985 | 0.000000 | Sí |
| small | 12 | 3 | -0.296264 | -0.374672 (infactible) | -0.296264 | -0.296264 | 0.000000 | Sí |
| medium | 26 | 5 | -0.311534 | -0.288934 (infactible) | -0.290423 | n/a | +0.021111 | Sí |
| large | 52 | 5 | -0.252702 | -0.202342 (infactible) | -0.201466 | n/a | +0.051236 | Sí |
| large | 52 | 7 | -0.252702 | -0.188453 (infactible) | -0.177643 | n/a | +0.075059 | Sí |

Notas:
- La **regla de umbral pura es infactible** por balance en todas las instancias (reduce más de lo
  permitido por `|Σu| ≤ η·ΣR_obs`); su SRS a veces parece mejor solo por violar esa restricción. El DP
  es el mejor **factible**. La variante clamped solo cambia en large (evita `R<0`), sin volverse factible.
- La **regla de umbral balanceada** (variante de ivan, `threshold/balanced`) detiene la reducción al
  llegar al presupuesto de balance -> es **factible en todas las instancias** y `≤ DP`. En **medium
  coincide con el DP** (`-0.290423`). SRS balanceada por instancia: debug -0.301065, small -0.320233,
  medium -0.290423, large L5 -0.206566, large L7 -0.200622.
- En debug/small el óptimo factible es `u=0` (ΔSRS 0): con `u_max=Δu` en L=3 la desviación es cara y
  no compensa; a mayor horizonte (medium/large) el DP sí mejora al histórico.
- `brute_force == dp` donde el brute es enumerable (debug, small), confirmando el óptimo.

## 3. Complejidad del DP: por qué es rápido (y correcto)

Estado del DP: `(t, C_t, k_prev)`, con `C_t = Σ_{j<t} k_j` (suma entera acumulada de niveles) y
`k_prev` para `C_smooth`. Como `S(t)=H_t - Δu·C_t`, el almacenamiento depende SOLO de `C_t`
(no del camino). `C_t ∈ [-half·t, half·t]`, así que el nº de estados es **polinomial**: `O(T²·L²)`.

**Estados realmente explorados vs espacio combinatorio `L^T`** (medido):

| Instancia | Espacio `L^T` | Estados DP | Transiciones DP | Runtime DP (s) |
|---|---:|---:|---:|---:|
| T12/L3 | 5.31e5 | 433 | 1,092 | 0.0009 |
| T26/L5 | **1.49e18** | **6,505** | **30,024** | 0.020 |
| T52/L5 | 2.22e36 | 25,086 | 116,676 | 0.076 |
| T52/L7 | 8.81e43 | 49,739 | 317,615 | 0.195 |

El DP para medium hace ~30 mil operaciones, no `1.49e18`. Por eso corre en ~0.02 s. **El runtime
sub-segundo es lo esperado**, no un error: es el colapso exponencial->polinomial por subestructura
óptima (muchísimos calendarios distintos llegan al mismo estado `(C_t, k_prev)`; el DP guarda solo el
de menor costo). Si tardara horas, ESO indicaría que estaríamos enumerando caminos por error.

## 4. La poda es sin pérdida (no descarta el óptimo)

Dos mecanismos:
1. **Fusión de estados** (el grueso del ahorro): válido por propiedad de Markov / subestructura óptima;
   el costo futuro depende solo de `(C_t, k_prev)`. `k_prev` está en el estado justo para `C_smooth`.
2. **Poda de infactibles duros**: descarta estados con `S ∉ [0,S_max]`, niveles con `R<0`, y balance
   terminal `|Σu|>B`. Esos estados nunca pueden ser parte de un óptimo factible, así que podarlos no
   pierde nada.

**Prueba de que es sin pérdida:** `brute_force == dp` en T5/L3 y T12/L3 (acuerdo exhaustivo). Un
algoritmo exacto validado contra fuerza bruta en casos chicos es confiable a escala.

## 5. Validación cruzada y un bug encontrado

- **Match exacto con ivan** (implementación independiente) en las 4 instancias
  (small/medium/large L5 y L7). Ver `results/comparacion_vs_ivan.md`.
- **Bug propio corregido**: `instance_params` fijaba `u_max=2·Δu` para todo `L`; correcto solo en L=5.
  Se corrigió a `u_max=half·Δu` (Δu en L3, 2Δu en L5, 3Δu en L7). La comparación con ivan lo destapó;
  el brute "confirmaba" el valor viejo porque usaba los pesos equivocados (valida al optimizador dado
  el objetivo, no al objetivo). Tras el fix, todo coincide.

## 6. Datos y caveat

Entradas en m^3 desde `data/processed/falcon_weekly_benchmark.csv`; constantes oficiales de
`FalconChallenge/data/falcon_reservoir_constants.json`. `ΔS_obs` es **derivado** (`S(t+1)-S(t)`),
falta el dataset oficial `Discharge.Total.Change-in-Storage@08461200`, por eso todo se marca
`preliminary`.

## 7. Pendiente

- **MILP** (opcional): óptimo exacto independiente para verificar el DP en medium (donde brute es
  imposible, `5^26≈1.5e18`). Ver spec `docs/SPEC_IMPLEMENTACION_QUBO.md` (Opcionales).
- **QUBO/QAOA** (Fase 1+): el DP es el ground truth para validar la energía y el óptimo.
- Métodos extra de ivan: regla de umbral balanceada **ya incorporada** (`threshold/balanced`);
  annealing quantum-inspired pendiente.
