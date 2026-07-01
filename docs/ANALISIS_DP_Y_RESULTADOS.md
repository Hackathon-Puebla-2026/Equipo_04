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

## 6.5. Aclaración: los solvers QUBO de Fase 1 son CLÁSICOS (no simulación cuántica)

> Importante para no malinterpretar los resultados: en `results/runs/` hay corridas `qubo_sa` y
> `qubo_exhaustive` con un campo `n_qubits` (p. ej. **135** en T26/L5). Eso **no** es una simulación
> cuántica. Son solvers **clásicos** (numpy) sobre la matriz `Q` del QUBO.

- **`qubo_sa`** (`scripts/falcon_solvers.py:98`): **simulated annealing clásico** (Metropolis). El
  estado son **26 enteros de nivel** (`lv`), no un statevector; cada paso cambia el nivel de una
  semana y evalúa la energía `H(x)=xᵀQx+const` (`scripts/falcon_qubo.py:42`), con `Q` una matriz
  135×135. Son ~135² flops por paso × ~60 000 pasos ≈ **sub-segundo**. Nada exponencial.
- **`qubo_exhaustive`** (`scripts/falcon_solvers.py:47`): enumeración **vectorizada** de las `L^T`
  combinaciones factibles (`X@Q` por lotes). También clásico; `simulator: "numpy_exhaustive"`.
- **"qubit" aquí = variable binaria del QUBO** (one-hot `T·L` + slacks de balance), no un qubit
  cuántico simulado. Un statevector de 135 qubits ocuparía `2^135·16 B` → físicamente imposible;
  por eso justamente NO se está simulando eso. El annealing clásico es polinomial por paso y esquiva
  la explosión `L^T=5^26≈1.5e18` con un paseo aleatorio sesgado (no la enumera).
- **Estado del repo: cero corridas cuánticas o quantum-inspired.** QAOA es **Fase 3** y sigue
  pendiente (`docs/SPEC_IMPLEMENTACION_QUBO.md` §"Fase 3"; los `*Done:*` de esa sección son
  criterios de terminación, no logros). Los `simulator` registrados son `numpy_sa`,
  `numpy_exhaustive` o `None`.
- **Contraste (repo hermano Georgia):** ahí **sí** se corrió cuántico *simulado* -QAOA en Qiskit Aer
  statevector (≤20 qubits) y MPS tensor-network (≤50 qubits, χ=32), más variantes CVaR/noisy/VQE-. Es
  la maquinaria que Fase 3 debe portar acá; hoy todavía no se ejecutó.

## 7. Pendiente

- **MILP** (opcional): óptimo exacto independiente para verificar el DP en medium (donde brute es
  imposible, `5^26≈1.5e18`). Ver spec `docs/SPEC_IMPLEMENTACION_QUBO.md` (Opcionales).
- **QUBO/QAOA** (Fase 1+): el DP es el ground truth para validar la energía y el óptimo.
- Métodos extra de ivan: regla de umbral balanceada **ya incorporada** (`threshold/balanced`);
  annealing quantum-inspired pendiente.

## 8. Escalamiento de qubits del QUBO y restricciones eliminables

> Snapshot 2026-06-30. Análisis para Fase 1+ (aún no construimos el QUBO). Motiva encodings compactos
> y el experimento de chunking (E2, ver spec Fase 2/3).

### 8.1 Qubits por encoding (bits de decisión)

| Instancia | one-hot `T·L` | binary `T·⌈log₂L⌉` | domain-wall `T·(L-1)` |
|---|---:|---:|---:|
| small T12/L3 | 36 | 24 | 24 |
| medium T26/L5 | 130 | 78 | 104 |
| large T52/L5 | 260 | 156 | 208 |
| large T52/L7 | 364 | 156 | 312 |

Eso es solo decisión. Los **slacks de restricciones exactas** son lo que dispara el conteo (ivan
T12/L5 = 278 vars, ~218 slacks). Contra el techo de statevector (~30 qubits, T4 16 GB), incluso el
medium lean (130) queda lejos -> hacen falta encodings compactos y/o chunking.

### 8.2 ¿Qué restricciones podemos quitar? (medido sobre nuestros datos)

Por instancia, ¿cada restricción llega a atar?

| Restricción (spec) | ¿Ata en nuestros datos? | Costo en qubits | Decisión |
|---|---|---|---|
| `0 ≤ S ≤ S_max` (ec. 10) | **NUNCA** (peor caso p.ej. `[55M,1125M]` en T52/L7 vs `[0,3289M]`) | slacks por semana (grande) | **QUITAR** |
| `\|Σu\| ≤ B` (ec. 11) | **Sí, al límite** (DP usa ~99% de B en medium/large) | 0 (soft) / pocos bits (exacta) | mantener |
| `R(t) ≥ 0` (ec. 8) | Parcial (1/8/20 niveles prohibidos en medium/large) | 0 (penalización lineal / nivel prohibido) | mantener (barato) |
| `\|u\| ≤ u_max` (ec. 9) | Automática (niveles ⊂ `[-u_max,u_max]`) | 0 | por construcción |
| one-hot (selección única) | Estructural | 0 (penalización) / mixer | mantener |

**Conclusión:** con post-selección + términos soft + **quitar las cotas de storage**, el QUBO queda en
**`T·L` qubits, cero slacks**. Solo el balance exacto o `C_crit` Opción B (déficit/superávit) agregan
unos pocos. Es decir, para ESTOS datos el temido "las restricciones explotan los qubits" es evitable.
Aun así, `T·L` en medium (130) supera el statevector -> ver encodings compactos (§8.1) y chunking (E2).

### 8.3 Balance exacto: casi gratis (preferirlo sobre el soft)

El balance es UNA restricción global, no por-semana. En unidades enteras `Σu = M·Δu`, `M=Σkₜ` entero,
y `|Σu|≤B ⟺ |M| ≤ M_cap=⌊B/Δu⌋`. Con **un solo slack acotado** (`M+s=M_cap`, `s∈[0,2·M_cap]`,
log-encoded) el costo es:

| Instancia | `M_cap` | slack exacto (1 var) |
|---|---:|---:|
| debug T5/L3 | 1 | 2 qubits |
| small T12/L3 | 4 | 4 qubits |
| medium T26/L5 | 11 | 5 qubits |
| large T52/L5 | 22 | 6 qubits |
| large T52/L7 | 22 | 6 qubits |

Crece **logarítmico** en T (`M_cap≈0.4·T`, bits `≈log₂(0.8·T)`): 6 qubits en el peor caso, `<3%` sobre
`T·L`. **Recomendación: usar el balance EXACTO de 1 slack** (`balance="slack"`), no el soft
`P_bal(Σu)²`. El soft cuesta 0 qubits pero es aproximado (sesga `Σu→0`, no impone `|Σu|≤B`); solo
conviene como atajo de MVP. El balance exacto es correcto y prácticamente gratis.

## 9. QAOA (Fase 3) - primer cuántico simulado (2026-07-01)

Primer algoritmo **cuántico simulado** del repo (contrasta con §6.5: SA/exhaustive eran clásicos).
`scripts/falcon_qaoa.py` + driver `falcon_run_qaoa.py`, con **Qiskit + Aer statevector** (correr con
`.venv-quantum/bin/python`, Python 3.12). Circuito QAOA manual (`|+>ⁿ`; por capa: costo RZ/RZZ desde
`qp.to_ising()`, mixer RX), COBYLA sobre `<H>`, seed=42 con ≥5 restarts, transpile 1 vez + bind de
params. `<H>` y muestreo se calculan desde el statevector + diagonal QUBO precomputada (mismo indexado,
sin reversos de endianness). Solo statevector (n≲26 qubits); instancias grandes → MPS/sampling.

**Validaciones (todas pasan):** diagonal precomputada == `qubo_energy` (exacto); `H_ising` del circuito
== `diag−const` a `1e-14` (el circuito optimiza exactamente el landscape que mido); gate energía
`·scale == −SRS` a `~1e-15`; y para el **binary encoding** el mínimo global sobre `2ⁿ` es válido +
factible y `== DP*` (0 estados por debajo del óptimo → los codewords inválidos quedan bien penalizados).

**Binary encoding (Fase 2):** `l∈{0..L-1}` en `⌈log₂L⌉` bits/semana (`u` lineal en bits con niveles
equiespaciados). Necesario para que **T12/L3 entre en statevector: 24 qubits binary vs 36 one-hot**.
Exacto en QUBO para `b≤2` (L≤4); codewords inválidos (`l≥L`) se penalizan con producto de bits. Cambiar
de encoding toca solo la capa `falcon_encodings.py` (los términos de costo no cambian).

**Resultados (ΔSRS = SRS_qaoa − SRS_baseline; `hist=DP*` en estas ventanas: u=0 ya es óptimo):**

| Instancia | encoding | p | n_qubits | SRS_qaoa | DP\* | AR | factible | runtime |
|---|---|---:|---:|---:|---:|---:|---|---:|
| debug T5/L3 | one-hot | 1 | 17 | −0.276985 | −0.276985 | **1.000** | sí | 2.5 s |
| debug T5/L3 | one-hot | 2 | 17 | −0.276985 | −0.276985 | **1.000** | sí | 9.5 s |
| small T12/L3 | binary | 1 | 24 | −0.389293 | −0.296264 | 1.314 | sí | 189 s |

- **debug (17 q): QAOA alcanza el óptimo** (AR 1.000) a p=1. (Requiere muestrear top-256 por
  probabilidad al decodificar: el óptimo tiene solo ~0.2% de probabilidad, un top-64 lo perdía.)
- **small (24 q): QAOA factible pero subóptimo** (AR 1.31, prob del óptimo ~5e-6): el límite de
  profundidad p=1 se nota al escalar. Mejorarlo: mayor p con init INTERP (naive p=2 no ayuda),
  o **XY-mixer** (restringe al subespacio factible y elimina la penalización one-hot).
- **Calibración de penalties (spec §8, `falcon_calibrate_penalties.py`):** barrer el multiplicador de
  penalties `∈{0.1..10}` **no mueve** el resultado en debug. Motivo: `normalize="maxabs"` reescala `Q`
  por su máximo (dominado por las penalties), así que multiplicarlas se cancela al normalizar. El
  cuello no es el condicionamiento por escala de penalty sino la **profundidad/mixer**.
- Consistente con el spec §7: no se busca ventaja cuántica; el objetivo es codificar/benchmarkear.
  Sigue todo `preliminary` (falta el `ΔS_obs` oficial). Registrado en `results/runs_summary.csv`
  (`method="qaoa"`, con `p_depth/beta/gamma/optimizer_iterations`).

