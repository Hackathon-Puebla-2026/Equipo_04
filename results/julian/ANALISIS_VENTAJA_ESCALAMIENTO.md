# Escalamiento, ventaja eventual y feasibilidad del enfoque cuántico

> Figuras en `results/julian/figures/` (generadas por `scripts/julian/falcon_advantage_figures.py`).
> Marco **honesto** (spec §7): el objetivo NO es demostrar ventaja cuántica. Para el problema oficial
> (**un solo embalse**) el DP clásico es exacto y **polinomial** → gana. La ventaja "eventual" vive en
> una extensión **multi-embalse** (proyección, fuera del benchmark oficial).

## TL;DR

- **Hoy, 1 embalse:** el DP resuelve el óptimo exacto en sub-segundo hasta T=52. No hay ventaja
  cuántica y no la reclamamos.
- **El valor cuántico entregado:** una implementación **factible, mitigada y benchmarkeada** en QPU
  real; y una arquitectura (**chunking**) que corre en hardware de pocos qubits para **cualquier**
  horizonte T.
- **Cuándo habría ventaja:** cuando el problema pierde la subestructura que hace barato al DP -
  típicamente **R embalses acoplados**: el estado del DP crece **exponencial en R**, los qubits
  crecen **lineal** (R·T·L). Ahí una (FT-)QPU maneja el QUBO full/por-batches que el DP ya no puede.

## Tabla-resumen de escalamiento (por enfoque)

| Enfoque | Qubits | Tiempo | Memoria | Exacto | Corre en HW hoy |
|---|---|---|---|---|---|
| Fuerza bruta | — | `O(Lᵀ)` | `O(1)` | sí | — |
| **DP exacto** (1 embalse) | — | `O(T²·L²)` | `O(T²·L)` | **sí** | — (clásico) |
| QAOA sim **full** | `T·L` | `O(2^{T·L}·it)` | `2^{T·L}·16B` | no (heurístico) | no (>30 qubits no simula) |
| QAOA sim **chunked** | `b·L` (const) | `O((T/b)·2^{b·L}·it)` (**lineal en T**) | `2^{b·L}·16B` (**const**) | no | sí (bloques chicos) |
| QAOA **hardware** | `T·L` o `b·L` | `O(depth·shots)` (poli en qubits) | — (QPU) | no | **sí** |
| DP conjunto (R embalses) | — | `O(T²·L²·((4T+1)L)^{R−1})` | `((4T+1)L)^R` (**exp. en R**) | sí | — |
| QUBO/QAOA (R embalses) | `R·T·L` (lineal) | poli en qubits (HW) | por-batch acotado | no | sí (FT para full) |

`T`=horizonte (semanas), `L`=niveles (5 oficial), `b`=weeks/chunk (5 → 25 qubits/chunk), `R`=nº embalses.

## Figuras

**figA_time_scaling.png** - Tiempo vs T (log). Brute `Lᵀ` y QAOA-sim-full `2^{T·L}` explotan; el **DP es
polinomial** (marcadores medidos, docs §3); **QAOA-chunked es lineal en T** y QAOA-hardware es
`depth·shots` (sin el factor `2ⁿ` de la simulación). Mensaje: chunking y hardware esquivan el muro
exponencial de simular el estado cuántico; el DP es el ganador clásico a 1 embalse.

**figB_memory_scaling.png** - Memoria vs T (log). El statevector full `2^{T·L}·16B` cruza los 16 GB
alrededor de ~30 qubits (T≈6 en one-hot L=5); el **chunked queda plano** (`2^{b·L}·16B`); el DP es
minúsculo. Es la razón física por la que la simulación clásica exacta se corta ~30 qubits y por la que
el hardware es la única vía para circuitos grandes.

**fig1_qubits_vs_T.png** - Qubits vs T por encoding (one-hot/binary/domain-wall) con el muro
statevector ~30. El full-run crece **lineal con T**; el **chunking fija el circuito en b·L=25** →
independiente de T. Un horizonte largo full-run necesita FT-QPU; el chunked corre "en batches" hoy.

**fig2_dp_vs_brute.png** - Un embalse: candidatos `Lᵀ` (exponencial) vs **estados DP ~T²·L**
(polinomial, medido: 433/6505/25086). El DP colapsa la explosión → **no hay ventaja cuántica** en el
problema oficial. Baseline honesto.

**fig3_multireservoir_crossover.png** *(proyección)* - Ventaja eventual. Con R embalses acoplados el
estado conjunto del DP `~((4T+1)L)^R` es **exponencial en R**; los qubits `R·T·L` son **lineales**. Hay
un **crossover R\*** donde el DP supera el muro clásico tratable (~1e12 estados) mientras los qubits
siguen en rango near-term/FT → régimen donde la QPU corre lo que el DP no puede. **Extensión fuera del
benchmark oficial**, con los supuestos escritos en la figura.

**figC_depth_vs_qubits.png** - Recurso en hardware: profundidad transpilada **medida** (FakeBrisbane,
opt=3) de nuestros bloques QAOA-XY. **Medido:** 10q→674, 15q→853, 20q→1257, **25q→1563** vs Ivan
~1900 @ 28q. Lectura honesta: **el driver de profundidad es el tamaño del bloque**, no solo el mixer -
bloques chicos (T12/L3, 12q) son **~3× menos profundos**; a 25q el XY iguala el orden de Ivan (el
`initialize` del W-state + RXX/RYY del XY-mixer no es gratis a escala). Conclusión: para minimizar ruido
conviene **bloques más chicos**, y el XY compra factibilidad (subespacio one-hot) a cambio de algo de
profundidad. Trade-off explícito.

**fig4_feasibility.png** - Datos **reales** de `results/runs_summary.csv` (solo ventana canónica).
SRS por instancia para histórico/DP/QAOA-sim/QAOA-chunked/QAOA-hardware; hatch ✗ = infactible. El
enfoque cuántico con **chunks + relajaciones** (balance soft, XY-mixer, `eta_local`, post-selección)
da soluciones **factibles y cerca del DP**; el hardware (T5/L3, AR=1.0) corre honesto.

## Resultados en HARDWARE real (chunked-QAOA, ibm_kingston, 2026-07-01)

Corridas recuperadas por `job_id` (`scripts/julian/falcon_hw_retrieve.py`); análisis en
`results/julian/hw_analysis_*.json`, figuras `figHW1/2/3`. Registrado `method="qaoa_chunked_hardware"`.

**small T12/L3** (3 bloques de 4 semanas = 12 qubits/bloque, XY-mixer, DD + gate/measure twirling,
4096 shots): **SRS = −0.296264 = óptimo DP, FACTIBLE, AR = 1.000.**

| bloque | one-hot válido en HW | muestras factibles | mejor SRS bloque |
|---|---:|---:|---:|
| 0 | 31.6% | 19 | −0.107 |
| 1 | 31.2% | 71 | −0.110 |
| 2 | 35.5% | 51 | −0.123 |

- **figHW1** - el ruido de hardware tira **~⅔ de los shots fuera del subespacio one-hot** (solo ~32%
  válido vs 100% en sim), pero el **XY-mixer mantiene esa fracción usable** y la **post-selección** se
  queda con la mejor factible → bloque factible fiable.
- **figHW2** - distribución de SRS de las muestras factibles por bloque (dispersa por ruido); se marca
  la que eligió la post-selección.
- **figHW3** - SRS HW-chunked vs DP/histórico: en small **el HW alcanza el óptimo** (en régimen de
  sequía el óptimo es `u=0`, y aun con ruido se recupera factible). **Supera la corrida previa del
  equipo** (T26/L5, −0.3322, **infactible**, sin mitigación, params heurísticos).
- medium T26/L5 y large T52/L5: enviadas vía `falcon_hw_submit.py` (manifest
  `results/julian/hw_jobs_manifest.json`); se recuperan con `falcon_hw_retrieve.py` al quedar DONE.

**Lectura honesta:** no es speedup ni mejor SRS que el DP; es una **ejecución cuántica real,
factible y mitigada**, con el análisis del impacto de ruido medido en la QPU.

## Por qué el chunking hace la diferencia (no el DP)

El chunking ataca la **simulabilidad/ejecutabilidad del circuito**, no la complejidad del DP:
- **Circuito acotado:** `b·L` qubits fijos → entra en statevector y en QPU de pocos qubits para
  cualquier T. Sin chunking, T26/L5 = 130 qubits one-hot (imposible de simular; profundo en HW).
- **Menos profundidad → menos ruido:** la profundidad la fija el tamaño del bloque (12q≈674, 25q≈1563
  transpilado); bloques más chicos = circuitos más superficiales (ver figC).
- **Costo:** rompe la optimalidad global (el balance/redistribución que cruza fronteras se pierde;
  ver `docs/ANALISIS_DP_Y_RESULTADOS.md` §10, gap medido ~0.021 en T26/L5). Trade-off explícito:
  factibilidad/ejecutabilidad en HW a cambio de un gap de optimalidad acotado.

## Conclusión honesta

1. **Problema oficial (1 embalse):** clásico gana; el DP es exacto y polinomial. Reportamos el
   cuántico como implementación correcta + benchmark, no como speedup.
2. **La QPU es *necesaria*** cuando el circuito supera lo simulable (~30 qubits) - hoy vía chunking en
   batches; a futuro full-run en FT-QPU.
3. **La ventaja *aparecería*** en la extensión multi-embalse (o objetivos no separables) donde el DP se
   vuelve exponencial y los qubits siguen lineales (fig 3, proyección).
4. **Lo demostrado ya:** el enfoque cuántico con chunks y relajaciones es **factible** en hardware real,
   con mitigación de error, mejorando la corrida previa infactible del equipo.
