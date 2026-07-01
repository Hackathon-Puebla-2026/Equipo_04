# Hallazgos clave del equipo (digest para decidir)

Conclusiones decisivas del trabajo de ivan y emilio. Detalle completo: `docs/RESUMEN_HALLAZGOS_EQUIPO.md`.
**Plan activo (qué construir y cómo):** `docs/SPEC_IMPLEMENTACION_QUBO.md`.

- **Datos listos:** 1 año, **53 semanas**, EDA hecha (0 faltantes/0 duplicados). Dataset semanal en
  `notebooks/ivan/Act2/falcon_weekly_benchmark_preliminary.csv`.

- ⚠️ **Constantes (lo más importante):** el equipo usó `S_max = máx observado = 648.84M m³` →
  `S_min = 162.2M m³` (**PRELIMINAR, NO oficial**). La oficial ya resuelta
  (`FalconChallenge/data/falcon_reservoir_constants.json`): **`S_max = 3,288.7M m³`,
  `S_min ≈ 822.18M m³`**. Como el storage máximo observado (648.84M) es **menor** que el `S_min`
  oficial, con la oficial el embalse está **siempre bajo el umbral** → `C_crit ≠ 0` y los baselines
  **dejan de ser triviales**. Hoy ambos dan SRS≈0 solo por usar el `S_min` preliminar bajo.
  **Decisión: usar las constantes oficiales.**

- ⚠️ **`Δu` depende de la ventana** (ivan 2.69M, emilio 2.95M, cálculo previo 2.0M m³).
  **Decisión: fijar una convención única** (qué ventana, año completo vs por instancia) o los SRS no
  comparan entre integrantes.

- **`ΔS_obs` es derivado** (`S(t+1)−S(t)`); falta el dataset oficial
  `Discharge.Total.Change-in-Storage@08461200` → todo resultado es **preliminar**.

- **Encoding:** ambos usan **one-hot**. ivan = restricciones con **slacks exactos** (278 vars en
  T12/L5); emilio = **post-selección/filtrado de factibilidad** (130 vars en T26/L5). **QAOA:
  pendiente por todos** (solo SA clásico hasta ahora).

- **Pesos de penalización sin calibrar** (one-hot/balance se violan en ivan; λ enorme con w~1e-18 en
  emilio). Pendiente barrido de sensibilidad antes de confiar en resultados.

- ⚠️ **`qubo_sa` y `qubo_exhaustive` son solvers CLÁSICOS sobre el QUBO** (numpy: Metropolis /
  enumeración `X@Q`), **NO** simulación cuántica. El campo `n_qubits` (p. ej. 135 en T26/L5) cuenta
  **variables binarias del QUBO** (one-hot 26×5 + slacks), no qubits cuánticos simulados; por eso
  135 "qubits" corren en <1 s (un statevector de 135 qubits, `2^135`, sería imposible). **Aún no se
  corrió ningún algoritmo cuántico ni quantum-inspired** (QAOA = Fase 3, pendiente). Detalle:
  `docs/ANALISIS_DP_Y_RESULTADOS.md`.
