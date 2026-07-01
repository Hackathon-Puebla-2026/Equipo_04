# Comparación de resultados: pipeline `scripts/` vs ivan

Comparación de nuestros baselines/DP (`scripts/falcon_*.py`, corrida en `results/runs_summary.csv`)
contra los resultados reportados por ivan en `docs/FalconChallenge_Report_para_Claude.md` (reporte
LaTeX->md, autoritativo) y su notebook `notebooks/ivan/Act2/FalconChallenge_Colab_Solution(1)(1).ipynb`.
Ambos usan las constantes oficiales (`S_max=3,288.726 Mm3`, `S_min=822.182 Mm3`, `eta=0.10`).

## Resumen del reporte de ivan

- Benchmark oficial **T=26, L=5** (ventana 2025-06-29 a 2025-12-28): `SRS_DP=-0.290423`,
  `SRS_hist=-0.311534`, `ΔSRS=+0.021111`. Política óptima: `-Δu` durante 11 semanas y luego 0
  (`sum_u=-30.848 Mm3`, dentro de `B=31.151`). `Δu=2.804 Mm3`.
- Métodos: histórico (`u=0`), regla de umbral (pura), **regla de umbral balanceada** (frena las
  reducciones cuando violarían el balance; aquí iguala al DP), **annealing quantum-inspired**
  (`-0.291188`) y **DP exacto**.
- Tabla de escalamiento (su Tabla 4): T12/L3, T26/L5, T52/L5, T52/L7.

## Tabla comparativa (SRS del óptimo DP)

| Instancia | ivan `SRS_DP` | nuestro `SRS_DP` | ivan `SRS_hist` | nuestro `SRS_hist` | ΔSRS (ambos) |
|---|---:|---:|---:|---:|---:|
| small T12/L3 | -0.296264 | **-0.296264** | -0.296264 | -0.296264 | 0.000000 |
| medium T26/L5 | -0.290423 | **-0.290423** | -0.311534 | -0.311534 | +0.021111 |
| large T52/L5 | -0.201466 | **-0.201466** | -0.252702 | -0.252702 | +0.051236 |
| large T52/L7 | -0.177643 | **-0.177643** | -0.252702 | -0.252702 | +0.075059 |

**Coinciden exactamente en las 4 instancias** (a 6 decimales). Validación mutua de todo el pipeline
(datos semanales, constantes, SRS, DP).

## Bug encontrado al comparar (en NUESTRO código) y su corrección

Inicialmente diferíamos en T12/L3 (nosotros `-0.293529` con `ΔSRS +0.0027`; ivan `-0.296264`,
`ΔSRS 0`) y en T52/L7 (nosotros dábamos lo mismo que L5, `-0.201466`; ivan `-0.177643`).

**Causa:** `instance_params` fijaba `u_max = 2·Δu` para TODO `L`. Eso solo es correcto para `L=5`
(oficial). La definición consistente es `u_max = half·Δu = max|niveles|`:
- `L=3` -> `u_max = Δu` (no `2Δu`)
- `L=5` -> `u_max = 2Δu` (igual, sin cambio)
- `L=7` -> `u_max = 3Δu` (no `2Δu`)

Como los pesos `w2 = 0.1/(T·u_max²)` y `w3 = 0.1/((T-1)(2u_max)²)` dependen de `u_max`, con el valor
equivocado en L=3 la desviación quedaba artificialmente barata y el DP "mejoraba" el histórico de
forma espuria; en L=7 el `u_max` chico impedía usar los niveles ±3, dando lo mismo que L=5.

**Fix:** `u_max = ((L-1)//2)·Δu` en `scripts/falcon_constants.py::instance_params`. Tras el fix,
nuestro DP coincide con ivan en todas las instancias (tabla de arriba). L=5 no cambió.

**Lección sobre el brute force:** antes del fix, el brute force "confirmaba" nuestro `-0.293529`
porque usaba NUESTROS pesos (los equivocados). Es decir, el brute force valida al optimizador **dado
el objetivo**, pero no puede detectar que el objetivo (los pesos vía `u_max`) estaba mal. La
comparación con ivan fue lo que expuso el bug. Con el `u_max` corregido, `brute == dp` sigue valiendo
(en debug T5/L3 y small T12/L3, `brute == dp == ivan`).

## Métodos de ivan que aún no tenemos

- **Regla de umbral balanceada** (detiene reducciones al llegar al límite de balance; factible, iguala
  al DP en T26). Similar a nuestra idea de post-selección/factibilidad; candidata a agregar.
- **Annealing quantum-inspired** (`-0.291188` en T26): heurística tipo QUBO/Ising, no supera al DP en
  esta ventana. Referencia para nuestra Fase 3 (QAOA).

## Conclusión

Tras corregir el `u_max`, **nuestro pipeline y el de ivan coinciden exactamente en small/medium/large
(L5 y L7)**. El DP exacto es ground truth compartido. La discrepancia inicial era un bug nuestro, no de
ivan. Todo sigue marcado `preliminary` (ΔS derivado, falta el dataset oficial de cambio de
almacenamiento).
