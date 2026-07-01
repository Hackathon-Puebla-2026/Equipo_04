# Falcon Challenge A — Análisis integral del notebook, resultados e IBM Quantum

**Propósito del documento:** entregar a Claude un contexto técnico completo, autocontenido e integrado con el análisis previo del proyecto. Este markdown resume qué se hizo, qué resultados se obtuvieron, cómo interpretar los resultados de IBM Quantum y qué ajustes conviene pedir en la siguiente iteración.

**Archivos base analizados:**

- `FalconChallenge_IBM_Qiskit_Runtime_30qubits(3).ipynb` — notebook ejecutado con resultados.
- `job-d92j4vd958jc73brsfg0-info.json` — metadatos del job IBM.
- `job-d92j4vd958jc73brsfg0-result.json` — salida cruda del Sampler de IBM.
- `FalconChallenge_analisis_output_IBM_para_Claude.md` — análisis previo del notebook.
- `FalconChallenge_analisis_job_IBM_complementario.md` — análisis previo del job IBM.
- Gráficas exportadas:
  - `falcon_ibm_markdown_assets/grafica_almacenamiento_observado_dp_ibm.png`
  - `falcon_ibm_markdown_assets/grafica_liberacion_observada_dp_ibm.png`

---

## 0. Instrucción sugerida para Claude

Claude debe usar este documento para integrar el resultado IBM dentro del análisis global del proyecto Falcon Challenge. La conclusión principal no debe redactarse como “ventaja cuántica” ni como “IBM superó a DP”. La conclusión correcta es:

> Se logró cerrar el ciclo metodológico completo: datos hidrológicos → benchmark clásico → DP exacta → QUBO → codificación domain-wall → QAOA por bloques ≤30 qubits → ejecución en IBM Quantum → reconstrucción de política semanal → evaluación con SRS. La ejecución cuántica fue funcional y respetó el límite de qubits, pero el resultado QPU actual no superó a DP exacta y no fue factible por balance acumulado.

Claude debe mantener la DP exacta como benchmark principal factible, y presentar IBM/QAOA como demostración de escalamiento cuántico/híbrido e integración experimental.

---

## 1. Resumen ejecutivo

El notebook ejecutado confirma que el problema del **Falcon Challenge A** fue escalado desde una formulación clásica hacia una formulación cuántica/híbrida ejecutable en IBM Quantum. El horizonte oficial de `T=26` semanas y `L=5` niveles de ajuste fue dividido en cuatro bloques `7 + 7 + 7 + 5`, usando codificación **domain-wall**.

Con `L=5`, cada semana requiere `L-1 = 4` qubits. Por tanto:

- Bloques de 7 semanas: `7 × 4 = 28 qubits`.
- Bloque final de 5 semanas: `5 × 4 = 20 qubits`.
- Máximo usado: `28 qubits`, por debajo del límite solicitado de `30 qubits`.

El job IBM se ejecutó correctamente en `ibm_fez`, backend de 156 qubits, con `2048 shots` por circuito. El job reportó estado `Completed`.

El resultado IBM/QAOA reconstruido obtuvo:

```text
SRS IBM/QAOA = -0.33216410584163464
SRS DP exacta = -0.293160
Factibilidad IBM/QAOA = False
Restricción fallida = balance acumulado
```

La interpretación global es:

- **DP exacta** sigue siendo la mejor solución factible y el benchmark técnico central.
- **IBM/QAOA** sí demuestra ejecución real y escalamiento a hardware cuántico bajo 30 qubits.
- **IBM/QAOA no es todavía competitivo** por ruido, profundidad de circuitos, ausencia de optimización de ángulos QAOA y debilidad de la restricción global de balance al resolver por bloques.

---

## 2. Contexto matemático del challenge

La variable de decisión es la secuencia de ajustes de liberación:

```text
u(t)
```

La liberación optimizada se define como:

```text
R(t) = R_obs(t) + u(t)
```

La dinámica simplificada de almacenamiento es:

```text
S_opt(t+1) = S_opt(t) + ΔS_obs(t) - u(t)
```

La métrica principal es el **Storage Resilience Score (SRS)**:

```text
SRS = - (w1 Ccrit + w2 Cdev + w3 Csmooth)
```

Donde:

```text
Ccrit   = Σ_t [max(0, Smin - Sopt(t))]^2
Cdev    = Σ_t u(t)^2
Csmooth = Σ_t [u(t) - u(t-1)]^2
```

Dado que el SRS está definido como negativo de penalizaciones, **un valor más cercano a cero es mejor**.

Restricciones principales:

```text
R(t) ≥ 0
|u(t)| ≤ umax
0 ≤ Sopt(t) ≤ Smax
|Σ_t u(t)| ≤ η Σ_t R_obs(t)
```

La restricción que falló en IBM/QAOA fue la última: **balance acumulado**. Esto significa que la política QPU redistribuyó o redujo liberaciones de forma globalmente incompatible con la restricción del benchmark.

---

## 3. Datos y configuración del notebook

El notebook cargó los datos IBWC y constantes del embalse Falcon. La ventana semanal usada fue:

| Elemento | Valor |
|---|---:|
| Semanas `T` | 26 |
| Fecha inicial | 2025-06-29 |
| Fecha final | 2025-12-28 |
| `Smax` | 3,288,726,000 m³ |
| `Smin` | 822,181,500 m³ |
| `η` | 0.10 |
| `L` | 5 niveles |

Los niveles discretos fueron:

```text
k ∈ {-2, -1, 0, 1, 2}
u(t) = k(t) Δu
```

El notebook calculó:

```text
Δu = 2,706,476 m³/semana
u niveles = [-5,412,952.42, -2,706,476.21, 0, 2,706,476.21, 5,412,952.42]
```

Escala combinatoria del problema oficial:

```text
Número de schedules = L^T = 5^26 = 1,490,116,119,384,765,625
```

Por eso la solución directa por enumeración no es práctica. La DP exacta reduce el costo al trabajar con estados y transiciones; QUBO/QAOA transforma el problema en optimización binaria por bloques.

---

## 4. Benchmark clásico del notebook

El notebook evaluó tres soluciones clásicas iniciales:

1. Histórico: `u(t)=0`.
2. Regla umbral: reduce liberación cuando el almacenamiento cae por debajo de `Smin`.
3. DP exacta: optimización discreta sobre los niveles permitidos.

| modelo | SRS | Ccrit | Cdev | Csmooth | sum_u_m3 | min_storage_pct_cap | factible | runtime_s |
|---|---:|---:|---:|---:|---:|---:|---|---:|
| Histórico | -0.311534 | 5.685972e18 | 0 | 0 | 0 | 9.653041 | True | 0.000000 |
| Regla umbral | -0.290505 | 4.845873e18 | 1.904504e14 | 0 | -7.036838e7 | 11.298953 | False | 0.000000 |
| DP exacta | -0.293160 | 5.148445e18 | 8.057515e13 | 1.465003e13 | -2.977124e7 | 10.558293 | True | 0.049052 |

### Interpretación

La regla umbral obtiene el SRS más cercano a cero, pero **no es factible**. En la comparación final se observa que genera una liberación mínima negativa:

```text
min_release_m3 = -1,847,126
```

Por eso no debe ser tratada como solución ganadora, sino como baseline heurístico que mejora almacenamiento a costa de violar restricciones.

La DP exacta es la mejor solución factible. La secuencia DP reportada fue:

```text
k_DP = [0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
```

Lectura operativa de DP:

- Reduce liberación moderadamente durante las primeras semanas críticas.
- Evita liberaciones negativas.
- Respeta balance acumulado.
- Mantiene la trayectoria más suave que IBM/QAOA.

Comparado con histórico, DP mejora el SRS en aproximadamente:

```text
ΔSRS_DP_vs_histórico = 0.018374
reducción aproximada de penalización = 5.90%
reducción de Ccrit ≈ 9.45%
```

---

## 5. Escalamiento cuántico: partición por bloques

Para cumplir el límite de 30 qubits, el notebook usó:

```text
Qubits por semana = L - 1 = 4
Máximo semanas por bloque = floor(30 / 4) = 7
```

Partición usada:

| block_id | start_week | end_week_exclusive | weeks | qubits_domain_wall |
|---:|---:|---:|---:|---:|
| 0 | 0 | 7 | 7 | 28 |
| 1 | 7 | 14 | 7 | 28 |
| 2 | 14 | 21 | 7 | 28 |
| 3 | 21 | 26 | 5 | 20 |

Esto confirma el cumplimiento de la restricción:

```text
max_qubits = 28 ≤ 30
```

---

## 6. Formulación QUBO por bloques

El QUBO se construyó con codificación domain-wall. Para cada semana local se usan cuatro variables binarias:

```text
y_{t,0}, y_{t,1}, y_{t,2}, y_{t,3}
```

El nivel se recupera con la suma de bits y reparación domain-wall. El notebook generó los siguientes QUBO:

| block_id | weeks | qubits | linear_terms | quadratic_terms | target_sum_k_DP |
|---:|---:|---:|---:|---:|---:|
| 0 | 7 | 28 | 28 | 378 | -6 |
| 1 | 7 | 28 | 28 | 378 | -5 |
| 2 | 7 | 28 | 28 | 378 | 0 |
| 3 | 5 | 20 | 20 | 190 | 0 |

Energía QUBO de la política DP codificada:

| Bloque | Energía DP codificada |
|---:|---:|
| 0 | 0.0758668 |
| 1 | 0.0771668 |
| 2 | 0.0780793 |
| 3 | 0.0511204 |

Lectura técnica:

- El QUBO cabe en IBM por qubits, pero la cantidad de términos cuadráticos es alta.
- Cada término cuadrático se traduce en interacciones de fase controlada dentro del circuito QAOA.
- Esto produce circuitos profundos después de transpilación.

---

## 7. Configuración QAOA usada

El notebook construyó circuitos QAOA con:

```text
QAOA_P = 1
gamma = 0.8
beta = 0.35
shots = 2048
```

Punto importante: el notebook **no ejecutó una optimización clásica iterativa de ángulos**. Usó parámetros fijos conservadores. Por tanto, el resultado IBM debe interpretarse como una primera prueba funcional de hardware, no como una búsqueda variacional plenamente optimizada.

Circuitos brutos antes de transpilación:

| Bloque | Qubits | Depth bruto | Gates brutas |
|---:|---:|---:|---|
| 0 | 28 | 60 | 378 `cp`, 28 `h`, 28 `p`, 28 `rx`, 28 `measure` |
| 1 | 28 | 60 | 378 `cp`, 28 `h`, 28 `p`, 28 `rx`, 28 `measure` |
| 2 | 28 | 60 | 378 `cp`, 28 `h`, 28 `p`, 28 `rx`, 28 `measure` |
| 3 | 20 | 44 | 190 `cp`, 20 `h`, 20 `p`, 20 `rx`, 20 `measure` |

---

## 8. Ejecución en IBM Quantum

Metadatos del job IBM:

| Campo | Valor |
|---|---:|
| Job ID | `d92j4vd958jc73brsfg0` |
| Backend | `ibm_fez` |
| Qubits backend | 156 |
| Programa | `sampler` |
| Estado | `Completed` |
| Fecha de creación UTC | 2026-07-01T15:30:05.467936Z |
| Runtime estimado | 7.151241653 s |
| Cost | 115 |
| Número de pubs/circuitos | 4 |
| Shots por circuito | 2048 |

El archivo `job-result.json` contiene:

```text
PrimitiveResult
└── pub_results
    ├── SamplerPubResult bloque 0
    ├── SamplerPubResult bloque 1
    ├── SamplerPubResult bloque 2
    └── SamplerPubResult bloque 3
```

Cada bloque contiene un `BitArray` con el registro clásico `c`.

Tamaño confirmado de resultados:

| Bloque | Bits medidos | Shots | Interpretación |
|---:|---:|---:|---|
| 0 | 28 | 2048 | 7 semanas × 4 qubits |
| 1 | 28 | 2048 | 7 semanas × 4 qubits |
| 2 | 28 | 2048 | 7 semanas × 4 qubits |
| 3 | 20 | 2048 | 5 semanas × 4 qubits |

---

## 9. Circuitos transpilados en IBM

Profundidad ISA después de transpilación con `optimization_level=3`:

| Bloque | Qubits | Depth ISA | Compuertas transpiladas |
|---:|---:|---:|---|
| 0 | 28 | 1731 | 4461 `sx`, 2133 `rz`, 2112 `cz`, 37 `x`, 28 `measure` |
| 1 | 28 | 2133 | 4464 `sx`, 2135 `cz`, 1996 `rz`, 31 `x`, 28 `measure` |
| 2 | 28 | 1793 | 4497 `sx`, 2142 `cz`, 2066 `rz`, 48 `x`, 28 `measure` |
| 3 | 20 | 1132 | 2062 `sx`, 1041 `rz`, 968 `cz`, 30 `x`, 20 `measure` |

### Diagnóstico

Aunque los circuitos caben por número de qubits, son muy profundos para hardware NISQ. La cantidad de compuertas `cz` es especialmente relevante porque las compuertas de dos qubits suelen tener más error que las de un qubit.

Esto explica por qué la distribución de bitstrings fue dispersa y por qué el bitstring top no necesariamente corresponde a una buena solución QUBO.

---

## 10. Resultado QPU por bloques

Bitstrings top reportados por IBM:

```text
Bloque 0 top bitstring: 0111010010101000001000000110
Bloque 1 top bitstring: 0110000000000000001010110011
Bloque 2 top bitstring: 0101010100010000010110000010
Bloque 3 top bitstring: 11111001011100000110
```

Luego el notebook reparó domain-wall y decodificó los niveles semanales:

| block_id | source | bitstring_top | bitstring_repaired | shots_top | k_block | sum_k_block | qubits |
|---:|---|---|---|---:|---|---:|---:|
| 0 | IBM_QPU | `0111010010101000001000000110` | `1110100011001000100000001100` | 1 | [1, -1, 0, -1, -1, -2, 0] | -4 | 28 |
| 1 | IBM_QPU | `0110000000000000001010110011` | `1100000000000000100011101100` | 1 | [0, -2, -2, -2, -1, 1, 0] | -6 | 28 |
| 2 | IBM_QPU | `0101010100010000010110000010` | `1100110010000000110010001000` | 1 | [0, 0, -1, -2, 0, -1, -1] | -5 | 28 |
| 3 | IBM_QPU | `11111001011100000110` | `11111100111000001100` | 2 | [2, 0, 1, -2, 0] | 1 | 20 |

### Lectura estadística

Los bitstrings más frecuentes tuvieron muy baja frecuencia:

```text
shots_top = 1, 1, 1, 2
```

Con 2048 shots por bloque, esto significa que el muestreo no concentró probabilidad en una solución dominante. Por tanto, **no conviene seleccionar solo el bitstring más frecuente**. Es necesario evaluar un conjunto `top-k` de candidatos por energía QUBO, factibilidad y SRS global.

---

## 11. Política reconstruida desde IBM/QAOA

La política completa reconstruida desde los cuatro bloques fue:

```text
k_QAOA_IBM = [1, -1, 0, -1, -1, -2, 0, 0, -2, -2, -2, -1, 1, 0, 0, 0, -1, -2, 0, -1, -1, 2, 0, 1, -2, 0]
```

Ajustes en volumen:

```text
u_QAOA_IBM(t) = k_QAOA_IBM(t) × 2,706,476.21 m³/semana
```

Evaluación final:

```text
SRS reconstruido = -0.33216410584163464
Factibilidad:
  release_nonnegative = True
  storage_bounds      = True
  balance             = False
  all                 = False
```

Interpretación:

- La solución QPU no genera liberaciones negativas.
- La solución QPU mantiene límites de almacenamiento.
- La solución QPU falla el balance acumulado.
- La solución QPU es más irregular que DP, lo que eleva `Csmooth`.

---

## 12. Comparación final de modelos

| modelo | SRS | Delta_vs_regla | Delta_vs_DP | Ccrit | Cdev | Csmooth | sum_u_m3 | min_storage_pct_cap | min_release_m3 | factible | runtime_s | max_qubits |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|
| Histórico | -0.311534 | -0.021029 | -0.018374 | 5.685972e18 | 0 | 0 | 0 | 9.653041 | 859350 | True | 0.000000 | 0 |
| Regla umbral | -0.290505 | 0.000000 | 0.002655 | 4.845873e18 | 1.904504e14 | 0 | -7.036838e7 | 11.298953 | -1.847126e6 | False | 0.000000 | 0 |
| DP exacta | -0.293160 | -0.002655 | 0.000000 | 5.148445e18 | 8.057515e13 | 1.465003e13 | -2.977124e7 | 10.558293 | 584091 | True | 0.049052 | 0 |
| QUBO/QAOA IBM por bloques | -0.332164 | -0.041659 | -0.039004 | 5.144660e18 | 2.783505e14 | 4.028757e14 | -3.789067e7 | 10.805179 | 584091 | False | NaN | 28 |

### Lectura comparativa

1. **Regla umbral**: tiene el SRS más alto, pero es no factible por liberación negativa. No debe considerarse solución final.
2. **DP exacta**: mejor solución factible. Debe tratarse como benchmark principal.
3. **Histórico**: factible, pero con menor desempeño que DP.
4. **QUBO/QAOA IBM**: ejecutado correctamente, pero peor que DP y no factible por balance.

Deterioro IBM respecto a DP:

```text
Delta_vs_DP = -0.039004
penalización absoluta IBM ≈ 13.30% mayor que DP
```

Componentes que explican el deterioro:

```text
Cdev_IBM / Cdev_DP       ≈ 3.45
Csmooth_IBM / Csmooth_DP ≈ 27.50
```

Aunque `Ccrit` IBM es ligeramente menor que DP, la política IBM es mucho más costosa en desviación y suavidad. Por eso el SRS final empeora.

---

## 13. Vista semanal: DP vs IBM/QAOA

Primeras diez semanas mostradas en el notebook:

| week | date | Robs_m3 | Sobs_m3 | k_DP | u_DP_m3 | k_QAOA_IBM | u_QAOA_IBM_m3 |
|---:|---|---:|---:|---:|---:|---:|---:|
| 0 | 2025-06-29 | 859350 | 3.756112e8 | 0 | 0 | 1 | 2.706476e6 |
| 1 | 2025-07-06 | 3.290567e6 | 3.813530e8 | -1 | -2.706476e6 | -1 | -2.706476e6 |
| 2 | 2025-07-13 | 8.086461e6 | 3.810606e8 | -1 | -2.706476e6 | 0 | 0 |
| 3 | 2025-07-20 | 7.562678e6 | 4.018978e8 | -1 | -2.706476e6 | -1 | -2.706476e6 |
| 4 | 2025-07-27 | 8.349581e6 | 4.023278e8 | -1 | -2.706476e6 | -1 | -2.706476e6 |
| 5 | 2025-08-03 | 1.511289e7 | 3.953748e8 | -1 | -2.706476e6 | -2 | -5.412952e6 |
| 6 | 2025-08-10 | 1.726360e7 | 3.810244e8 | -1 | -2.706476e6 | 0 | 0 |
| 7 | 2025-08-17 | 1.429395e7 | 3.692016e8 | -1 | -2.706476e6 | 0 | 0 |
| 8 | 2025-08-24 | 1.360895e7 | 3.621061e8 | -1 | -2.706476e6 | -2 | -5.412952e6 |
| 9 | 2025-08-31 | 1.017766e7 | 3.542194e8 | -1 | -2.706476e6 | -2 | -5.412952e6 |

Lectura:

- DP aplica una política conservadora y relativamente suave.
- IBM/QAOA alterna entre aumentos, reducciones moderadas, reducciones fuertes y ceros.
- Esa irregularidad explica el alto `Csmooth` de IBM/QAOA.

---

## 14. Gráficas integradas

### 14.1 Almacenamiento observado vs DP vs IBM/QAOA

![Almacenamiento observado vs DP vs IBM](falcon_ibm_markdown_assets/grafica_almacenamiento_observado_dp_ibm.png)

Lectura:

- DP e IBM/QAOA pueden elevar el almacenamiento respecto a la trayectoria observada.
- IBM/QAOA no necesariamente mejora el SRS porque el objetivo penaliza también desviación y suavidad.
- Visualmente IBM puede parecer competitivo en almacenamiento, pero globalmente falla balance y suavidad.

### 14.2 Liberación histórica vs DP vs IBM/QAOA

![Liberación histórica vs DP vs IBM](falcon_ibm_markdown_assets/grafica_liberacion_observada_dp_ibm.png)

Lectura:

- DP modifica la liberación de forma controlada.
- IBM/QAOA genera una política más oscilante.
- El incremento de oscilación penaliza `Csmooth`, deteriorando el SRS final.

---

## 15. Diagnóstico integral

### 15.1 Lo que sí funcionó

El pipeline completo funcionó:

```text
Datos IBWC
→ limpieza y agregación semanal
→ cálculo de niveles discretos
→ benchmark histórico
→ regla umbral
→ DP exacta
→ QUBO por bloques
→ codificación domain-wall
→ QAOA p=1
→ transpilación a ibm_fez
→ ejecución Sampler en IBM Quantum
→ reconstrucción de política semanal
→ evaluación SRS y restricciones
```

Esto es una contribución sólida porque demuestra que el problema puede llevarse a hardware cuántico real respetando el límite de 30 qubits.

### 15.2 Lo que no funcionó todavía

IBM/QAOA no produjo una política final competitiva:

- SRS peor que DP.
- No factible por balance acumulado.
- Bitstrings top con frecuencia casi nula.
- Circuitos muy profundos.
- Ángulos QAOA no optimizados.

### 15.3 Causa probable 1: profundidad de circuito

Los circuitos tienen entre `1132` y `2133` de profundidad ISA. Esto es alto para hardware NISQ. Los bloques de 28 qubits requieren más de `2100` compuertas `cz`, lo que probablemente degrada la distribución de medición.

### 15.4 Causa probable 2: top-1 bitstring no es suficiente

En IBM, el bitstring más frecuente apareció solo 1 o 2 veces. En ese contexto, seleccionar el top-1 equivale a elegir una muestra casi aleatoria entre muchas. La solución debe reconstruirse evaluando muchos candidatos.

### 15.5 Causa probable 3: balance global debilitado por bloques

La restricción de balance es global sobre las 26 semanas. Al dividir en bloques, cada QUBO resuelve un subproblema local. Si los bloques no están acoplados por una cuota global fuerte, la política concatenada puede fallar el balance.

### 15.6 Causa probable 4: QAOA no optimizado

El notebook usó `p=1`, `gamma=0.8`, `beta=0.35`. Esto valida ejecución, pero no explora el paisaje variacional. La comparación contra DP será más justa cuando se optimicen ángulos por bloque.

---

## 16. Cómo debe integrarse con el análisis previo del proyecto

El análisis previo ya contenía:

- Formulación clásica del challenge.
- Justificación de SRS.
- DP exacta como benchmark.
- QUBO/domain-wall.
- Estrategia de bloques para ≤30 qubits.
- Preparación para IBM Qiskit Runtime.

Este análisis añade la evidencia experimental:

1. El job realmente se envió a IBM.
2. El backend fue `ibm_fez`.
3. El job se completó.
4. Se ejecutaron cuatro circuitos, uno por bloque.
5. Cada circuito usó 2048 shots.
6. Los bloques respetaron 28, 28, 28 y 20 bits/qubits medidos.
7. La reconstrucción QPU generó una política semanal completa.
8. La política fue evaluada con las mismas métricas SRS y restricciones.
9. El resultado IBM/QAOA quedó peor que DP y no factible.
10. El diagnóstico técnico apunta a profundidad, ruido, top-1 y balance global.

Narrativa recomendada:

> El trabajo evolucionó de una solución clásica validada a una prueba cuántica/híbrida real. El resultado IBM no busca demostrar ventaja cuántica, sino demostrar codificación, escalamiento, ejecución y evaluación bajo un benchmark común. La DP exacta se conserva como solución de referencia; IBM/QAOA se reporta como primera iteración experimental en hardware real.

---

## 17. Recomendaciones técnicas para la siguiente iteración

### 17.1 Reemplazar selección top-1 por evaluación top-k

En lugar de:

```python
bs = max(counts, key=counts.get)
```

usar:

```python
top_candidates = sorted(counts.items(), key=lambda x: x[1], reverse=True)[:K]
```

Para cada candidato:

1. Reparar domain-wall.
2. Decodificar `k_block`.
3. Calcular energía QUBO.
4. Concatenar candidatos entre bloques con búsqueda beam search.
5. Calcular SRS global.
6. Filtrar factibilidad.
7. Seleccionar la mejor política factible.

### 17.2 Beam search entre bloques

La mejor combinación global no necesariamente usa el top-1 de cada bloque. Se recomienda:

```text
Para cada bloque, conservar top K candidatos.
Hacer beam search sobre combinaciones de bloques.
Evaluar factibilidad global y SRS final.
```

Esto ataca directamente el fallo de balance acumulado.

### 17.3 Penalización global de balance

Aumentar la penalización asociada a:

```text
(Σ_t u(t))^2
```

o asignar cuotas por bloque derivadas de DP:

```text
target_sum_k_block = [-6, -5, 0, 0]
```

Pero no como referencia débil; debe convertirse en restricción o penalización más fuerte.

### 17.4 Reducir tamaño de bloque

Probar:

| Semanas por bloque | Qubits | Ventaja |
|---:|---:|---|
| 7 | 28 | Mayor horizonte local, pero circuitos profundos |
| 6 | 24 | Menor profundidad |
| 5 | 20 | Más estable en QPU y simulable localmente |
| 4 | 16 | Más robusto, aunque menos expresivo |

Una prueba de 5 semanas puede ser más defendible experimentalmente porque el circuito transpilado debería reducir considerablemente profundidad y ruido.

### 17.5 Optimizar ángulos QAOA

Comparar:

- `p=1` con ángulos fijos.
- `p=1` con COBYLA/SPSA.
- `p=2` con pocos pasos de optimización.
- Warm start desde DP.
- Warm start desde política de regla factibilizada.

### 17.6 Simuladores como control

Antes de QPU real, correr:

1. Simulador ideal.
2. Simulador con ruido del backend.
3. QPU real.

Esto permite separar errores de formulación, errores de ansatz y errores físicos.

### 17.7 Mitigación de ruido en IBM Runtime

Solicitar o activar, según disponibilidad:

- `optimization_level=3` ya usado.
- `resilience_level` si está disponible.
- Dynamical decoupling.
- Twirling.
- Readout mitigation.
- Comparación con backend alternativo de menor error de compuertas.

---

## 18. Párrafo listo para reporte ejecutivo

> Se implementó una solución híbrida para el Falcon Challenge A, escalando el problema desde benchmarks clásicos hacia una formulación QUBO/QAOA ejecutada en IBM Quantum. El horizonte oficial de 26 semanas se particionó en bloques de 7, 7, 7 y 5 semanas mediante codificación domain-wall, lo que permitió respetar el límite de 30 qubits, usando un máximo de 28 qubits por circuito. El job se ejecutó correctamente en el backend `ibm_fez` con 2048 shots por bloque. La DP exacta obtuvo el mejor resultado factible con SRS de `-0.293160`, mientras que la política reconstruida desde IBM/QAOA obtuvo SRS de `-0.332164` y falló la restricción de balance acumulado. Por tanto, el resultado cuántico actual debe interpretarse como una validación funcional de escalamiento e integración con hardware real, no como una mejora sobre el benchmark clásico. La siguiente iteración debe enfocarse en selección top-k de bitstrings, penalización global de balance, reducción de profundidad de circuitos y optimización variacional de ángulos QAOA.

---

## 19. Párrafo listo para presentación oral

> En esta etapa, el valor principal del proyecto fue demostrar que el problema sí puede llevarse a IBM Quantum respetando el límite de 30 qubits. La formulación clásica se validó con DP exacta, luego se transformó a QUBO con codificación domain-wall y finalmente se ejecutó por bloques en `ibm_fez`. Sin embargo, los resultados muestran que el hardware cuántico todavía no supera al modelo clásico: la DP exacta es factible y obtiene mejor SRS, mientras que IBM/QAOA falla el balance acumulado. Esto es consistente con las limitaciones actuales de hardware NISQ: los circuitos transpilados fueron profundos y los bitstrings medidos estuvieron muy dispersos. La contribución no es una ventaja cuántica, sino una ruta reproducible de escalamiento clásico-cuántico y una base clara para iteraciones futuras.

---

## 20. Qué pedirle a Claude a partir de este markdown

Prompt sugerido:

```text
Con base en este análisis integral, actualiza el reporte global del Falcon Challenge. Integra la ejecución real en IBM Quantum como evidencia experimental, mantén DP exacta como benchmark factible principal y explica que QAOA/IBM fue una prueba de escalamiento, no una mejora final. Reestructura la narrativa para una presentación técnica de 9 minutos y una sección de resultados para reporte académico. No afirmes ventaja cuántica. Incluye recomendaciones para la siguiente iteración: top-k bitstrings, beam search, penalización global de balance, reducción de bloques y optimización de ángulos QAOA.
```

---

## 21. Conclusión final integrada

El proyecto ya cuenta con una ruta metodológica completa y verificable:

```text
Datos hidrológicos IBWC
→ SRS oficial
→ restricciones operativas
→ benchmark histórico
→ regla clásica
→ DP exacta
→ QUBO por bloques
→ domain-wall ≤30 qubits
→ QAOA p=1
→ IBM Quantum Runtime Sampler
→ reconstrucción de política
→ evaluación comparativa
```

La conclusión técnica es:

1. El escalamiento a IBM Quantum fue exitoso desde el punto de vista de ejecución.
2. El límite de 30 qubits se cumplió.
3. La DP exacta sigue siendo el mejor resultado factible.
4. IBM/QAOA produjo una política completa, pero no factible por balance.
5. El bajo desempeño QPU se explica por profundidad, ruido, selección top-1 y falta de optimización variacional.
6. La siguiente fase debe mejorar post-procesamiento y formulación antes de buscar mejores resultados en hardware real.

Esta es la narrativa más sólida para integrar el resultado IBM con todo el análisis previo del Falcon Challenge.
