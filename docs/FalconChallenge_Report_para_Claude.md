---
title: "Reporte de resultados en LaTeX: Guided Challenge A"
subtitle: "Programación resiliente de descargas para el International Falcon Reservoir"
project: "Inteligencia Artificial"
implementation: "Reproducible en Google Colab/Python"
date: "2026-06-30"
source_pdf: "FalconChallenge_Report_LaTeX(1).pdf"
conversion_note: "Conversión a Markdown optimizada para lectura en Claude. Las figuras del PDF se representan con captions y descripciones textuales."
---

# Reporte de resultados en LaTeX: Guided Challenge A

## Programación resiliente de descargas para el International Falcon Reservoir

**Proyecto:** Inteligencia Artificial  
**Implementación:** Reproducible en Google Colab/Python  
**Fecha:** 30 de junio de 2026

## Resumen

Este documento reporta, con formato de entrega técnica, la solución del **Guided Challenge A: Resilient Release Scheduling for the International Falcon Reservoir**. El reto pide diseñar una política de ajustes semanales de descarga $u(t)$ que mejore la resiliencia del almacenamiento del embalse sin violar restricciones operativas.

La métrica oficial es el **Storage Resilience Score (SRS)**, que penaliza almacenamiento crítico, desviaciones respecto a la operación histórica y cambios bruscos entre semanas. Para el benchmark oficial de $T = 26$ semanas y $L = 5$ niveles discretos, la mejor política factible encontrada por programación dinámica exacta obtuvo:

$$
SRS = -0.290423
$$

frente al histórico:

$$
SRS_{hist} = -0.311534
$$

con una mejora:

$$
\Delta SRS = +0.021111
$$

El resultado no debe interpretarse como política operativa real; es una solución de benchmark sobre el modelo simplificado del challenge.

## Índice

1. [Contexto del challenge y alineación con la documentación anexa](#1-contexto-del-challenge-y-alineación-con-la-documentación-anexa)
2. [Información requerida por el challenge](#2-información-requerida-por-el-challenge)
   - [2.1. Variables de decisión y dinámica del sistema](#21-variables-de-decisión-y-dinámica-del-sistema)
   - [2.2. Métrica oficial: Storage Resilience Score](#22-métrica-oficial-storage-resilience-score)
   - [2.3. Restricciones de factibilidad](#23-restricciones-de-factibilidad)
   - [2.4. Benchmark oficial y entregables esperados](#24-benchmark-oficial-y-entregables-esperados)
3. [Datos y preprocesamiento](#3-datos-y-preprocesamiento)
4. [Métodos implementados](#4-métodos-implementados)
5. [Resultados del benchmark oficial: T = 26, L = 5](#5-resultados-del-benchmark-oficial-t--26-l--5)
6. [Visualización de resultados](#6-visualización-de-resultados)
7. [Escalamiento computacional](#7-escalamiento-computacional)
8. [Resultados esperados para la entrega del challenge](#8-resultados-esperados-para-la-entrega-del-challenge)
9. [Conclusiones](#9-conclusiones)
10. [Apéndice A. Calendario semanal optimizado](#apéndice-a-calendario-semanal-optimizado)
11. [Referencias](#referencias)

---

# 1. Contexto del challenge y alineación con la documentación anexa

El challenge se ubica en la operación de sistemas hídricos transfronterizos del Río Grande/Río Bravo, específicamente en el **International Falcon Reservoir**. La documentación del reto lo alinea con **ODS 6.4**, **ODS 6.5** y **ODS 13.1**: uso eficiente y sostenible del agua, gestión integrada del recurso hídrico con cooperación transfronteriza, y resiliencia ante estrés climático.

La pregunta concreta es si una optimización basada en datos puede identificar calendarios de descarga que mejoren la resiliencia del almacenamiento durante periodos de bajo almacenamiento.

La documentación del curso de Inteligencia Artificial sirve como marco metodológico: se plantea una función objetivo, se evalúan restricciones y se comparan métodos de optimización. Aunque aquí no se entrena una red neuronal, la lógica es análoga al diseño de modelos visto en el proyecto: definir una función de costo, evaluar el desempeño, comparar alternativas y justificar la estrategia de optimización. En este caso, la función a minimizar es el costo negativo asociado al SRS.

# 2. Información requerida por el challenge

## 2.1. Variables de decisión y dinámica del sistema

La variable de decisión es el ajuste de descarga semanal $u(t)$. La descarga optimizada se calcula como:

$$
R(t) = R^{obs}(t) + u(t)
\tag{1}
$$

donde $R^{obs}(t)$ es la descarga histórica observada. La operación histórica se reproduce con $u(t) = 0$.

La dinámica de almacenamiento simplificada del challenge es:

$$
S^{opt}(t + 1) = S^{opt}(t) + \Delta S^{obs}(t) - u(t)
\tag{2}
$$

donde $\Delta S^{obs}(t)$ es el cambio observado de almacenamiento. En la implementación usada en Colab, esto se calcula de forma equivalente como:

$$
S^{opt}(t) = S^{obs}(t) - \sum_{i=0}^{t-1} u(i)
\tag{3}
$$

## 2.2. Métrica oficial: Storage Resilience Score

El score oficial es:

$$
SRS = -\left(w_1 C_{crit} + w_2 C_{dev} + w_3 C_{smooth}\right)
\tag{4}
$$

con:

$$
C_{crit} = \sum_{t=0}^{T}\left[\max\left(0, S_{min} - S^{opt}(t)\right)\right]^2
\tag{5}
$$

$$
C_{dev} = \sum_{t=0}^{T-1} u(t)^2
\tag{6}
$$

$$
C_{smooth} = \sum_{t=1}^{T-1}\left[u(t) - u(t-1)\right]^2
\tag{7}
$$

La interpretación es directa:

- $w_1$ prioriza evitar almacenamiento críticamente bajo.
- $w_2$ evita soluciones que se alejen demasiado de la operación histórica.
- $w_3$ penaliza cambios abruptos de una semana a otra.

## 2.3. Restricciones de factibilidad

La secuencia $u(t)$ debe respetar:

$$
R(t) \ge 0
\tag{8}
$$

$$
|u(t)| \le u_{max}
\tag{9}
$$

$$
0 \le S^{opt}(t) \le S_{max}
\tag{10}
$$

$$
\left|\sum_{t=0}^{T-1} u(t)\right| \le \eta \sum_{t=0}^{T-1} R^{obs}(t)
\tag{11}
$$

La restricción de balance acumulado es clave: evita mejorar artificialmente el almacenamiento mediante una reducción sistemática de todas las descargas. Por diseño, fuerza al optimizador a redistribuir descargas en el tiempo.

## 2.4. Benchmark oficial y entregables esperados

El benchmark oficial usa $L = 5$ niveles de ajuste:

$$
u(t) \in \{-2\Delta u, -\Delta u, 0, \Delta u, 2\Delta u\}, \qquad
\Delta u = 0.25\,\tilde{R}^{obs}_{week}, \qquad
u_{max} = 2\Delta u
\tag{12}
$$

Los pesos oficiales son:

$$
w_1 = \frac{1}{(T+1)S_{scale}^2}, \qquad
w_2 = \frac{0.1}{T u_{max}^2}, \qquad
w_3 = \frac{0.1}{(T-1)(2u_{max})^2}, \qquad
S_{scale} = S_{min}
\tag{13}
$$

Los entregables esperados son:

1. La política $u_{opt}(t)$.
2. El cálculo de $R^{opt}(t)$ y $S^{opt}(t)$.
3. El SRS comparado contra baselines.
4. $\Delta SRS$.
5. Runtime.
6. Análisis de escalamiento para instancias pequeñas, medianas y grandes.

# 3. Datos y preprocesamiento

## 3.1. Fuentes utilizadas

Se usaron los archivos anexos del challenge:

- `DataSetExport-Discharge.Best Available@08461300-Instantaneous-m^3 s-...csv`: descarga observada bajo Falcon Dam, usada como $R^{obs}$.
- `DataSetExport-Total Storage.Web-Daily-tcm@08461200-Instantaneous-m^3-...csv`: almacenamiento observado del International Falcon Reservoir.
- `falcon_reservoir_constants.json`: constantes del embalse, incluyendo $S_{max}$.

La descarga reportada en $m^3/s$ se integró a volumen semanal en $m^3$. Las series se agregaron en semanas de domingo a domingo para obtener $T$ decisiones semanales y $T+1$ fronteras de almacenamiento.

## 3.2. Parámetros numéricos de la corrida oficial

**Tabla 1. Parámetros usados en el benchmark oficial reproducido.**

| Concepto | Valor |
|---|---:|
| Ventana oficial usada | 2025-06-29 a 2025-12-28 |
| Horizonte | $T = 26$ semanas |
| Niveles de ajuste | $L = 5$, $u(t) \in \{-2\Delta u, -\Delta u, 0, \Delta u, 2\Delta u\}$ |
| Mediana semanal de $R^{obs}$ | 11.217 Mm³ |
| $\Delta u$ | 2.804 Mm³ por semana |
| $u_{max}$ | 5.609 Mm³ por semana |
| $S_{max}$ | 3288.726 Mm³ |
| $S_{min} = 0.25S_{max}$ | 822.182 Mm³ |
| $\eta$ balance acumulado | 0.10 |
| Presupuesto de balance $\eta \sum R^{obs}$ | 31.151 Mm³ |
| Volumen observado total $\sum R^{obs}$ | 311.506 Mm³ |
| Volumen optimizado total $\sum R^{opt}$ | 280.658 Mm³ |

# 4. Métodos implementados

## 4.1. Histórico

La referencia histórica usa $u(t) = 0$ para toda semana. Por tanto:

$$
R^{hist}(t) = R^{obs}(t)
$$

$$
S^{opt}(t) = S^{obs}(t)
$$

## 4.2. Regla clásica de umbral

La regla de conservación indicada por el challenge reduce una unidad de ajuste cuando el almacenamiento está por debajo del umbral crítico:

$$
u_{rule}(t) =
\begin{cases}
-\Delta u, & S_{rule}(t) < S_{min} \\
0, & S_{rule}(t) \ge S_{min}
\end{cases}
\tag{14}
$$

Además de la regla literal, se reporta una versión balanceada que detiene reducciones cuando la suma acumulada de $u(t)$ violaría la restricción oficial de balance.

## 4.3. Programación dinámica exacta

La programación dinámica usa el hecho de que, una vez discretizados los niveles, el almacenamiento depende de la suma acumulada de niveles enteros. Para:

$$
u(t) = k_t \Delta u, \qquad k_t \in \{-2, -1, 0, 1, 2\}
$$

el estado se define como:

$$
(t, K_t, k_{t-1}), \qquad K_t = \sum_{i=0}^{t-1} k_i
\tag{15}
$$

Este estado conserva la información necesaria para calcular almacenamiento, balance acumulado y suavidad. En lugar de enumerar $L^T$ calendarios completos, se propagan sólo estados factibles y se poda cualquier trayectoria que ya no pueda cumplir las restricciones restantes.

## 4.4. Recocido simulado quantum-inspired

También se implementó una heurística de recocido simulado como aproximación **quantum-inspired**. Esta heurística representa cada semana por un nivel discreto de descarga, evalúa la misma función de costo y aplica penalizaciones grandes cuando viola restricciones. Sirve como comparación híbrida/heurística, pero no garantiza óptimo global. El benchmark factible se reporta con la DP exacta.

## 4.5. Estructura QUBO/Ising sugerida

Para una formulación cuántica o híbrida, se puede definir una variable binaria $x_{t,\ell} \in \{0, 1\}$ que indica si en la semana $t$ se selecciona el nivel $q_\ell \in \{-2, -1, 0, 1, 2\}$. Entonces:

$$
u(t) = \Delta u \sum_{\ell=1}^{L} q_\ell x_{t,\ell}, \qquad
\sum_{\ell=1}^{L} x_{t,\ell} = 1
\tag{16}
$$

La desviación $C_{dev}$ y la suavidad $C_{smooth}$ son cuadráticas en variables binarias. Las restricciones de selección única, almacenamiento, descarga no negativa y balance acumulado se pueden incorporar como penalizaciones. El término de almacenamiento crítico con $\max(0, \cdot)$ puede manejarse mediante variables auxiliares/slacks o aproximarse con penalización cuadrática por déficit.

Esto justifica que el problema sea compatible con QUBO/Ising, QAOA, quantum annealing o métodos híbridos.

# 5. Resultados del benchmark oficial: T = 26, L = 5

## 5.1. Tabla comparativa de políticas

**Tabla 2. Resumen de resultados oficiales para la ventana 2025-06-29 a 2025-12-28. Volúmenes en Mm³.**

| Política | SRS | Factible | R ≥ 0 | 0 ≤ S ≤ Smax | Balance | Σu | mín S %cap | mín R | Runtime (s) | ΔSRS vs hist. |
|---|---:|:---:|:---:|:---:|:---:|---:|---:|---:|---:|---:|
| Regla umbral | -0.288934 | No | Sí | Sí | No | -72.913 | 11.358 | 0.501 | 0.000 | 0.022600 |
| Óptimo DP exacto | -0.290423 | Sí | Sí | Sí | Sí | -30.848 | 10.591 | 0.501 | 0.117 | 0.021111 |
| Regla umbral balanceada | -0.290423 | Sí | Sí | Sí | Sí | -30.848 | 10.591 | 0.501 | 0.000 | 0.021111 |
| Annealing quantum-inspired | -0.291188 | Sí | Sí | Sí | Sí | -30.848 | 10.591 | 0.501 | 6.594 | 0.020346 |
| Histórico $u=0$ | -0.311534 | Sí | Sí | Sí | Sí | 0.000 | 9.653 | 3.305 | 0.000 | 0.000000 |

## 5.2. Resultado principal

El mejor resultado factible es la política de programación dinámica exacta:

$$
SRS_{DP} = -0.290423, \qquad SRS_{hist} = -0.311534, \qquad \Delta SRS = +0.021111
\tag{17}
$$

La mejora relativa frente al histórico ocurre porque la política reduce descargas al inicio de la ventana, elevando el almacenamiento optimizado respecto al observado, pero manteniendo descarga no negativa y respetando el balance acumulado.

La regla de umbral literal alcanza un SRS levemente mayor ($-0.288934$), pero no es admisible porque viola el balance acumulado: reduce 72.913 Mm³, por encima del presupuesto permitido de 31.151 Mm³. Por esa razón no debe reportarse como solución oficial factible.

## 5.3. Descomposición del costo

**Tabla 3. Contribuciones normalizadas al costo del SRS para histórico y DP exacta.**

| Política | $w_1C_{crit}$ | $w_2C_{dev}$ | $w_3C_{smooth}$ | Costo total | SRS |
|---|---:|---:|---:|---:|---:|
| Histórico $u = 0$ | 0.311534 | 0.000000 | 0.000000 | 0.311534 | -0.311534 |
| Óptimo DP exacto | 0.279596 | 0.010577 | 0.000250 | 0.290423 | -0.290423 |

La DP reduce el componente crítico de 0.311534 a 0.279596. Paga una penalización de desviación 0.010577 y una penalización de suavidad 0.000250, pero el balance neto sigue siendo favorable: el costo total baja de 0.311534 a 0.290423.

## 5.4. Interpretación operacional de la política óptima

La política óptima aplica:

$$
u(t) = -\Delta u = -2.804 \text{ Mm}^3
$$

durante 11 semanas consecutivas y luego regresa a $u(t) = 0$ durante las 15 semanas restantes. La secuencia de niveles es:

$$
(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
\tag{18}
$$

El ajuste acumulado es:

$$
\sum u = -30.848 \text{ Mm}^3
$$

valor que queda dentro del límite permitido $\pm 31.151$ Mm³. El volumen total optimizado liberado es 280.658 Mm³, frente a 311.506 Mm³ en el histórico.

El almacenamiento mínimo optimizado sube de 317.462 Mm³ observado a 348.310 Mm³, equivalente a pasar de 9.653% a 10.591% de la capacidad de conservación.

# 6. Visualización de resultados

## Figura 1. Trayectoria semanal de almacenamiento observado y optimizado

**Caption original:** Trayectoria semanal de almacenamiento observado y optimizado. La línea punteada marca $S_{min} = 25\%S_{max}$.

**Descripción para Markdown/Claude:** El gráfico compara dos series temporales entre julio y diciembre de 2025: almacenamiento observado y almacenamiento optimizado por DP exacta. La trayectoria optimizada se mantiene por encima de la observada durante la mayor parte de la ventana, especialmente en el tramo crítico de bajo almacenamiento. La línea horizontal punteada de $S_{min}$ está alrededor de 822 Mm³, por encima de ambas trayectorias, mostrando que toda la ventana opera bajo el umbral crítico de 25% de capacidad.

## Figura 2. Descargas semanales observadas, optimizadas y ajuste $u(t)$

**Caption original:** Descargas semanales observadas, descargas optimizadas y ajuste $u(t)$.

**Descripción para Markdown/Claude:** El gráfico muestra las descargas observadas históricas, las descargas optimizadas y barras del ajuste $u(t)$. En las primeras 11 semanas, $u(t)$ es negativo y reduce la descarga optimizada respecto a la histórica; después, el ajuste vuelve a cero y la descarga optimizada coincide con la histórica.

# 7. Escalamiento computacional

## 7.1. Resultados de escalamiento

**Tabla 4. Instancias de escalamiento requeridas por el challenge.**

| T | L | $L^T$ | $\log_{10}(L^T)$ | SRS hist. | SRS DP | ΔSRS | Factible | Runtime (s) | Estados finales |
|---:|---:|---:|---:|---:|---:|---:|:---:|---:|---:|
| 12 | 3 | $3^{12}$ | 5.725 | -0.296264 | -0.296264 | 0.000000 | Sí | 0.004 | 27 |
| 26 | 5 | $5^{26}$ | 18.173 | -0.311534 | -0.290423 | 0.021111 | Sí | 0.093 | 115 |
| 52 | 5 | $5^{52}$ | 36.346 | -0.252702 | -0.201466 | 0.051236 | Sí | 0.514 | 135 |
| 52 | 7 | $7^{52}$ | 43.945 | -0.252702 | -0.177643 | 0.075059 | Sí | 1.056 | 180 |

## Figura 3. Tiempo observado de DP exacta frente al tamaño combinatorio teórico

**Caption original:** Tiempo observado de DP exacta frente al tamaño combinatorio teórico del espacio de búsqueda.

**Descripción para Markdown/Claude:** El gráfico presenta el runtime de DP exacta contra $\log_{10}(L^T)$ para cuatro instancias: $T=12,L=3$; $T=26,L=5$; $T=52,L=5$; y $T=52,L=7$. El tiempo observado aumenta con el tamaño combinatorio teórico, pero permanece bajo para estas instancias gracias a la representación de estados acumulados y la poda de trayectorias no factibles.

El espacio bruto de calendarios crece como $L^T$. Por ejemplo, el benchmark oficial $5^{26}$ equivale a:

$$
\log_{10}(5^{26}) = 18.173
$$

mientras que la instancia $7^{52}$ llega a:

$$
\log_{10}(7^{52}) = 43.945
$$

La DP exacta no enumera todas esas combinaciones; usa la estructura acumulativa del almacenamiento y poda estados no factibles. Por eso el runtime observado se mantiene bajo en estas instancias. En problemas más realistas con más variables de estado, restricciones legales o objetivos multiusuario, el crecimiento de estados puede ser mucho mayor.

# 8. Resultados esperados para la entrega del challenge

Para esta corrida, una entrega técnicamente completa debe reportar:

1. **Definición de datos y ventana:** estación 08461300 para descarga histórica, estación 08461200 para almacenamiento, conversión de $m^3/s$ a volumen semanal, ventana 2025-06-29 a 2025-12-28.
2. **Benchmark oficial:** $T = 26$, $L = 5$, $\Delta u = 2.804$ Mm³, $u_{max} = 5.609$ Mm³, $S_{max} = 3288.726$ Mm³, $S_{min} = 822.182$ Mm³, $\eta = 0.10$.
3. **Política optimizada:** secuencia $u_{opt}(t)$, descarga $R^{opt}(t)$ y almacenamiento $S^{opt}(t)$ para cada semana.
4. **Comparación de SRS:** histórico, regla de umbral, regla balanceada, DP exacta y annealing quantum-inspired.
5. **Factibilidad:** descarga no negativa, almacenamiento dentro de límites, y restricción de balance acumulado.
6. **Métrica principal:** $\Delta SRS$ frente a la referencia seleccionada. Frente al histórico, el resultado esperado es $+0.021111$ para la DP exacta.
7. **Runtime:** DP exacta alrededor de 0.117 s en la corrida oficial guardada; annealing quantum-inspired alrededor de 6.594 s.
8. **Escalamiento:** resultados para $T = 12, L = 3$; $T = 26, L = 5$; $T = 52, L = 5$; y $T = 52, L = 7$.
9. **Limitaciones:** el modelo es simplificado y no representa una política oficial de operación del embalse.

# 9. Conclusiones

El benchmark muestra que una redistribución temporal sencilla de descargas puede mejorar el SRS sin violar restricciones. La política factible óptima reduce descargas en las primeras 11 semanas, cuando el almacenamiento está críticamente bajo, y después vuelve a la operación histórica. Esta estrategia aumenta el almacenamiento mínimo optimizado y reduce el costo crítico, con una penalización pequeña por desviación y suavidad.

La principal conclusión computacional es que el problema tiene una estructura combinatoria adecuada para métodos exactos, heurísticos y quantum-inspired. La DP exacta es el mejor método para validar el benchmark oficial en esta instancia. El annealing quantum-inspired ofrece una implementación comparable en estructura a una formulación QUBO/Ising, pero para esta ventana no supera a la DP exacta.

La comparación de escalamiento confirma que el espacio bruto $L^T$ crece rápidamente, lo que justifica explorar técnicas híbridas para extensiones con más niveles, horizontes más largos o restricciones operativas adicionales.

# Apéndice A. Calendario semanal optimizado

**Tabla 5. Calendario detallado de la política DP exacta. Volúmenes en Mm³.**

| Semana | Inicio | Fin | $R^{obs}$ | $u^*$ | $R^{opt}$ | $S^{opt}_{inicio}$ | $S^{opt}_{fin}$ | Nivel |
|---:|---|---|---:|---:|---:|---:|---:|---:|
| 1 | 2025-06-29 | 2025-07-06 | 3.305 | -2.804 | 0.501 | 375.611 | 384.157 | -1 |
| 2 | 2025-07-06 | 2025-07-13 | 8.212 | -2.804 | 5.407 | 384.157 | 386.669 | -1 |
| 3 | 2025-07-13 | 2025-07-20 | 6.268 | -2.804 | 3.463 | 386.669 | 410.311 | -1 |
| 4 | 2025-07-20 | 2025-07-27 | 9.475 | -2.804 | 6.670 | 410.311 | 413.545 | -1 |
| 5 | 2025-07-27 | 2025-08-03 | 12.063 | -2.804 | 9.259 | 413.545 | 409.397 | -1 |
| 6 | 2025-08-03 | 2025-08-10 | 18.426 | -2.804 | 15.622 | 409.397 | 397.851 | -1 |
| 7 | 2025-08-10 | 2025-08-17 | 15.911 | -2.804 | 13.106 | 397.851 | 388.832 | -1 |
| 8 | 2025-08-17 | 2025-08-24 | 11.267 | -2.804 | 8.463 | 388.832 | 384.541 | -1 |
| 9 | 2025-08-24 | 2025-08-31 | 12.183 | -2.804 | 9.379 | 384.541 | 379.459 | -1 |
| 10 | 2025-08-31 | 2025-09-07 | 11.759 | -2.804 | 8.954 | 379.459 | 377.013 | -1 |
| 11 | 2025-09-07 | 2025-09-14 | 10.699 | -2.804 | 7.894 | 377.013 | 393.341 | -1 |
| 12 | 2025-09-14 | 2025-09-21 | 10.382 | 0.000 | 10.382 | 393.341 | 390.029 | 0 |
| 13 | 2025-09-21 | 2025-09-28 | 6.379 | 0.000 | 6.379 | 390.029 | 409.198 | 0 |
| 14 | 2025-09-28 | 2025-10-05 | 12.323 | 0.000 | 12.323 | 409.198 | 405.471 | 0 |
| 15 | 2025-10-05 | 2025-10-12 | 12.445 | 0.000 | 12.445 | 405.471 | 399.626 | 0 |
| 16 | 2025-10-12 | 2025-10-19 | 11.168 | 0.000 | 11.168 | 399.626 | 393.697 | 0 |
| 17 | 2025-10-19 | 2025-10-26 | 18.948 | 0.000 | 18.948 | 393.697 | 382.089 | 0 |
| 18 | 2025-10-26 | 2025-11-02 | 20.578 | 0.000 | 20.578 | 382.089 | 368.846 | 0 |
| 19 | 2025-11-02 | 2025-11-09 | 22.681 | 0.000 | 22.681 | 368.846 | 356.477 | 0 |
| 20 | 2025-11-09 | 2025-11-16 | 15.771 | 0.000 | 15.771 | 356.477 | 348.310 | 0 |
| 21 | 2025-11-16 | 2025-11-23 | 10.733 | 0.000 | 10.733 | 348.310 | 357.964 | 0 |
| 22 | 2025-11-23 | 2025-11-30 | 8.417 | 0.000 | 8.417 | 357.964 | 382.421 | 0 |
| 23 | 2025-11-30 | 2025-12-07 | 8.793 | 0.000 | 8.793 | 382.421 | 384.779 | 0 |
| 24 | 2025-12-07 | 2025-12-14 | 10.689 | 0.000 | 10.689 | 384.779 | 383.409 | 0 |
| 25 | 2025-12-14 | 2025-12-21 | 10.392 | 0.000 | 10.392 | 383.409 | 388.355 | 0 |
| 26 | 2025-12-21 | 2025-12-28 | 12.242 | 0.000 | 12.242 | 388.355 | 418.515 | 0 |

# Referencias

1. Guided Quantum Computing Challenges for Transboundary Water Systems. **Guided Challenge A: Resilient Release Scheduling for the International Falcon Reservoir**, 24 de junio de 2026.
2. Manifest del dataset Falcon: estaciones IBWC 08461200 y 08461300, archivos descargados y notas de uso.
3. Constantes del International Falcon Reservoir: capacidad total de conservación $S_{max} = 3288.726$ Mm³ y umbral $S_{min} = 822.1815$ Mm³.
4. Documentación anexa del curso Inteligencia Artificial: optimización de funciones de costo, comparación de métricas, uso de Google Colab y diseño reproducible de modelos.
