# CLAUDE.md — Falcon Challenge: EDA, Benchmark Clásico y Formulación QUBO/QAOA

Este archivo es el **input de trabajo para Claude** dentro del proyecto `FalconChallenge`. Resume lo más importante del PDF del reto, de los CSV anexos y de los notebooks de referencia sobre QUBO:

- `QUBO_Mathematical_Definition.ipynb`
- `QUBO_PenaltyMethod.ipynb`
- `QUBO_Examples_MaximumCut.ipynb`

El objetivo es que Claude implemente un pipeline reproducible para:

1. analizar los datos del embalse Falcon;
2. construir el dataset semanal requerido por el benchmark;
3. implementar baselines clásicos;
4. formular el problema como QUBO;
5. resolverlo primero con métodos clásicos y después con QAOA / optimización híbrida.

---

## 1. Reglas de seguridad sobre datos

Claude debe seguir estas reglas estrictamente.

1. **No modificar archivos originales dentro de `data/` o `data/source_csv/`.**
2. No borrar, renombrar, mover ni sobrescribir archivos originales.
3. No descargar datos de internet sin autorización explícita.
4. No usar archivos `.htm`, `.gif`, `.json`, `.txt`, `.pdf` como sustitutos de las series numéricas principales.
5. Los archivos procesados deben guardarse únicamente en:
   - `data/processed/`
   - `results/`
6. El código fuente nuevo debe guardarse en:
   - `src/`
7. Si falta un dataset oficial, debe reportarse claramente y usar una aproximación sólo como **preliminar / no oficial**.
8. No ejecutar `git add`, `git commit`, `git push` ni comandos destructivos.
9. No instalar paquetes nuevos sin permiso.
10. No presentar resultados sintéticos o derivados como benchmark oficial.

---

## 2. Objetivo del reto

El reto busca optimizar la política de liberaciones del **International Falcon Reservoir**, parte del sistema binacional México–Estados Unidos del Río Grande / Río Bravo.

La variable de decisión es el ajuste semanal de liberación:

\[
u(t)
\]

La liberación optimizada se define como:

\[
R(t) = R_{obs}(t) + u(t)
\]

donde:

- \(R_{obs}(t)\): descarga histórica observada.
- \(u(t)\): ajuste decidido por el optimizador.
- Si \(u(t)=0\), se reproduce la operación histórica.

La dinámica simplificada del almacenamiento es:

\[
S_{opt}(t+1)=S_{opt}(t)+\Delta S_{obs}(t)-u(t)
\]

donde:

- \(S_{opt}(t)\): almacenamiento optimizado.
- \(\Delta S_{obs}(t)\): cambio observado de almacenamiento.

El objetivo oficial es maximizar el **Storage Resilience Score**:

\[
SRS = -\left(w_1 C_{crit} + w_2 C_{dev} + w_3 C_{smooth}\right)
\]

Equivalentemente, para QUBO se minimiza:

\[
J = w_1 C_{crit} + w_2 C_{dev} + w_3 C_{smooth} + \text{penalizaciones}
\]

---

## 3. Términos de costo del benchmark

### 3.1 Penalización por almacenamiento crítico

\[
C_{crit}=\sum_{t=0}^{T}\left[\max\left(0,S_{min}-S_{opt}(t)\right)\right]^2
\]

Este término penaliza únicamente las semanas donde el almacenamiento cae por debajo del umbral crítico \(S_{min}\).

### 3.2 Penalización por desviación respecto a la operación histórica

\[
C_{dev}=\sum_{t=0}^{T-1}u(t)^2
\]

Este término evita cambios excesivos respecto a la descarga histórica.

### 3.3 Penalización por cambios bruscos entre semanas

\[
C_{smooth}=\sum_{t=1}^{T-1}\left[u(t)-u(t-1)\right]^2
\]

Este término favorece ajustes operativamente suaves.

---

## 4. Restricciones oficiales

La secuencia \(u(t)\) debe respetar:

\[
R(t) \ge 0
\]

\[
|u(t)| \le u_{max}
\]

\[
0 \le S_{opt}(t) \le S_{max}
\]

\[
\left|\sum_{t=0}^{T-1}u(t)\right| \le \eta \sum_{t=0}^{T-1}R_{obs}(t)
\]

La última restricción evita soluciones que mejoran el almacenamiento simplemente reduciendo sistemáticamente las liberaciones. La intención es redistribuir liberaciones en el tiempo, no retener agua de forma artificial.

---

## 5. Parámetros oficiales del benchmark

Para el benchmark oficial:

\[
L=5
\]

\[
u(t) \in \{-2\Delta u,-\Delta u,0,\Delta u,2\Delta u\}
\]

\[
\Delta u = 0.25 R^{obs}_{week}
\]

\[
u_{max}=2\Delta u
\]

\[
\eta=0.10
\]

\[
S_{min}=0.25S_{max}
\]

\[
S_{scale}=S_{min}
\]

Pesos oficiales:

\[
w_1=\frac{1}{(T+1)S_{scale}^2}
\]

\[
w_2=\frac{0.1}{T u_{max}^2}
\]

\[
w_3=\frac{0.1}{(T-1)(2u_{max})^2}
\]

---

## 6. Instancias de escalamiento

Implementar y reportar al menos:

| Instancia | T | L | Uso |
|---|---:|---:|---|
| Small | 12 semanas | 3 niveles | Validación / debugging |
| Medium | 26 semanas | 5 niveles | Benchmark oficial |
| Large | 52 semanas | 5 o 7 niveles | Escalamiento |

El número de agendas candidatas crece como:

\[
N_{schedules}=L^T
\]

---

## 7. Datasets requeridos

### 7.1 Datasets oficiales mínimos

Para construir el benchmark oficial se necesitan estas series:

| Variable | Dataset oficial | Estación | Uso |
|---|---|---|---|
| \(S_{obs}(t)\) | `Total Storage.Web-Daily-tcm@08461200` | 08461200 | almacenamiento observado |
| \(\Delta S_{obs}(t)\) | `Discharge.Total.Change-in-Storage@08461200` | 08461200 | cambio observado de almacenamiento |
| \(R_{obs}(t)\) | `Discharge.Best Available@08461300` | 08461300 | descarga histórica / release |
| \(S_{max}\) | Falcon total conservation storage capacity | reservoir overview | capacidad máxima para definir \(S_{min}\) |

### 7.2 CSV anexos disponibles en este proyecto

Los archivos CSV anexos actualmente son:

1. `DataSetExport-Discharge.Best Available@08461300-Instantaneous-m^3 s-20260629185451.csv`
   - Contiene `Timestamp (UTC-06:00)` y `Value (m^3/s)`.
   - Sirve para \(R_{obs}(t)\).
   - Debe convertirse de caudal a volumen semanal:

   \[
   R_{week,m^3}=R_{m^3/s}\times604800
   \]

2. `DataSetExport-Total Storage.Web-Daily-tcm@08461200-Instantaneous-m^3-20260629185416.csv`
   - Contiene `Timestamp (UTC-06:00)` y `Value (m^3)`.
   - Sirve para \(S_{obs}(t)\).

3. `DataSetExport-Reservoir Elevation.Web-Daily-m@08461200-Instantaneous-m-20260629185508.csv`
   - Contiene elevación del embalse en metros.
   - Útil para EDA y validación, no es indispensable para el QUBO mínimo.

4. `DataSetExport-Lake Area.Best Available@08461200-Instantaneous-m^2-20260629185344.csv`
   - Contiene área del lago en \(m^2\).
   - Útil para EDA, evaporación o extensiones, no es indispensable para el QUBO mínimo.

### 7.3 Faltante importante

No se detectó el CSV oficial:

```text
Discharge.Total.Change-in-Storage@08461200
```

Por tanto, si este archivo no está disponible, Claude debe construir una columna preliminar:

\[
\Delta S_{obs}(t) \approx S_{obs}(t+1)-S_{obs}(t)
\]

Debe etiquetarse como:

```text
derived_from_total_storage_not_official
```

Esta aproximación permite validar el pipeline, pero **no es benchmark oficial**.

---

## 8. Dataset semanal esperado

Claude debe construir:

```text
data/processed/falcon_weekly_benchmark.csv
```

con columnas mínimas:

```text
week
week_start
week_end
S_obs_m3
DeltaS_obs_m3
R_obs_m3_week
DeltaS_source
```

Reglas de agregación sugeridas:

- `S_obs_m3`: último valor semanal o promedio semanal. Preferir último valor si se interpreta como almacenamiento al cierre de semana.
- `R_obs_m3_week`: promedio semanal del caudal \(m^3/s\) multiplicado por 604800 segundos.
- `DeltaS_obs_m3`: si existe dataset oficial, usarlo convertido a volumen semanal; si no, derivarlo como diferencia semanal de `S_obs_m3`.
- Alinear las series por semana y eliminar semanas con valores críticos faltantes.

---

## 9. Baselines clásicos que debe implementar Claude

### 9.1 Histórico

\[
u_{hist}(t)=0
\]

\[
R_{hist}(t)=R_{obs}(t)
\]

### 9.2 Regla de umbral

\[
u_{rule}(t)=
\begin{cases}
-\Delta u, & S_{rule}(t)<S_{min} \\
0, & S_{rule}(t)\ge S_{min}
\end{cases}
\]

\[
R_{rule}(t)=R_{obs}(t)+u_{rule}(t)
\]

### 9.3 Clásico fuerte opcional

Implementar al menos uno:

- búsqueda exhaustiva para `T=12, L=3`;
- simulated annealing;
- dynamic programming;
- MILP / MIQP;
- evolutionary algorithm.

Para `T=12, L=3`, búsqueda exhaustiva es factible:

\[
3^{12}=531441
\]

Para `T=26, L=5`, búsqueda exhaustiva no es práctica:

\[
5^{26}\approx1.49\times10^{18}
\]

---

# 10. Principios QUBO tomados de los notebooks de referencia

## 10.1 Definición matemática

Un problema QUBO minimiza:

\[
\min_{x\in\{0,1\}^n} x^TQx
\]

Equivalente a:

\[
f(x)=\sum_i Q_{ii}x_i+\sum_{i<j}Q_{ij}x_ix_j
\]

Notas para implementación:

- Los términos lineales van en la diagonal de \(Q\).
- Los términos cuadráticos van fuera de la diagonal.
- Puede usarse matriz triangular superior.
- Para variables binarias:

\[
x_i^2=x_i
\]

Por eso un término lineal puede colocarse en la diagonal.

Si el problema original maximiza una función, convertirlo a minimización:

\[
\max f(x)=\min -f(x)
\]

## 10.2 Método de penalización

Para convertir restricciones en QUBO, agregar penalizaciones:

\[
f(x)+\sum_iP_i g_i(x)
\]

Cada \(g_i(x)\) debe ser:

- cero si la restricción se cumple;
- positivo si la restricción se viola.

Para igualdad lineal:

\[
\sum_i a_iy_i=b
\]

usar:

\[
P\left(\sum_i a_iy_i-b\right)^2
\]

Para desigualdad:

\[
\sum_i a_iy_i\le b
\]

introducir slack no negativo \(\eta\):

\[
\sum_i a_iy_i+\eta=b
\]

penalizar:

\[
P\left(\sum_i a_iy_i+\eta-b\right)^2
\]

El slack debe convertirse a binario.

## 10.3 Transformación de variables enteras a binarias

Para una variable entera acotada:

\[
\underline y_i\le y_i\le\overline y_i
\]

usar expansión binaria:

\[
y_i=\underline y_i+\sum_{j=0}^{N-2}2^jx_j^i+
\left(\overline y_i-\underline y_i-\sum_{j=0}^{N-2}2^j\right)x_{N-1}^i
\]

con:

\[
N=\left\lceil\log_2(\overline y_i-\underline y_i+1)\right\rceil
\]

## 10.4 Analogía con Max-Cut

En el ejemplo de Max-Cut se define una variable binaria por decisión y se traduce la contribución de cada decisión al costo QUBO.

Para una arista \((i,j)\), el conteo de corte es:

\[
x_i+x_j-2x_ix_j
\]

Como Max-Cut maximiza, se minimiza el negativo:

\[
-x_i-x_j+2x_ix_j
\]

Este patrón sirve como analogía para Falcon:

1. definir variables binarias para cada decisión semanal;
2. escribir \(u(t)\), \(S(t)\) y restricciones como expresiones lineales en esas binarias;
3. elevar al cuadrado para obtener términos cuadráticos;
4. cargar coeficientes en la matriz \(Q\).

---

# 11. Formulación QUBO recomendada para Falcon

## 11.1 Variable de decisión discreta

Para cada semana \(t\) y nivel \(\ell\), definir una variable binaria one-hot:

\[
x_{t,\ell}\in\{0,1\}
\]

\[
x_{t,\ell}=1 \Longleftrightarrow u(t)=a_\ell
\]

con niveles:

### Instancia oficial L = 5

\[
a_\ell\in\{-2\Delta u,-\Delta u,0,\Delta u,2\Delta u\}
\]

### Instancia small L = 3

\[
a_\ell\in\{-\Delta u,0,\Delta u\}
\]

La condición one-hot es:

\[
\sum_{\ell=0}^{L-1}x_{t,\ell}=1
\]

Entonces:

\[
u_t=\sum_{\ell=0}^{L-1}a_\ell x_{t,\ell}
\]

Número de bits de decisión:

\[
N_{decision}=T\times L
\]

Ejemplos:

- `T=12, L=3` → 36 bits de decisión.
- `T=26, L=5` → 130 bits de decisión.

## 11.2 Penalización one-hot

Agregar:

\[
P_{onehot}\sum_{t=0}^{T-1}\left(\sum_{\ell=0}^{L-1}x_{t,\ell}-1\right)^2
\]

Expansión para cada semana:

\[
\left(\sum_\ell x_{t,\ell}-1\right)^2
=1-\sum_\ell x_{t,\ell}+2\sum_{\ell<m}x_{t,\ell}x_{t,m}
\]

El término constante se puede ignorar en \(Q\).

Agregar a \(Q\):

- diagonal: \(-P_{onehot}\) para cada \(x_{t,\ell}\);
- off-diagonal: \(2P_{onehot}\) para cada par \((x_{t,\ell},x_{t,m})\), \(\ell<m\).

## 11.3 Expresión lineal de almacenamiento

Definir almacenamiento histórico acumulado sin ajustes:

\[
H_t=S_0+\sum_{k=0}^{t-1}\Delta S_{obs}(k)
\]

Entonces:

\[
S_t=H_t-\sum_{k=0}^{t-1}u_k
\]

Sustituyendo \(u_k\):

\[
S_t=H_t-\sum_{k=0}^{t-1}\sum_{\ell=0}^{L-1}a_\ell x_{k,\ell}
\]

Esto es lineal en las variables binarias, por lo que cualquier término cuadrático en \(S_t\) puede convertirse en QUBO.

## 11.4 Término QUBO para \(C_{dev}\)

\[
C_{dev}=\sum_{t=0}^{T-1}u_t^2
\]

con:

\[
u_t=\sum_\ell a_\ell x_{t,\ell}
\]

Agregar al QUBO:

\[
w_2\sum_t\left(\sum_\ell a_\ell x_{t,\ell}\right)^2
\]

Si se aplica one-hot, esto se simplifica a:

\[
w_2\sum_t\sum_\ell a_\ell^2x_{t,\ell}
\]

Pero para robustez puede implementarse con un helper general `add_square_of_linear_expression`.

## 11.5 Término QUBO para \(C_{smooth}\)

\[
C_{smooth}=\sum_{t=1}^{T-1}(u_t-u_{t-1})^2
\]

con:

\[
u_t-u_{t-1}=\sum_\ell a_\ell x_{t,\ell}-\sum_\ell a_\ell x_{t-1,\ell}
\]

Agregar:

\[
w_3\sum_{t=1}^{T-1}
\left(
\sum_\ell a_\ell x_{t,\ell}-\sum_\ell a_\ell x_{t-1,\ell}
\right)^2
\]

Esto genera acoplamientos entre semanas consecutivas.

## 11.6 Término \(C_{crit}\): dos opciones

El término oficial usa una función hinge:

\[
\max(0,S_{min}-S_t)^2
\]

Esto no es directamente cuadrático por el operador `max`. Se proponen dos rutas.

---

### Opción A — MVP preliminar, más simple, no estrictamente oficial

Usar un costo cuadrático suave alrededor del umbral:

\[
C_{storage}^{soft}=\sum_t(S_t-S_{target})^2
\]

con \(S_{target}=S_{min}\) o un valor operativo superior.

Ventaja:

- QUBO simple.
- Sin variables auxiliares.
- Útil para validar construcción de \(Q\), QAOA y pipeline.

Desventaja:

- Penaliza también estar por arriba de \(S_{min}\).
- No reproduce exactamente \(C_{crit}\).
- Debe reportarse como aproximación no oficial.

---

### Opción B — Formulación más fiel al benchmark con déficit y superávit

Introducir variables auxiliares no negativas:

\[
D_t\ge0
\]

\[
E_t\ge0
\]

con:

\[
S_t+D_t-E_t=S_{min}
\]

Interpretación:

- \(D_t\): déficit bajo \(S_{min}\).
- \(E_t\): superávit sobre \(S_{min}\).

Si la igualdad se cumple y se minimiza \(D_t^2\), entonces:

\[
D_t=\max(0,S_{min}-S_t)
\]

Agregar al QUBO:

\[
w_1\sum_tD_t^2+
P_{crit}\sum_t(S_t+D_t-E_t-S_{min})^2
\]

Codificar \(D_t\) y \(E_t\) con bits:

\[
D_t=q_S\sum_r2^rd_{t,r}
\]

\[
E_t=q_S\sum_r2^re_{t,r}
\]

donde:

- \(q_S\) es la resolución de almacenamiento en \(m^3\);
- el número de bits debe cubrir el rango esperado de déficit/superávit;
- esta formulación aumenta el número de qubits, pero es más fiel al benchmark.

Recomendación:

1. Implementar primero Opción A para validar.
2. Implementar Opción B para la versión final de QUBO/QAOA.
3. Comparar ambas y reportar diferencias.

---

## 11.7 Restricción de descarga no negativa

La restricción es:

\[
R_t=R_{obs,t}+u_t\ge0
\]

Como \(u_t\) toma niveles discretos, hay dos formas.

### Forma simple

Si para una semana \(t\) y nivel \(\ell\):

\[
R_{obs,t}+a_\ell<0
\]

entonces prohibir ese nivel con penalización lineal:

\[
P_R x_{t,\ell}
\]

### Forma con slack

Usar:

\[
R_{obs,t}+u_t-r_t=0
\]

con \(r_t\ge0\) codificado en binario. Es más costoso en qubits y no se recomienda en el MVP.

---

## 11.8 Restricción de balance acumulado de liberaciones

La restricción oficial es:

\[
\left|\sum_tu_t\right|\le B
\]

con:

\[
B=\eta\sum_tR_{obs,t}
\]

### Aproximación práctica QUBO

Agregar:

\[
P_{bal}\left(\sum_tu_t\right)^2
\]

Esto no implementa exactamente la desigualdad, pero favorece balance cercano a cero y evita retención sistemática.

### Formulación exacta con slacks

Transformar en dos desigualdades:

\[
\sum_tu_t\le B
\]

\[
-\sum_tu_t\le B
\]

Agregar slacks no negativos \(\eta_+\), \(\eta_-\):

\[
\sum_tu_t+\eta_+=B
\]

\[
-\sum_tu_t+\eta_-=B
\]

Penalizar:

\[
P_{bal}\left(\sum_tu_t+\eta_+-B\right)^2+
P_{bal}\left(-\sum_tu_t+\eta_--B\right)^2
\]

Codificar \(\eta_+\), \(\eta_-\) en binario.

Recomendación:

- MVP: usar \(P_{bal}(\sum u)^2\) y luego filtrar factibilidad al decodificar muestras.
- Versión final: slacks exactos si el número de qubits lo permite.

---

## 11.9 Restricción de almacenamiento físico

\[
0\le S_t\le S_{max}
\]

### MVP

1. No incluir en QUBO.
2. Después de obtener muestras, simular \(S_t\).
3. Rechazar muestras infeasibles.
4. Reportar tasa de factibilidad.

### Versión penalizada

Para \(S_t\le S_{max}\):

\[
S_t+s^{upper}_t=S_{max},\quad s^{upper}_t\ge0
\]

Para \(S_t\ge0\):

\[
-S_t+s^{lower}_t=0,\quad s^{lower}_t\ge0
\]

Agregar penalizaciones cuadráticas con slacks binarios. Esta versión requiere muchos qubits auxiliares.

---

# 12. Helper recomendado para construir QUBO

Implementar un helper general:

```python
def add_square_of_linear_expression(Q, offset, expr, weight):
    """
    Adds weight * (c + sum_i a_i x_i)^2 to Q.

    expr:
        dict with:
        - 'constant': c
        - 'linear': {var_index: coefficient}

    For binary x_i, x_i^2 = x_i.
    Expansion:
        weight * [
            c^2
            + sum_i (2*c*a_i + a_i^2) x_i
            + sum_{i<j} 2*a_i*a_j x_i*x_j
        ]
    """
```

Expansion:

\[
\left(c+\sum_i a_ix_i\right)^2
=c^2+\sum_i(2ca_i+a_i^2)x_i+2\sum_{i<j}a_ia_jx_ix_j
\]

Use this helper for:

- one-hot penalties;
- \(C_{dev}\);
- \(C_{smooth}\);
- soft storage cost;
- deficit/surplus equality;
- release-balance penalty;
- slack equality constraints.

---

# 13. QUBO objective recommended for implementation

## 13.1 MVP QUBO

Use:

\[
J_{MVP}=w_2C_{dev}+w_3C_{smooth}
+\lambda_S\sum_t(S_t-S_{target})^2
+P_{onehot}\sum_t\left(\sum_\ell x_{t,\ell}-1\right)^2
+P_{bal}\left(\sum_tu_t\right)^2
+P_R\sum_{t,\ell\in\mathcal I_R}x_{t,\ell}
\]

where \(\mathcal I_R\) is the set of invalid release choices with \(R_{obs,t}+a_\ell<0\).

This is the fastest route to a working QUBO/QAOA prototype.

## 13.2 Official-compatible QUBO

Use:

\[
J_{official-like}=w_1\sum_tD_t^2+w_2C_{dev}+w_3C_{smooth}
+P_{crit}\sum_t(S_t+D_t-E_t-S_{min})^2
+P_{onehot}\sum_t\left(\sum_\ell x_{t,\ell}-1\right)^2
+P_{bal}\cdot \text{balance penalties}
+P_R\cdot \text{release penalties}
+P_S\cdot \text{storage bound penalties}
\]

This should be labeled **official-compatible** only if:

- \(\Delta S_{obs}\) comes from `Discharge.Total.Change-in-Storage@08461200`;
- \(S_{max}\) is the official Falcon total conservation storage capacity;
- official weights are used;
- constraints are enforced or rejected consistently;
- all units are compatible.

---

# 14. Elección de coeficientes de penalización

Los notebooks de penalización advierten que \(P\) debe ser suficientemente grande para evitar violaciones, pero no tan grande que opaque diferencias entre soluciones factibles.

Estrategia práctica:

1. Calcular escala típica del objetivo:

\[
J_{scale}\approx |w_1C_{crit}^{hist}|+|w_2C_{dev}^{rule}|+|w_3C_{smooth}^{rule}|
\]

2. Probar:

```text
P_onehot = 10 * J_scale
P_bal    = 1  * J_scale / max(B^2, eps)
P_R      = 10 * J_scale
P_crit   = 10 * J_scale / storage_scale^2
```

3. Hacer barrido de sensibilidad:

```text
P in {0.1, 1, 10, 100} * base_penalty
```

4. Reportar:

- SRS;
- factibilidad;
- violaciones one-hot;
- violaciones de balance;
- almacenamiento mínimo;
- runtime.

---

# 15. Conversión QUBO → Ising / QAOA

Para QAOA, el QUBO puede convertirse a Ising usando:

\[
x_i=\frac{1-z_i}{2}
\]

con \(z_i\in\{-1,1\}\).

Recomendación técnica:

- Construir primero `Q` como matriz o diccionario triangular superior.
- Convertir a `QuadraticProgram` de Qiskit Optimization.
- Usar `QuadraticProgram.to_ising()` o el conversor de Qiskit correspondiente.
- Ejecutar QAOA con simulador para `T=12, L=3`.
- Comparar contra búsqueda exhaustiva y simulated annealing.

No empezar con hardware real.

---

# 16. Estructura de archivos esperada

Claude debe crear o mantener esta estructura:

```text
FalconChallenge/
│
├── data/
│   ├── source_csv/
│   │   ├── DataSetExport-Discharge.Best Available@08461300-...
│   │   ├── DataSetExport-Total Storage.Web-Daily-tcm@08461200-...
│   │   ├── DataSetExport-Reservoir Elevation.Web-Daily-m@08461200-...
│   │   └── DataSetExport-Lake Area.Best Available@08461200-...
│   └── processed/
│       └── falcon_weekly_benchmark.csv
│
├── src/
│   ├── data_loader.py
│   ├── eda_falcon.py
│   ├── storage_model.py
│   ├── srs_score.py
│   ├── classical_baselines.py
│   ├── qubo_utils.py
│   ├── qubo_falcon.py
│   ├── solve_classical.py
│   └── solve_qaoa.py
│
├── results/
│   ├── data_quality_report.md
│   ├── benchmark_summary.csv
│   ├── benchmark_timeseries.csv
│   ├── qubo_summary.md
│   ├── qaoa_results.csv
│   └── figures/
│
└── CLAUDE.md
```

---

# 17. Implementación mínima esperada

## 17.1 `src/data_loader.py`

Funciones:

```python
load_ibwc_export_csv(path) -> pd.DataFrame
identify_dataset(path) -> dict
build_weekly_benchmark(input_dir, output_path) -> pd.DataFrame
```

Debe manejar que los CSV tienen encabezado especial:

```text
#Data Set Export - ...
Timestamp (UTC-06:00), Value (...)
```

La primera fila real contiene nombres de columnas. El loader debe renombrar a:

```text
timestamp
value
```

## 17.2 `src/storage_model.py`

Funciones:

```python
simulate_storage(S0, deltaS_obs, u) -> np.ndarray
compute_release(R_obs, u) -> np.ndarray
check_constraints(R_obs, u, S, Smax, umax, eta) -> dict
```

## 17.3 `src/srs_score.py`

Funciones:

```python
compute_costs(S, u, Smin) -> dict
compute_weights(T, Smin, umax) -> dict
compute_srs(costs, weights) -> float
```

## 17.4 `src/classical_baselines.py`

Funciones:

```python
historical_policy(T) -> np.ndarray
threshold_policy(R_obs, deltaS_obs, S0, Smin, delta_u) -> np.ndarray
exhaustive_search_small_instance(levels, R_obs, deltaS_obs, S0, params) -> dict
```

## 17.5 `src/qubo_utils.py`

Funciones:

```python
add_linear(Q, i, coeff)
add_quadratic(Q, i, j, coeff)
add_square_of_linear_expression(Q, offset, expr, weight)
qubo_energy(Q, x, offset=0.0)
make_upper_triangular(Q)
```

## 17.6 `src/qubo_falcon.py`

Funciones:

```python
build_variable_index(T, L, aux_config=None) -> dict
linear_expr_u(t, levels, var_index) -> dict
linear_expr_storage(t, S0, deltaS_obs, levels, var_index) -> dict
build_falcon_qubo_mvp(data, params, penalties) -> dict
build_falcon_qubo_official_like(data, params, penalties, aux_config) -> dict
decode_solution(x, var_index, levels, data, params) -> dict
```

## 17.7 `src/solve_qaoa.py`

Debe quedar preparado para:

- convertir QUBO a Qiskit `QuadraticProgram`;
- correr QAOA para small instance;
- guardar muestras y soluciones;
- comparar con baselines.

Si Qiskit no está instalado, el script debe fallar con mensaje claro, no romper todo el pipeline.

---

# 18. Salidas esperadas

## 18.1 `results/data_quality_report.md`

Debe incluir:

- archivos usados;
- columnas detectadas;
- unidades detectadas;
- rango de fechas;
- faltantes;
- si existe o no \(\Delta S_{obs}\) oficial;
- si el benchmark es oficial o preliminar.

## 18.2 `results/benchmark_summary.csv`

Columnas:

```text
mode
method
T
L
SRS
DeltaSRS_vs_historical
DeltaSRS_vs_threshold
Ccrit
Cdev
Csmooth
min_storage
weeks_below_Smin
release_balance_error
feasible
runtime_seconds
```

## 18.3 `results/qubo_summary.md`

Debe incluir:

- número de variables binarias;
- número de variables auxiliares;
- número de términos lineales;
- número de términos cuadráticos;
- penalties usados;
- tipo de formulación: `MVP_SOFT_STORAGE` u `OFFICIAL_LIKE_WITH_DEFICIT_SLACKS`;
- advertencia si usa \(\Delta S_{obs}\) derivado.

## 18.4 Figuras

Guardar en `results/figures/`:

- `storage_comparison.png`
- `release_adjustments.png`
- `release_comparison.png`
- `srs_comparison.png`
- `qubo_energy_distribution.png` si hay muestras QAOA/annealing.

---

# 19. Orden recomendado de trabajo para Claude

1. Inspeccionar CSV anexos.
2. Construir `data/processed/falcon_weekly_benchmark.csv`.
3. Generar EDA y reporte de calidad.
4. Implementar simulación clásica.
5. Implementar SRS.
6. Implementar baseline histórico.
7. Implementar baseline de regla de umbral.
8. Implementar búsqueda exhaustiva para `T=12, L=3`.
9. Implementar QUBO MVP.
10. Validar QUBO comparando energía QUBO vs costo clásico decodificado.
11. Implementar QUBO official-like con déficit/superávit si el tamaño lo permite.
12. Implementar QAOA en simulador.
13. Comparar contra baselines.
14. Guardar reportes reproducibles.

---

# 20. Criterios de aceptación

Una implementación se considera aceptable si:

1. No modifica los datos originales.
2. Produce un dataset semanal procesado.
3. Reporta claramente si el dataset es oficial o preliminar.
4. Calcula SRS histórico y SRS de regla de umbral.
5. Implementa al menos una optimización clásica para `T=12, L=3`.
6. Construye una matriz QUBO o diccionario QUBO reproducible.
7. Decodifica soluciones binarias a secuencias \(u(t)\).
8. Verifica restricciones después de decodificar.
9. Compara resultados contra baselines.
10. Reporta runtime y escalabilidad.
11. Distingue claramente entre:
    - resultados oficiales;
    - resultados preliminares;
    - resultados sintéticos;
    - resultados aproximados por QUBO MVP.

---

# 21. Prompt operativo para Claude

Usa este prompt dentro de Claude/Codex si se quiere arrancar la implementación:

```text
Lee este CLAUDE.md completo antes de modificar archivos.

Primero implementa sólo las etapas 1 a 8 del orden recomendado:
1. inspección de CSV,
2. dataset semanal procesado,
3. EDA,
4. simulación clásica,
5. SRS,
6. baseline histórico,
7. baseline por umbral,
8. búsqueda exhaustiva para T=12, L=3.

Reglas:
- No modifiques datos originales.
- Guarda procesados en data/processed/.
- Guarda resultados en results/.
- No implementes QAOA todavía.
- Antes de QAOA, valida que el costo clásico y la energía QUBO coincidan para varias soluciones binarias decodificadas.

Cuando termines, entrega:
- data_quality_report.md,
- benchmark_summary.csv,
- benchmark_timeseries.csv,
- figuras principales,
- y una explicación clara de si el resultado es oficial o preliminar.
```

---

# 22. Notas críticas para evitar errores

1. La descarga `Discharge.Best Available@08461300` suele estar en \(m^3/s\). No se puede restar directamente de almacenamiento en \(m^3\). Convertir a volumen semanal.
2. El archivo de `Total Storage` puede tener nombre con `tcm`, pero el CSV anexo indica `Value (m^3)`. El loader debe leer la unidad real del encabezado.
3. Si \(\Delta S_{obs}\) se deriva de almacenamiento, etiquetar el resultado como preliminar.
4. Si se usa \(S_{max}=\max(S_{obs})\) por falta de capacidad oficial, etiquetar como preliminar.
5. El QUBO MVP con almacenamiento suave no reproduce exactamente \(C_{crit}\). Debe reportarse como aproximación.
6. El QUBO con slacks aumenta mucho el número de qubits. Usarlo primero en `T=12`.
7. Siempre decodificar y verificar factibilidad después de resolver QUBO/QAOA.
8. No comparar QAOA contra histórico solamente; también comparar contra regla de umbral y búsqueda exhaustiva/simulated annealing cuando sea posible.

---

# 23. Resumen ejecutivo para implementación

La ruta más factible es:

```text
CSV anexos
→ limpieza semanal
→ benchmark clásico
→ QUBO one-hot L=3, T=12
→ validación contra búsqueda exhaustiva
→ QUBO L=5, T=26
→ QAOA / simulated annealing
→ comparación ΔSRS
```

La formulación QUBO recomendada para iniciar es one-hot:

\[
x_{t,\ell}=1 \iff u_t=a_\ell
\]

con:

\[
u_t=\sum_\ell a_\ell x_{t,\ell}
\]

\[
S_t=H_t-\sum_{k<t}\sum_\ell a_\ell x_{k,\ell}
\]

Minimizar:

\[
J=w_1C_{crit}+w_2C_{dev}+w_3C_{smooth}+\text{penalties}
\]

Primero implementar `Ccrit` con aproximación suave para validar. Después implementar déficit/superávit con slacks para acercarse al benchmark oficial.
