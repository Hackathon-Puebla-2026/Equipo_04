# Guía de trabajo del repositorio (Equipo 4)

Convenciones para que los tres integrantes (**julian**, **emilio**, **ivan**) trabajemos en paralelo sin pisarnos ni generar conflictos de git.

## 1. `FalconChallenge/` está congelado

La carpeta `FalconChallenge/` (el PDF del spec) **no se modifica por nadie**. Es común a todos; cualquier cambio genera conflictos de merge para todo el equipo.

- Usá `docs/FalconChallenge_V6.md` (transcripción en markdown) como referencia de trabajo.
- El PDF sigue siendo la fuente de verdad ante cualquier diferencia.

## 2. Estructura de carpetas y propiedad

```
scripts/            # Python compartido en la raíz (lo importa todo el equipo)
  julian/           # scripts personales / experimentales
  emilio/
  ivan/
notebooks/          # notebooks, una subcarpeta por integrante
  julian/
  emilio/
  ivan/
results/            # salidas de corridas; se commitea todo
  julian/
  emilio/
  ivan/
docs/               # referencia compartida (spec transcripto, esta guía)
FalconChallenge/    # CONGELADO — no tocar
data/               # ignorado por git (datasets)
```

**Solo editá tu propia subcarpeta** (`julian/` / `emilio/` / `ivan/`). No toques la carpeta de un compañero.

## 3. Código compartido en la raíz

El código reutilizable (cálculo del **SRS**, carga de datos, **constantes oficiales del benchmark**, codificación QUBO/Ising) vive en la **raíz de `scripts/`** para que todos importen lo mismo y el benchmark sea idéntico entre integrantes.

- Mantené el SRS y las constantes oficiales en un solo lugar; no hagas copias divergentes.
- Cambiar código compartido afecta a todos: **coordiná antes** de modificarlo.

## 4. `results/`: se commitea todo

Todas las salidas de corridas van versionadas. Nombrá cada corrida de forma identificable, incluyendo dueño, tamaño de instancia y parámetros clave:

```
results/julian/medium_T26_L5_qaoa_2026-06-29.json
```

## 5. `data/` ignorado por git

Los datasets **no se commitean** (ya está en `.gitignore`). Se comparten por la carpeta de SharePoint del benchmark (URL en `docs/FalconChallenge_V6.md`, sección 9).

## 6. Higiene de notebooks

Antes de commitear, *restart & run-all* para que las salidas guardadas coincidan con el código. Mantiene los diffs coherentes aun dentro de tu propia carpeta.

## 7. Reproducibilidad

Anotá el entorno para que las corridas se reproduzcan. `requirements.txt` en la raíz tiene las librerías base; las variantes GPU (`lightning.gpu`, Aer-GPU, `cudaq`) son opcionales y dependen de la máquina (ver fallback CPU en `CLAUDE.md`).

## 8. Idioma

Español para prosa, comentarios y mensajes de commit, consistente con el resto del repo.

## 9. Cómo programar con Claude Code

Este repo está armado para trabajo **guiado por el spec**. Estas prácticas hacen que Claude sea mucho más preciso. No hace falta seguirlas al pie de la letra; sumá las que puedas.

### 9.1. Jerarquía: fuente de verdad vs guía

Claude trata por defecto todo lo que tiene en contexto como igual de autoritativo. Por eso `CLAUDE.md` deja explícita la jerarquía:

- **Fuente de verdad:** el spec (`FalconChallenge_V6.pdf` + `docs/FalconChallenge_V6.md`) y las constantes oficiales del benchmark en `CLAUDE.md`. **Ante cualquier conflicto, gana el spec.**
- **Guía / ejemplos:** los docs de QUBO, el handoff de Georgia, los hallazgos del equipo. Se adaptan, no se obedecen.

Cuando crees un doc nuevo, marcá arriba si es *autoritativo* o *guía/ejemplo*. Así Claude resuelve conflictos bien en vez de copiar una implementación de referencia como si fuera el spec.

### 9.2. Presupuesto de contexto: auto-import vs on-demand

`CLAUDE.md` importa con `@docs/...` solo lo que se necesita en **cada** turno (spec, esta guía, los hallazgos clave). Los docs grandes (guías de implementación, snippets) quedan **on-demand**: se leen solo cuando hacen falta.

- Si un doc está auto-importado pero rara vez cambia lo que hace Claude, **bajalo a on-demand**.
- Contexto gastado en un doc que no se usa es contexto que no se usa para razonar.

### 9.3. Hacé el spec verificable (lo más importante)

Un spec sirve más cuando se puede chequear. Ejemplo ya hecho: `Smax_search.ipynb` deriva `S_max` por dos caminos independientes y coinciden exactamente. Eso convierte un número del spec en un **invariante verificable**.

Generalizá esto:
- Poné las constantes oficiales (`L`, `w1/w2/w3`, `η`, `S_min = 0.25·S_max`) como constantes con nombre en `scripts/` (raíz, compartido), citando la ecuación del spec en un comentario.
- Escribí asserts chiquitos: `u=0` debe reproducir `SRS_hist`; la regla umbral debe dar `u ∈ {-Δu, 0}`; la factibilidad (`R(t)≥0`, `|Σu| ≤ η·ΣR_obs`) debe cumplirse. Así "¿coincide con el spec?" es correr un test, no una opinión.

### 9.4. El loop de trabajo

1. **Apuntá a la sección concreta del spec**, no "leé el PDF". Ej.: *"Implementá `C_crit` según `docs/FalconChallenge_V6.md`, usando las constantes de `falcon_reservoir_constants.json`."* Referencias precisas > "arreglátelas".
2. **Usá plan mode (`/plan`)** para cualquier cosa no trivial. Saca a la luz las decisiones que son **tuyas** (encoding, dónde va el código) antes de escribir nada.
3. **Decí la verificación de entrada:** *"está bien cuando `u=0` reproduce el SRS histórico."* Claude construye hacia un objetivo chequeable, no hacia algo que *parece* bien.
4. **Respetá los límites congelados:** `FalconChallenge/` es read-only y no se tocan carpetas de otros (secciones 1 y 2). Está en `CLAUDE.md`, así que Claude lo respeta.

### 9.5. Que el spec y el código no se desincronicen

El riesgo del trabajo guiado por spec es la deriva: el código dice una cosa, el spec otra, nadie lo nota.

- Cuando se resuelve una **interpretación** del spec (ej.: la discrepancia `S_max` preliminar vs oficial), escribilo en un digest de decisiones (como `docs/HALLAZGOS_CLAVE.md`).
- Si Claude propone algo que **se desvía** del spec (otra normalización, otros pesos), debe marcarlo como *"alternativa, reportar por separado"* - que es justo lo que pide el spec. Exigíselo.

### 9.6. A vigilar

- **Docs de referencia con rutas equivocadas:** el handoff de Georgia asume `FalconChallenge/src/` y `data/processed/` que no son nuestro layout. `CLAUDE.md` ya lo aclara; la sección 2 de esta guía manda para dónde van los archivos.
- **No sobre-importar:** a medida que sumemos docs de hallazgos, la lista de auto-import crece. Podá: si un digest reemplaza a un doc detallado, pasá el detallado a on-demand.
