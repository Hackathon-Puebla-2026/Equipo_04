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
