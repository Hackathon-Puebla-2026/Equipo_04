# Guía de rendimiento: MacBook Pro M4 (CPU-first)

Cómo correr las simulaciones de este proyecto de forma performante en la laptop. **Objetivo primario =
CPU en el M4.** El servidor GPU (WCentroid / T4) tiene problemas de disponibilidad, así que **no
dependemos de GPU**. Basado en las particularidades del código del repo Georgia
(`docs/georgia_qubo_snippets.md`).

## 1. Backend cuántico (elegir el correcto)

- **Preferido en M4: PennyLane `lightning.qubit`** (statevector en C++, wheels nativas ARM). Rápido y
  se instala sin problemas en Apple Silicon.
- **`qiskit-aer`** (statevector C++) también sirve, pero en Apple Silicon + Python 3.14 puede no haber
  wheels -> usar un **venv con Python 3.11/3.12**. La ruta rápida de Georgia depende de Aer.
- **Nada de GPU**: no usar `qml.device("lightning.gpu")`, `AerSimulator(device="GPU")` ni cuStateVec.
- Poner el device **detrás de un flag/config** para que el mismo código corra CPU aquí y (si algún día
  está disponible) GPU en WCentroid. Default en M4: `lightning.qubit`.

## 2. Threads (evitar oversubscription)

- **Statevector**: usar todos los cores. Con Aer: `max_parallel_threads=n_cores`,
  `statevector_parallel_threshold=12` (paraleliza solo desde ~12 qubits; abajo el overhead de hilos
  perjudica).
- **MPS o varios restarts en paralelo (multiprocessing)**: **1 hilo por proceso** para no sobre-suscribir
  (`cores × hilos` explota). Georgia fija `OMP_NUM_THREADS=1` en MPS.
- Snippet de env vars (ajustar al conteo de cores del M4):

```python
import multiprocessing, os
n_cores = multiprocessing.cpu_count()
for v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(v, str(n_cores))   # statevector single-run
# Para restarts en paralelo: fijar estas a "1" DENTRO de cada worker.
```

## 3. Techo de memoria (lo que decide el tamaño de instancia)

Statevector denso usa `2ⁿ × 16 bytes`; el arreglo `precompute_diagonal` usa `2ⁿ × 8 bytes`. En una
laptop (RAM ~24-48 GB) el tope práctico es **~28-30 qubits**.

| n qubits | statevector | diagonal | ¿entra en M4? |
|---:|---:|---:|:--:|
| 24 | 256 MB | 128 MB | sí (holgado) |
| 28 | 4 GB | 2 GB | sí |
| 30 | 16 GB | 8 GB | límite |
| 36 | 1 TB | 0.5 TB | no |

Qubits por instancia/encoding (bits de decisión, ver `docs/ANALISIS_DP_Y_RESULTADOS.md` §8):

| Instancia | one-hot `T·L` | binary | domain-wall | ¿statevector en M4? |
|---|---:|---:|---:|:--:|
| debug T5/L3 | 15 | 10 | 10 | sí (cualquiera) |
| small T12/L3 | 36 | **24** | **24** | solo compacto (binary/domain-wall) |
| medium T26/L5 | 130 | 78 | 104 | **no** (usar MPS/sampling o chunking) |
| large T52 | 260-364 | 156 | 208-312 | **no** |

**Conclusión:** en M4, statevector exacto solo para instancias chicas con **encoding compacto**
(small binary = 24q). Para medium/large: **MPS**, **sampling (shots)** o **chunking temporal (E2)**.

## 4. Trucos de rendimiento a usar

- **Encodings compactos** (binary/domain-wall) no son solo "más lindos": son lo que hace que la
  instancia chica entre en statevector.
- **Chunking temporal (E2)**: partir T en n sub-ventanas -> cada QUBO usa `(T/n)·L` qubits y entra.
- **`precompute_diagonal`** (energías `2ⁿ` una vez, luego `dot(|ψ|², diag)` por iteración de COBYLA):
  gran aceleración, pero **solo hasta n≈26** por memoria. Más allá, MPS/shots.
- **Circuito QAOA manual** (H/Rz/CX-Rz-CX/Rx) en vez de `PauliEvolutionGate` (que cuelga la simulación).
  Ver `docs/georgia_qubo_snippets.md`.
- **Clásico vectorizado**: evaluar energías QUBO en lote con `X @ Q` sobre arreglos `(batch, n)` en vez
  de un loop Python por bitstring (~100x más rápido). Aplicar en brute/exhaustive.
- **Warmup**: siempre una corrida de calentamiento antes de medir tiempos (como el notebook de smoke).

## 5. Checklist de setup del entorno

1. `.venv` con `pandas numpy pytest` (ya hecho) para Fase 0/1 (numpy puro, sin backend cuántico).
2. Para Fase 3 (QAOA): crear venv con Python 3.11/3.12 e instalar `pennylane pennylane-lightning`
   (y opcional `qiskit qiskit-aer` si hay wheels). Verificar `import pennylane; qml.device("lightning.qubit", wires=4)`.
3. `dimod`/`neal` (SA de D-Wave) es opcional; SA propio en numpy también sirve y evita dependencia.
4. Nunca instalar variantes GPU (`pennylane-lightning-gpu`, `qiskit-aer-gpu`, `cudaq`) en la laptop.
