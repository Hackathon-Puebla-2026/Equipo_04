"""Configuracion del builder QUBO (FalconConfig).

Un solo dataclass que hace configurable al builder: instancia, encoding, que
terminos incluir y como tratar cada restriccion. Los defaults son los DECIDIDOS
segun el analisis de datos (docs/SPEC_IMPLEMENTACION_QUBO.md §4,
docs/ANALISIS_DP_Y_RESULTADOS.md §8):

- encoding one-hot para el MVP.
- storage_bounds = "drop": la cota 0<=S<=S_max nunca ata en estos datos.
- balance = "slack" (exacto de 1 slack, ~log2 qubits) preferido; "soft" es atajo MVP.
- c_crit = "soft" (S_target=S_min): EXACTO aqui porque el storage siempre < S_min.
- release_nonneg = "prohibit" (nivel prohibido, sin slacks). onehot = "penalty".
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class FalconConfig:
    # Instancia
    T: int
    L: int

    # Encoding de la variable de decision
    encoding: str = "onehot"            # "onehot" | "domainwall" | "binary"

    # Terminos de costo a incluir
    c_dev: bool = True
    c_smooth: bool = True
    c_crit: str = "soft"               # "soft" (S_target=S_min) | "deficit_slack" (Opcion B)

    # Tratamiento de restricciones
    onehot: str = "penalty"            # "penalty" | "xy_mixer"
    balance: str = "slack"             # "slack" (exacto 1-slack, preferido) | "soft" | "off"
    release_nonneg: str = "prohibit"   # "prohibit" | "slack" | "off"
    storage_bounds: str = "drop"       # "drop" | "postselect" | "slack"

    # Pesos de penalizacion (None -> auto-escala desde J_scale en el builder)
    penalties: dict = field(default_factory=dict)

    # Normalizacion interna del Q para condicionamiento (m³ deja pesos ~1e-18)
    normalize: str | None = "maxabs"   # "maxabs" | "delta_u2" | None

    def __post_init__(self):
        if self.L % 2 == 0:
            raise ValueError(f"L debe ser impar (1-of-L simetrico); L={self.L}")
        valid = {
            "encoding": {"onehot", "domainwall", "binary"},
            "c_crit": {"soft", "deficit_slack"},
            "onehot": {"penalty", "xy_mixer"},
            "balance": {"slack", "soft", "off"},
            "release_nonneg": {"prohibit", "slack", "off"},
            "storage_bounds": {"drop", "postselect", "slack"},
        }
        for attr, allowed in valid.items():
            val = getattr(self, attr)
            if val not in allowed:
                raise ValueError(f"{attr}={val!r} invalido; opciones: {sorted(allowed)}")
        if self.normalize not in (None, "maxabs", "delta_u2"):
            raise ValueError(f"normalize={self.normalize!r} invalido")
