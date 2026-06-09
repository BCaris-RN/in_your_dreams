"""Compatibility exports for the model's core equations."""

from .solubility import (
    effective_dissolution_enthalpy,
    van_t_hoff_henry,
    warner_weiss_henry,
)

__all__ = [
    "effective_dissolution_enthalpy",
    "van_t_hoff_henry",
    "warner_weiss_henry",
]
