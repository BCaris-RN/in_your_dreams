# ============================================================================
# Copyright (c) 2026 Brandon W. Caris. All rights reserved.
# Part of the Coupled Marine Engine Suite (in_your_dreams/Coupled_Engine)
# Licensed under CC BY 4.0 with Explicit Academic Citation Mandate.
# See master LICENSE file at repository root for full terms and string.
# ============================================================================
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
