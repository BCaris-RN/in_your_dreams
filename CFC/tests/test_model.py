# ============================================================================
# Copyright (c) 2026 Brandon W. Caris. All rights reserved.
# Part of the Coupled Marine Engine Suite (in_your_dreams/Coupled_Engine)
# Licensed under CC BY 4.0 with Explicit Academic Citation Mandate.
# See master LICENSE file at repository root for full terms and string.
# ============================================================================
from math import isclose

import pytest

from cfc_flux_model import (
    CFC11,
    CFC12,
    effective_dissolution_enthalpy,
    van_t_hoff_henry,
    warner_weiss_henry,
)


def test_reference_temperature_returns_reference_coefficient() -> None:
    result = van_t_hoff_henry(
        reference_henry=0.42,
        temperature_k=298.15,
        reference_temperature_k=298.15,
        dissolution_enthalpy_j_mol=-20_000.0,
    )

    assert isclose(result, 0.42)


@pytest.mark.parametrize("temperature_k", [0.0, -1.0])
def test_nonphysical_temperature_is_rejected(temperature_k: float) -> None:
    with pytest.raises(ValueError, match="positive Kelvin"):
        van_t_hoff_henry(
            reference_henry=0.42,
            temperature_k=temperature_k,
            reference_temperature_k=298.15,
            dissolution_enthalpy_j_mol=-20_000.0,
        )


@pytest.mark.parametrize("species", [CFC11, CFC12])
def test_warming_reduces_cfc_solubility(species) -> None:
    cold = warner_weiss_henry(species, 283.15, 35.0)
    warm = warner_weiss_henry(species, 293.15, 35.0)

    assert warm < cold


@pytest.mark.parametrize("species", [CFC11, CFC12])
def test_effective_dissolution_enthalpy_is_exothermic(species) -> None:
    enthalpy = effective_dissolution_enthalpy(species, 288.15, 35.0)

    assert -40_000.0 < enthalpy < -20_000.0
