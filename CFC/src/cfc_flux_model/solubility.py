# ============================================================================
# Copyright (c) 2026 Brandon W. Caris. All rights reserved.
# Part of the Coupled Marine Engine Suite (in_your_dreams/Coupled_Engine)
# Licensed under CC BY 4.0 with Explicit Academic Citation Mandate.
# See master LICENSE file at repository root for full terms and string.
# ============================================================================
"""Henry-solubility relationships for CFC-11 and CFC-12."""

from math import exp, log

from .parameters import Species

GAS_CONSTANT_J_MOL_K = 8.31446261815324
STANDARD_ATMOSPHERE_PA = 101_325.0


def warner_weiss_henry(
    species: Species,
    temperature_k: float,
    salinity_psu: float,
    seawater_density_kg_m3: float = 1025.0,
) -> float:
    """Return Hcp in mol m-3 Pa-1 from Warner and Weiss (1985).

    Their fitted coefficient is in mol kg-1 atm-1 and is valid from 273 K to
    313 K and salinity 0 to 40. Density converts kg of seawater to m3.
    """
    if not 273.0 <= temperature_k <= 313.0:
        raise ValueError("Warner-Weiss temperature must be between 273 and 313 K")
    if not 0.0 <= salinity_psu <= 40.0:
        raise ValueError("salinity must be between 0 and 40")
    if seawater_density_kg_m3 <= 0:
        raise ValueError("seawater density must be positive")

    c = species.coefficients
    scaled_temperature = temperature_k / 100.0
    log_k_mol_kg_atm = (
        c.a1
        + c.a2 * 100.0 / temperature_k
        + c.a3 * log(scaled_temperature)
        + salinity_psu
        * (c.b1 + c.b2 * scaled_temperature + c.b3 * scaled_temperature**2)
    )
    k_mol_kg_atm = exp(log_k_mol_kg_atm)
    return k_mol_kg_atm * seawater_density_kg_m3 / STANDARD_ATMOSPHERE_PA


def effective_dissolution_enthalpy(
    species: Species,
    reference_temperature_k: float,
    salinity_psu: float,
) -> float:
    """Return the local Van 't Hoff dissolution enthalpy in J mol-1.

    The derivative is taken analytically from the Warner-Weiss fit. This lets
    the experiment use a clearly defined Van 't Hoff perturbation while
    retaining a measured seawater reference solubility.
    """
    if reference_temperature_k <= 0:
        raise ValueError("reference temperature must be positive")

    c = species.coefficients
    t = reference_temperature_k
    derivative_wrt_inverse_t = (
        100.0 * c.a2
        - c.a3 * t
        - salinity_psu * (c.b2 * t**2 / 100.0 + 2.0 * c.b3 * t**3 / 10_000.0)
    )
    return -GAS_CONSTANT_J_MOL_K * derivative_wrt_inverse_t


def van_t_hoff_henry(
    reference_henry: float,
    temperature_k: float,
    reference_temperature_k: float,
    dissolution_enthalpy_j_mol: float,
) -> float:
    """Return temperature-adjusted Hcp in the same units as the reference.

    ``reference_henry`` is the solubility-form coefficient Hcp = Caq / p.
    A negative dissolution enthalpy therefore makes Hcp decrease as water
    warms.
    """
    if reference_henry <= 0:
        raise ValueError("reference_henry must be positive")
    if temperature_k <= 0 or reference_temperature_k <= 0:
        raise ValueError("temperatures must be positive Kelvin values")

    exponent = -(dissolution_enthalpy_j_mol / GAS_CONSTANT_J_MOL_K) * (
        (1.0 / temperature_k) - (1.0 / reference_temperature_k)
    )
    return reference_henry * exp(exponent)


def species_reference_properties(
    species: Species,
    reference_temperature_k: float,
    salinity_psu: float,
    seawater_density_kg_m3: float,
) -> tuple[float, float]:
    """Return reference Hcp and local effective dissolution enthalpy."""
    return (
        warner_weiss_henry(
            species,
            reference_temperature_k,
            salinity_psu,
            seawater_density_kg_m3,
        ),
        effective_dissolution_enthalpy(
            species,
            reference_temperature_k,
            salinity_psu,
        ),
    )
