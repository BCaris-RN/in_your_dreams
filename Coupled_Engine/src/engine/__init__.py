# ============================================================================
# Copyright (c) 2026 Brandon W. Caris. All rights reserved.
# Part of the Coupled Marine Engine Suite (in_your_dreams/Coupled_Engine)
# Licensed under the MIT Open-Source License.
# ============================================================================
"""Unified marine biogeochemical box model for Paper 3.

The model combines stiff gas-phase halogen chemistry with time-dependent
CFC-11/CFC-12 air-sea exchange. All derivatives use SI seconds as the time
base. Gas-phase chemical concentrations use molecule cm-3, dissolved CHBr3
uses mol m-3, and CFC reservoirs use mol.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum
from math import exp, log
from typing import Any, Iterable

import numpy as np
from scipy.integrate import solve_ivp

AVOGADRO_MOL = 6.02214076e23
GAS_CONSTANT_J_MOL_K = 8.31446261815324
STANDARD_ATMOSPHERE_PA = 101_325.0
SECONDS_PER_YEAR = 365.25 * 86_400.0


class StateIndex(IntEnum):
    """Positions in the unified 17-element state vector."""

    CHBR3_GAS = 0
    BR = 1
    BRO = 2
    HO2 = 3
    OH = 4
    HOBR = 5
    CHBR3_AQUEOUS = 6

    CFC11_AIR = 7
    CFC11_MIXED = 8
    CFC11_DEEP = 9
    CFC11_CUMULATIVE_AIR_TO_SEA = 10
    CFC11_CUMULATIVE_EXTERNAL = 11

    CFC12_AIR = 12
    CFC12_MIXED = 13
    CFC12_DEEP = 14
    CFC12_CUMULATIVE_AIR_TO_SEA = 15
    CFC12_CUMULATIVE_EXTERNAL = 16


STATE_NAMES = (
    "CHBr3_gas_molecule_cm3",
    "Br_molecule_cm3",
    "BrO_molecule_cm3",
    "HO2_molecule_cm3",
    "OH_molecule_cm3",
    "HOBr_molecule_cm3",
    "CHBr3_aqueous_mol_m3",
    "CFC11_air_mol",
    "CFC11_mixed_mol",
    "CFC11_deep_mol",
    "CFC11_cumulative_air_to_sea_mol",
    "CFC11_cumulative_external_mol",
    "CFC12_air_mol",
    "CFC12_mixed_mol",
    "CFC12_deep_mol",
    "CFC12_cumulative_air_to_sea_mol",
    "CFC12_cumulative_external_mol",
)
STATE_SIZE = len(STATE_NAMES)


@dataclass(frozen=True)
class WarnerWeissCoefficients:
    """Warner-Weiss coefficients in mol kg-1 atm-1."""

    a1: float
    a2: float
    a3: float
    b1: float
    b2: float
    b3: float


@dataclass(frozen=True)
class TraceGas:
    """CFC properties needed by the physical exchange model."""

    name: str
    molecular_weight_g_mol: float
    coefficients: WarnerWeissCoefficients
    air_index: StateIndex
    mixed_index: StateIndex
    deep_index: StateIndex
    cumulative_flux_index: StateIndex
    cumulative_external_index: StateIndex


CFC11 = TraceGas(
    name="CFC-11",
    molecular_weight_g_mol=137.368,
    coefficients=WarnerWeissCoefficients(
        a1=-136.2685,
        a2=206.1150,
        a3=57.2805,
        b1=-0.148598,
        b2=0.095114,
        b3=-0.0163396,
    ),
    air_index=StateIndex.CFC11_AIR,
    mixed_index=StateIndex.CFC11_MIXED,
    deep_index=StateIndex.CFC11_DEEP,
    cumulative_flux_index=StateIndex.CFC11_CUMULATIVE_AIR_TO_SEA,
    cumulative_external_index=StateIndex.CFC11_CUMULATIVE_EXTERNAL,
)

CFC12 = TraceGas(
    name="CFC-12",
    molecular_weight_g_mol=120.913,
    coefficients=WarnerWeissCoefficients(
        a1=-124.4395,
        a2=185.4299,
        a3=51.6383,
        b1=-0.149779,
        b2=0.094668,
        b3=-0.0160043,
    ),
    air_index=StateIndex.CFC12_AIR,
    mixed_index=StateIndex.CFC12_MIXED,
    deep_index=StateIndex.CFC12_DEEP,
    cumulative_flux_index=StateIndex.CFC12_CUMULATIVE_AIR_TO_SEA,
    cumulative_external_index=StateIndex.CFC12_CUMULATIVE_EXTERNAL,
)

CFC_SPECIES = (CFC11, CFC12)


def _as_1d_float(name: str, values: Iterable[float]) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    if not np.isfinite(array).all():
        raise ValueError(f"{name} must contain only finite values")
    return array


@dataclass(frozen=True)
class CoupledForcing:
    """Time-varying physical and biogeochemical forcing arrays."""

    time_s: np.ndarray
    sst_k: np.ndarray
    air_temperature_k: np.ndarray
    wind_speed_m_s: np.ndarray
    salinity_psu: np.ndarray
    gradient_richardson_number: np.ndarray
    ph: np.ndarray
    chlorophyll_a_mg_m3: np.ndarray
    cfc11_ppt: np.ndarray
    cfc12_ppt: np.ndarray

    def __post_init__(self) -> None:
        fields = (
            "time_s",
            "sst_k",
            "air_temperature_k",
            "wind_speed_m_s",
            "salinity_psu",
            "gradient_richardson_number",
            "ph",
            "chlorophyll_a_mg_m3",
            "cfc11_ppt",
            "cfc12_ppt",
        )
        arrays: dict[str, np.ndarray] = {}
        for name in fields:
            arrays[name] = _as_1d_float(name, getattr(self, name))
            object.__setattr__(self, name, arrays[name])

        size = arrays["time_s"].size
        if size < 2:
            raise ValueError("forcing requires at least two time points")
        if any(array.size != size for array in arrays.values()):
            raise ValueError("all forcing arrays must have the same length")
        if np.any(np.diff(arrays["time_s"]) <= 0.0):
            raise ValueError("time_s must be strictly increasing")
        if np.any((arrays["sst_k"] < 273.0) | (arrays["sst_k"] > 313.0)):
            raise ValueError("sst_k must remain within 273-313 K")
        if np.any(arrays["air_temperature_k"] <= 0.0):
            raise ValueError("air_temperature_k must be positive")
        if np.any(arrays["wind_speed_m_s"] < 0.0):
            raise ValueError("wind_speed_m_s cannot be negative")
        if np.any((arrays["salinity_psu"] < 0.0) | (arrays["salinity_psu"] > 40.0)):
            raise ValueError("salinity_psu must remain within 0-40")
        if np.any((arrays["ph"] < 0.0) | (arrays["ph"] > 14.0)):
            raise ValueError("ph must remain within 0-14")
        if np.any(arrays["chlorophyll_a_mg_m3"] < 0.0):
            raise ValueError("chlorophyll_a_mg_m3 cannot be negative")
        if np.any(arrays["cfc11_ppt"] < 0.0) or np.any(arrays["cfc12_ppt"] < 0.0):
            raise ValueError("CFC atmospheric histories cannot be negative")

    def value(self, name: str, time_s: float) -> float:
        """Linearly interpolate one forcing field at ``time_s``."""
        return float(np.interp(time_s, self.time_s, getattr(self, name)))


@dataclass(frozen=True)
class UnifiedConfig:
    """Geometry, kinetics, and numerical settings for the coupled model."""

    ocean_area_m2: float = 3.61e14
    mixed_layer_depth_m: float = 50.0
    deep_ocean_depth_m: float = 3650.0
    marine_boundary_layer_height_m: float = 1000.0
    seawater_density_kg_m3: float = 1025.0
    surface_pressure_pa: float = STANDARD_ATMOSPHERE_PA
    tropospheric_dry_air_mol: float = 1.51e20
    vertical_exchange_timescale_s: float = 10.0 * SECONDS_PER_YEAR
    atmospheric_nudging_timescale_s: float = SECONDS_PER_YEAR / 12.0
    initial_deep_saturation_fraction: float = 0.15

    wanninkhof_coefficient_cm_hr: float = 0.251
    reference_schmidt_number: float = 660.0
    richardson_scaling_gamma: float = 1.0
    richardson_scaling_exponent: float = 0.35

    chbr3_henry_reference_mol_m3_pa: float = 1.9e-2
    chbr3_henry_reference_temperature_k: float = 298.15
    chbr3_henry_temperature_parameter_k: float = 4700.0
    chbr3_hydrolysis_arrhenius_prefactor_l_mol_s: float = 1.23e17
    chbr3_hydrolysis_activation_energy_j_mol: float = 107_300.0
    chbr3_hydrolysis_reference_ph: float = 8.1
    chlorophyll_emission_coefficient_molecule_cm2_s_per_mg_m3: float = 1.127e5
    chbr3_species_fraction: float = 2.0
    extra_tropical_coastal_modifier: float = 2.5
    initial_chbr3_aqueous_mol_m3: float = 5.0e-9

    chbr3_photolysis_s: float = 2.0e-6
    hobr_photolysis_s: float = 5.0e-3
    hox_production_molecule_cm3_s: float = 1.0e6
    chlorine_molecule_cm3: float = 1.0e4
    ozone_molecule_cm3: float = 7.0e11
    br_o3_rate_cm3_molecule_s: float = 1.2e-12
    bro_ho2_activation_parameter_k: float = 580.0
    ho2_background_loss_s: float = 5.0e-2
    oh_background_loss_s: float = 4.0
    hobr_termination_s: float = 1.0e-4

    initial_chbr3_gas_molecule_cm3: float = 1.0e9
    initial_br_molecule_cm3: float = 0.0
    initial_bro_molecule_cm3: float = 0.0
    initial_ho2_molecule_cm3: float = 1.0e8
    initial_oh_molecule_cm3: float = 1.0e6
    initial_hobr_molecule_cm3: float = 0.0

    relative_tolerance: float = 1.0e-7
    maximum_step_s: float = 86_400.0

    def __post_init__(self) -> None:
        if self.richardson_scaling_gamma < 0.0:
            raise ValueError("richardson_scaling_gamma cannot be negative")
        if self.richardson_scaling_exponent != 0.35:
            raise ValueError("richardson_scaling_exponent must be 0.35")

    @property
    def mixed_layer_volume_m3(self) -> float:
        return self.ocean_area_m2 * self.mixed_layer_depth_m

    @property
    def deep_ocean_volume_m3(self) -> float:
        return self.ocean_area_m2 * self.deep_ocean_depth_m


@dataclass(frozen=True)
class ModelSwitches:
    """Enable or disable the three Paper 3 dynamic indices."""

    ph_hydrolysis: bool = True
    chlorophyll_source: bool = True
    haline_stratification: bool = True


@dataclass(frozen=True)
class UnifiedResult:
    """Full trajectories and the two requested diagnostic matrices."""

    time_s: np.ndarray
    state: np.ndarray
    no_stratification_state: np.ndarray
    no_ph_biology_state: np.ndarray
    cfc_inversion_discrepancy_matrix: np.ndarray
    methane_lifetime_distortion_matrix: np.ndarray

    cfc_inversion_columns = (
        "time_s",
        "gradient_richardson_number",
        "kw_stratification_multiplier",
        "CFC11_discrepancy_Gg",
        "CFC12_discrepancy_Gg",
    )
    methane_lifetime_columns = (
        "time_s",
        "ph",
        "chlorophyll_a_mg_m3",
        "OH_dynamic_molecule_cm3",
        "OH_no_ph_biology_molecule_cm3",
        "methane_lifetime_distortion_percent",
    )


def warner_weiss_henry(
    species: TraceGas,
    temperature_k: float,
    salinity_psu: float,
    seawater_density_kg_m3: float,
) -> float:
    """Return the CFC solubility-form Henry coefficient in mol m-3 Pa-1."""
    if not 273.0 <= temperature_k <= 313.0:
        raise ValueError("Warner-Weiss temperature must be between 273 and 313 K")
    if not 0.0 <= salinity_psu <= 40.0:
        raise ValueError("salinity must be between 0 and 40")
    if seawater_density_kg_m3 <= 0.0:
        raise ValueError("seawater density must be positive")

    coefficients = species.coefficients
    scaled_temperature = temperature_k / 100.0
    log_k_mol_kg_atm = (
        coefficients.a1
        + coefficients.a2 * 100.0 / temperature_k
        + coefficients.a3 * log(scaled_temperature)
        + salinity_psu
        * (
            coefficients.b1
            + coefficients.b2 * scaled_temperature
            + coefficients.b3 * scaled_temperature**2
        )
    )
    return (
        exp(log_k_mol_kg_atm)
        * seawater_density_kg_m3
        / STANDARD_ATMOSPHERE_PA
    )


def chbr3_henry(temperature_k: float, config: UnifiedConfig) -> float:
    """Return CHBr3 Hcp in mol m-3 Pa-1."""
    return config.chbr3_henry_reference_mol_m3_pa * exp(
        config.chbr3_henry_temperature_parameter_k
        * (
            1.0 / temperature_k
            - 1.0 / config.chbr3_henry_reference_temperature_k
        )
    )


def water_ion_product_mol2_l2(temperature_k: float) -> float:
    """Return Kw(T) using the IAPWS 2024 equation at liquid density 1 g cm-3."""
    if not 273.15 <= temperature_k <= 1273.15:
        raise ValueError("water ion-product temperature must be 273.15-1273.15 K")

    density_g_cm3 = 1.0
    alpha_0 = -0.702132
    alpha_1_k = 8681.05
    alpha_2_k2 = -24_145.1
    beta_0_cm3_g = 0.813876
    beta_1_k_cm3_g = -51.4471
    beta_2_cm6_g2 = -0.469920
    coordination_number = 6.0
    water_molar_mass_g_mol = 18.015268

    z = density_g_cm3 * exp(
        alpha_0
        + alpha_1_k / temperature_k
        + (
            alpha_2_k2
            / temperature_k**2
            * density_g_cm3 ** (2.0 / 3.0)
        )
    )
    ideal_gas_pkw = (
        0.61415
        + 48_251.33 / temperature_k
        - 67_707.93 / temperature_k**2
        + 10_102_100.0 / temperature_k**3
    )
    density_term = (
        log(1.0 + z) / log(10.0)
        - z
        / (z + 1.0)
        * density_g_cm3
        * (
            beta_0_cm3_g
            + beta_1_k_cm3_g / temperature_k
            + beta_2_cm6_g2 * density_g_cm3
        )
    )
    pkw = (
        -2.0 * coordination_number * density_term
        + ideal_gas_pkw
        + 2.0 * log(water_molar_mass_g_mol / 1000.0) / log(10.0)
    )
    return 10.0**-pkw


def hydroxide_concentration_mol_l(ph: float, temperature_k: float) -> float:
    """Return liquid OH- from Kw(T) / H+."""
    if not 0.0 <= ph <= 14.0:
        raise ValueError("ph must remain within 0-14")
    hydrogen_mol_l = 10.0**-ph
    return water_ion_product_mol2_l2(temperature_k) / hydrogen_mol_l


def ph_scaled_hydrolysis_rate(
    ph: float,
    temperature_k: float,
    config: UnifiedConfig,
) -> float:
    """Return first-order aqueous CHBr3 loss from Arrhenius OH- hydrolysis."""
    hydroxide_mol_l = hydroxide_concentration_mol_l(ph, temperature_k)
    bimolecular_l_mol_s = (
        config.chbr3_hydrolysis_arrhenius_prefactor_l_mol_s
        * exp(
            -config.chbr3_hydrolysis_activation_energy_j_mol
            / (GAS_CONSTANT_J_MOL_K * temperature_k)
        )
    )
    return bimolecular_l_mol_s * hydroxide_mol_l


def chlorophyll_chbr3_emission_flux_molecule_cm2_s(
    chlorophyll_a_mg_m3: float,
    config: UnifiedConfig,
) -> float:
    """Return E = 1.127e5 * f * r * Chl-a for extra-tropical coasts."""
    if chlorophyll_a_mg_m3 < 0.0:
        raise ValueError("chlorophyll_a_mg_m3 cannot be negative")
    return (
        config.chlorophyll_emission_coefficient_molecule_cm2_s_per_mg_m3
        * config.chbr3_species_fraction
        * config.extra_tropical_coastal_modifier
        * chlorophyll_a_mg_m3
    )


def stratification_multiplier(
    gradient_richardson_number: float,
    config: UnifiedConfig,
) -> float:
    """Return power-law Gradient Richardson Number transfer suppression."""
    if gradient_richardson_number <= 0.0:
        return 1.0
    return (
        1.0
        + config.richardson_scaling_gamma * gradient_richardson_number
    ) ** -config.richardson_scaling_exponent


def wanninkhof_kw_m_s(
    wind_speed_m_s: float,
    gradient_richardson_number: float,
    config: UnifiedConfig,
) -> float:
    """Return Wanninkhof gas-transfer velocity with haline damping."""
    k_cm_hr = (
        config.wanninkhof_coefficient_cm_hr
        * wind_speed_m_s**2
        * (config.reference_schmidt_number / 660.0) ** -0.5
    )
    return (
        k_cm_hr
        * 0.01
        / 3600.0
        * stratification_multiplier(gradient_richardson_number, config)
    )


def _kinetic_rates(air_temperature_k: float, config: UnifiedConfig) -> dict[str, float]:
    return {
        "OH_CHBr3": 9.94e-13 * exp(-387.0 / air_temperature_k),
        "Cl_CHBr3": 0.43e-11 * exp(-809.0 / air_temperature_k),
        "BrO_HO2": 4.5e-12
        * exp(config.bro_ho2_activation_parameter_k / air_temperature_k),
    }


def _cfc_tendency(
    state: np.ndarray,
    species: TraceGas,
    target_ppt: float,
    temperature_k: float,
    salinity_psu: float,
    kw_m_s: float,
    config: UnifiedConfig,
) -> tuple[float, float, float]:
    air_mol = state[species.air_index]
    mixed_mol = state[species.mixed_index]
    deep_mol = state[species.deep_index]
    henry = warner_weiss_henry(
        species,
        temperature_k,
        salinity_psu,
        config.seawater_density_kg_m3,
    )
    partial_pressure_pa = (
        air_mol / config.tropospheric_dry_air_mol * config.surface_pressure_pa
    )
    equilibrium_mol_m3 = henry * partial_pressure_pa
    mixed_mol_m3 = mixed_mol / config.mixed_layer_volume_m3
    deep_mol_m3 = deep_mol / config.deep_ocean_volume_m3

    air_to_sea_mol_s = (
        kw_m_s
        * config.ocean_area_m2
        * (equilibrium_mol_m3 - mixed_mol_m3)
    )
    mixed_to_deep_mol_s = (
        config.mixed_layer_volume_m3
        / config.vertical_exchange_timescale_s
        * (mixed_mol_m3 - deep_mol_m3)
    )
    target_air_mol = target_ppt * 1.0e-12 * config.tropospheric_dry_air_mol
    external_mol_s = (
        target_air_mol - air_mol
    ) / config.atmospheric_nudging_timescale_s
    return air_to_sea_mol_s, mixed_to_deep_mol_s, external_mol_s


def initial_state(
    forcing: CoupledForcing,
    config: UnifiedConfig | None = None,
) -> np.ndarray:
    """Construct the common initial state for all matched integrations."""
    config = config or UnifiedConfig()
    state = np.zeros(STATE_SIZE, dtype=float)
    state[:7] = (
        config.initial_chbr3_gas_molecule_cm3,
        config.initial_br_molecule_cm3,
        config.initial_bro_molecule_cm3,
        config.initial_ho2_molecule_cm3,
        config.initial_oh_molecule_cm3,
        config.initial_hobr_molecule_cm3,
        config.initial_chbr3_aqueous_mol_m3,
    )

    temperature_k = float(forcing.sst_k[0])
    salinity_psu = float(forcing.salinity_psu[0])
    targets = (float(forcing.cfc11_ppt[0]), float(forcing.cfc12_ppt[0]))
    for species, target_ppt in zip(CFC_SPECIES, targets, strict=True):
        air_mol = target_ppt * 1.0e-12 * config.tropospheric_dry_air_mol
        partial_pressure_pa = (
            air_mol / config.tropospheric_dry_air_mol * config.surface_pressure_pa
        )
        equilibrium_mol_m3 = warner_weiss_henry(
            species,
            temperature_k,
            salinity_psu,
            config.seawater_density_kg_m3,
        ) * partial_pressure_pa
        state[species.air_index] = air_mol
        state[species.mixed_index] = (
            equilibrium_mol_m3 * config.mixed_layer_volume_m3
        )
        state[species.deep_index] = (
            equilibrium_mol_m3
            * config.initial_deep_saturation_fraction
            * config.deep_ocean_volume_m3
        )
    return state


def unified_derivatives(
    time_s: float,
    state: np.ndarray,
    forcing: CoupledForcing,
    config: UnifiedConfig | None = None,
    switches: ModelSwitches | None = None,
) -> np.ndarray:
    """Return all 17 state derivatives at one solver time."""
    config = config or UnifiedConfig()
    switches = switches or ModelSwitches()

    sst_k = forcing.value("sst_k", time_s)
    air_temperature_k = forcing.value("air_temperature_k", time_s)
    wind_speed_m_s = forcing.value("wind_speed_m_s", time_s)
    salinity_psu = forcing.value("salinity_psu", time_s)
    gradient_richardson_number = (
        forcing.value("gradient_richardson_number", time_s)
        if switches.haline_stratification
        else 0.0
    )
    ph = (
        forcing.value("ph", time_s)
        if switches.ph_hydrolysis
        else config.chbr3_hydrolysis_reference_ph
    )
    chlorophyll_a = (
        forcing.value("chlorophyll_a_mg_m3", time_s)
        if switches.chlorophyll_source
        else 0.0
    )
    kw_m_s = wanninkhof_kw_m_s(
        wind_speed_m_s,
        gradient_richardson_number,
        config,
    )

    chbr3, br, bro, ho2, oh, hobr, chbr3_aqueous = state[:7]
    rates = _kinetic_rates(air_temperature_k, config)
    hydrolysis_s = ph_scaled_hydrolysis_rate(ph, sst_k, config)
    biological_emission_molecule_cm2_s = (
        chlorophyll_chbr3_emission_flux_molecule_cm2_s(
            chlorophyll_a,
            config,
        )
    )
    biological_source_mol_m3_s = (
        biological_emission_molecule_cm2_s
        * 1.0e4
        / AVOGADRO_MOL
        / config.mixed_layer_depth_m
    )

    chbr3_air_mol_m3 = chbr3 * 1.0e6 / AVOGADRO_MOL
    chbr3_partial_pressure_pa = chbr3_air_mol_m3 * GAS_CONSTANT_J_MOL_K * sst_k
    chbr3_equilibrium_mol_m3 = (
        chbr3_henry(sst_k, config) * chbr3_partial_pressure_pa
    )
    chbr3_outgassing_mol_m2_s = kw_m_s * (
        chbr3_aqueous - chbr3_equilibrium_mol_m3
    )
    chbr3_ocean_source_molecule_cm3_s = (
        chbr3_outgassing_mol_m2_s
        / config.marine_boundary_layer_height_m
        * AVOGADRO_MOL
        / 1.0e6
    )

    chbr3_loss = (
        config.chbr3_photolysis_s
        + rates["OH_CHBr3"] * oh
        + rates["Cl_CHBr3"] * config.chlorine_molecule_cm3
    ) * chbr3
    br_to_bro = (
        config.br_o3_rate_cm3_molecule_s
        * config.ozone_molecule_cm3
        * br
    )
    hobr_formation = rates["BrO_HO2"] * bro * ho2
    hobr_photolysis = config.hobr_photolysis_s * hobr
    hobr_termination = config.hobr_termination_s * hobr
    br_production = 3.0 * config.chbr3_photolysis_s * chbr3 + hobr_photolysis

    derivative = np.zeros(STATE_SIZE, dtype=float)
    derivative[:7] = (
        chbr3_ocean_source_molecule_cm3_s - chbr3_loss,
        br_production - br_to_bro,
        br_to_bro - hobr_formation,
        (
            config.hox_production_molecule_cm3_s
            - hobr_formation
            - config.ho2_background_loss_s * ho2
        ),
        (
            config.hox_production_molecule_cm3_s
            + hobr_photolysis
            - config.oh_background_loss_s * oh
        ),
        hobr_formation - hobr_photolysis - hobr_termination,
        (
            biological_source_mol_m3_s
            - hydrolysis_s * chbr3_aqueous
            - chbr3_outgassing_mol_m2_s / config.mixed_layer_depth_m
        ),
    )

    target_ppt = (
        forcing.value("cfc11_ppt", time_s),
        forcing.value("cfc12_ppt", time_s),
    )
    for species, target in zip(CFC_SPECIES, target_ppt, strict=True):
        air_to_sea, mixed_to_deep, external = _cfc_tendency(
            state,
            species,
            target,
            sst_k,
            salinity_psu,
            kw_m_s,
            config,
        )
        derivative[species.air_index] = -air_to_sea + external
        derivative[species.mixed_index] = air_to_sea - mixed_to_deep
        derivative[species.deep_index] = mixed_to_deep
        derivative[species.cumulative_flux_index] = air_to_sea
        derivative[species.cumulative_external_index] = external
    return derivative


def _absolute_tolerances() -> np.ndarray:
    tolerances = np.full(STATE_SIZE, 1.0e-3)
    tolerances[StateIndex.CHBR3_AQUEOUS] = 1.0e-16
    return tolerances


def _integrate(
    forcing: CoupledForcing,
    initial: np.ndarray,
    config: UnifiedConfig,
    switches: ModelSwitches,
) -> Any:
    solution = solve_ivp(
        unified_derivatives,
        (float(forcing.time_s[0]), float(forcing.time_s[-1])),
        initial,
        args=(forcing, config, switches),
        method="Radau",
        t_eval=forcing.time_s,
        rtol=config.relative_tolerance,
        atol=_absolute_tolerances(),
        max_step=config.maximum_step_s,
    )
    if not solution.success:
        raise RuntimeError(f"Radau integration failed: {solution.message}")
    if not np.isfinite(solution.y).all():
        raise RuntimeError("Radau integration produced a non-finite state")

    chemical_minima = np.min(solution.y[:6], axis=1)
    if np.any(chemical_minima < -1.0e-3):
        raise RuntimeError("Radau integration produced negative gas concentrations")
    if np.min(solution.y[StateIndex.CHBR3_AQUEOUS]) < -1.0e-15:
        raise RuntimeError("Radau integration produced negative aqueous CHBr3")
    inventory_indices = (
        StateIndex.CFC11_AIR,
        StateIndex.CFC11_MIXED,
        StateIndex.CFC11_DEEP,
        StateIndex.CFC12_AIR,
        StateIndex.CFC12_MIXED,
        StateIndex.CFC12_DEEP,
    )
    if np.min(solution.y[list(inventory_indices)]) < -1.0e-3:
        raise RuntimeError("Radau integration produced negative CFC inventories")
    return solution


def integrate_unified_state_vector(
    forcing: CoupledForcing,
    initial: np.ndarray | None = None,
    config: UnifiedConfig | None = None,
) -> UnifiedResult:
    """Integrate the coupled model and return both requested effect matrices.

    The full trajectory is paired with two controlled counterfactuals:

    * ``no_stratification_state`` keeps pH and chlorophyll forcing but sets the
      Gradient Richardson Number to zero in ``kw``. Its cumulative CFC uptake
      difference from the full run is the inversion discrepancy.
    * ``no_ph_biology_state`` retains physical stratification but fixes pH at
      the hydrolysis reference and removes the chlorophyll-scaled CHBr3 source.
      The OH ratio gives the methane-lifetime distortion, assuming the methane
      loss rate is proportional to OH.
    """
    config = config or UnifiedConfig()
    y0 = initial_state(forcing, config) if initial is None else np.asarray(
        initial,
        dtype=float,
    ).copy()
    if y0.shape != (STATE_SIZE,):
        raise ValueError(f"initial must have shape ({STATE_SIZE},)")
    if not np.isfinite(y0).all():
        raise ValueError("initial must contain only finite values")

    actual = _integrate(forcing, y0, config, ModelSwitches())
    no_stratification = _integrate(
        forcing,
        y0,
        config,
        ModelSwitches(haline_stratification=False),
    )
    no_ph_biology = _integrate(
        forcing,
        y0,
        config,
        ModelSwitches(ph_hydrolysis=False, chlorophyll_source=False),
    )

    cfc_discrepancies: list[np.ndarray] = []
    for species in CFC_SPECIES:
        discrepancy_mol = (
            no_stratification.y[species.cumulative_flux_index]
            - actual.y[species.cumulative_flux_index]
        )
        cfc_discrepancies.append(
            discrepancy_mol * species.molecular_weight_g_mol / 1.0e9
        )
    gradient_richardson_number = forcing.gradient_richardson_number
    kw_multiplier = np.array(
        [
            stratification_multiplier(value, config)
            for value in gradient_richardson_number
        ]
    )
    cfc_matrix = np.column_stack(
        (
            forcing.time_s,
            gradient_richardson_number,
            kw_multiplier,
            cfc_discrepancies[0],
            cfc_discrepancies[1],
        )
    )

    dynamic_oh = actual.y[StateIndex.OH]
    baseline_oh = no_ph_biology.y[StateIndex.OH]
    lifetime_distortion_percent = (
        baseline_oh / np.maximum(dynamic_oh, np.finfo(float).tiny) - 1.0
    ) * 100.0
    methane_matrix = np.column_stack(
        (
            forcing.time_s,
            forcing.ph,
            forcing.chlorophyll_a_mg_m3,
            dynamic_oh,
            baseline_oh,
            lifetime_distortion_percent,
        )
    )

    return UnifiedResult(
        time_s=forcing.time_s.copy(),
        state=actual.y.copy(),
        no_stratification_state=no_stratification.y.copy(),
        no_ph_biology_state=no_ph_biology.y.copy(),
        cfc_inversion_discrepancy_matrix=cfc_matrix,
        methane_lifetime_distortion_matrix=methane_matrix,
    )


__all__ = [
    "CFC11",
    "CFC12",
    "CoupledForcing",
    "ModelSwitches",
    "STATE_NAMES",
    "STATE_SIZE",
    "StateIndex",
    "UnifiedConfig",
    "UnifiedResult",
    "chlorophyll_chbr3_emission_flux_molecule_cm2_s",
    "chbr3_henry",
    "hydroxide_concentration_mol_l",
    "initial_state",
    "integrate_unified_state_vector",
    "ph_scaled_hydrolysis_rate",
    "stratification_multiplier",
    "unified_derivatives",
    "wanninkhof_kw_m_s",
    "warner_weiss_henry",
    "water_ion_product_mol2_l2",
]
