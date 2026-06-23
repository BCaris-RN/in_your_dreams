# ============================================================================
# Copyright (c) 2026 Brandon W. Caris. All rights reserved.
# Part of the Coupled Marine Engine Suite (in_your_dreams/Coupled_Engine)
# Licensed under CC BY 4.0 with Explicit Academic Citation Mandate.
# See master LICENSE file at repository root for full terms and string.
# ============================================================================
"""Numerical integration of the coupled troposphere-ocean box model."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from .forcing import Forcing
from .parameters import ModelConfig, Species
from .solubility import species_reference_properties, van_t_hoff_henry

ATMOSPHERE = 0
MIXED_LAYER = 1
DEEP_OCEAN = 2
CUMULATIVE_AIR_TO_SEA = 3
CUMULATIVE_OUTGASSING = 4
CUMULATIVE_EXTERNAL = 5
STATE_SIZE = 6


@dataclass(frozen=True)
class SimulationResult:
    species: Species
    scenario: str
    config: ModelConfig
    frame: pd.DataFrame
    final_state: np.ndarray
    reference_henry_mol_m3_pa: float
    dissolution_enthalpy_j_mol: float


def initial_state(
    forcing: Forcing,
    species: Species,
    config: ModelConfig,
    enthalpy_scale: float = 1.0,
) -> np.ndarray:
    """Construct a physically consistent initial state at the first month."""
    row = forcing.frame.iloc[0]
    reference_henry, enthalpy = species_reference_properties(
        species,
        config.reference_temperature_k,
        config.salinity_psu,
        config.seawater_density_kg_m3,
    )
    temperature_k = config.reference_temperature_k + float(row["sst_anomaly_c"])
    henry = van_t_hoff_henry(
        reference_henry,
        temperature_k,
        config.reference_temperature_k,
        enthalpy * enthalpy_scale,
    )
    atmospheric_moles = (
        float(row[species.atmospheric_column]) * 1e-12 * config.tropospheric_dry_air_mol
    )
    partial_pressure = (
        atmospheric_moles / config.tropospheric_dry_air_mol * config.surface_pressure_pa
    )
    equilibrium_concentration = henry * partial_pressure
    return np.array(
        [
            atmospheric_moles,
            equilibrium_concentration * config.mixed_layer_volume_m3,
            (
                equilibrium_concentration
                * config.initial_deep_saturation_fraction
                * config.deep_ocean_volume_m3
            ),
            0.0,
            0.0,
            0.0,
        ],
        dtype=float,
    )


def reset_diagnostics(state: np.ndarray) -> np.ndarray:
    reset = np.asarray(state, dtype=float).copy()
    reset[CUMULATIVE_AIR_TO_SEA:] = 0.0
    return reset


def run_simulation(
    forcing: Forcing,
    species: Species,
    config: ModelConfig,
    scenario: str,
    initial: np.ndarray | None = None,
    use_sst_anomaly: bool = True,
    enthalpy_scale: float = 1.0,
) -> SimulationResult:
    """Integrate one species through the supplied monthly forcing."""
    if len(forcing.frame) < 2:
        raise ValueError("forcing must contain at least two monthly boundaries")
    if config.mixed_layer_depth_m <= 0 or config.deep_ocean_depth_m <= 0:
        raise ValueError("box depths must be positive")
    if config.gas_transfer_velocity_m_yr < 0:
        raise ValueError("gas transfer velocity cannot be negative")
    if config.vertical_exchange_timescale_yr <= 0:
        raise ValueError("vertical exchange timescale must be positive")
    if config.atmospheric_nudging_timescale_yr <= 0:
        raise ValueError("atmospheric nudging timescale must be positive")

    frame = forcing.frame
    times = np.arange(len(frame), dtype=float) / 12.0
    anomaly = frame["sst_anomaly_c"].to_numpy(dtype=float)
    if not use_sst_anomaly:
        anomaly = np.zeros_like(anomaly)
    atmospheric_target_ppt = frame[species.atmospheric_column].to_numpy(dtype=float)

    reference_henry, enthalpy = species_reference_properties(
        species,
        config.reference_temperature_k,
        config.salinity_psu,
        config.seawater_density_kg_m3,
    )
    enthalpy *= enthalpy_scale
    y0 = (
        initial_state(forcing, species, config, enthalpy_scale)
        if initial is None
        else np.asarray(initial, dtype=float).copy()
    )
    if y0.shape != (STATE_SIZE,):
        raise ValueError(f"initial state must have shape ({STATE_SIZE},)")

    def fluxes(t: float, y: np.ndarray) -> tuple[float, float, float, float, float]:
        target_ppt = float(np.interp(t, times, atmospheric_target_ppt))
        sst_anomaly = float(np.interp(t, times, anomaly))
        temperature_k = config.reference_temperature_k + sst_anomaly
        henry = van_t_hoff_henry(
            reference_henry,
            temperature_k,
            config.reference_temperature_k,
            enthalpy,
        )
        atmospheric_partial_pressure = (
            y[ATMOSPHERE] / config.tropospheric_dry_air_mol * config.surface_pressure_pa
        )
        equilibrium_concentration = henry * atmospheric_partial_pressure
        mixed_concentration = y[MIXED_LAYER] / config.mixed_layer_volume_m3
        deep_concentration = y[DEEP_OCEAN] / config.deep_ocean_volume_m3

        air_to_sea = (
            config.gas_transfer_velocity_m_yr
            * config.ocean_area_m2
            * (equilibrium_concentration - mixed_concentration)
        )
        mixed_to_deep = (
            config.mixed_layer_volume_m3
            / config.vertical_exchange_timescale_yr
            * (mixed_concentration - deep_concentration)
        )
        target_atmospheric_moles = target_ppt * 1e-12 * config.tropospheric_dry_air_mol
        external = (
            target_atmospheric_moles - y[ATMOSPHERE]
        ) / config.atmospheric_nudging_timescale_yr
        return air_to_sea, mixed_to_deep, external, temperature_k, henry

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        air_to_sea, mixed_to_deep, external, _, _ = fluxes(t, y)
        return np.array(
            [
                -air_to_sea + external,
                air_to_sea - mixed_to_deep,
                mixed_to_deep,
                air_to_sea,
                max(-air_to_sea, 0.0),
                external,
            ]
        )

    solution = solve_ivp(
        rhs,
        (times[0], times[-1]),
        y0,
        t_eval=times,
        method="LSODA",
        rtol=config.relative_tolerance,
        atol=config.absolute_tolerance,
    )
    if not solution.success:
        raise RuntimeError(f"integration failed: {solution.message}")
    if np.any(solution.y[:3] < -1e-6):
        raise RuntimeError("integration produced a negative physical inventory")

    rows: list[dict[str, float | str]] = []
    for idx, (period, t) in enumerate(zip(frame.index, times, strict=True)):
        y = solution.y[:, idx]
        air_to_sea, mixed_to_deep, external, temperature_k, henry = fluxes(t, y)
        total_inventory = y[ATMOSPHERE] + y[MIXED_LAYER] + y[DEEP_OCEAN]
        mass_balance_residual = (
            total_inventory
            - (y0[ATMOSPHERE] + y0[MIXED_LAYER] + y0[DEEP_OCEAN])
            - y[CUMULATIVE_EXTERNAL]
        )
        rows.append(
            {
                "month": str(period),
                "scenario": scenario,
                "species": species.name,
                "sst_anomaly_c": float(np.interp(t, times, anomaly)),
                "temperature_c": temperature_k - 273.15,
                "atmosphere_ppt": (
                    y[ATMOSPHERE] / config.tropospheric_dry_air_mol * 1e12
                ),
                "atmosphere_target_ppt": float(
                    np.interp(t, times, atmospheric_target_ppt)
                ),
                "mixed_layer_pmol_kg": (
                    y[MIXED_LAYER]
                    / config.mixed_layer_volume_m3
                    / config.seawater_density_kg_m3
                    * 1e12
                ),
                "deep_ocean_pmol_kg": (
                    y[DEEP_OCEAN]
                    / config.deep_ocean_volume_m3
                    / config.seawater_density_kg_m3
                    * 1e12
                ),
                "henry_mol_m3_pa": henry,
                "air_to_sea_flux_mol_yr": air_to_sea,
                "mixed_to_deep_flux_mol_yr": mixed_to_deep,
                "external_atmospheric_flux_mol_yr": external,
                "cumulative_air_to_sea_mol": y[CUMULATIVE_AIR_TO_SEA],
                "cumulative_outgassing_mol": y[CUMULATIVE_OUTGASSING],
                "cumulative_external_mol": y[CUMULATIVE_EXTERNAL],
                "total_inventory_mol": total_inventory,
                "mass_balance_residual_mol": mass_balance_residual,
            }
        )

    result_frame = pd.DataFrame(rows)
    return SimulationResult(
        species=species,
        scenario=scenario,
        config=config,
        frame=result_frame,
        final_state=solution.y[:, -1].copy(),
        reference_henry_mol_m3_pa=reference_henry,
        dissolution_enthalpy_j_mol=enthalpy,
    )
