# ============================================================================
# Copyright (c) 2026 Brandon W. Caris. All rights reserved.
# Part of the Coupled Marine Engine Suite (in_your_dreams/Coupled_Engine)
# Licensed under CC BY 4.0 with Explicit Academic Citation Mandate.
# See master LICENSE file at repository root for full terms and string.
# ============================================================================
"""Stiff marine boundary-layer box model for CHBr3 and reactive bromine."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class MarineState:
    """Environmental forcing for one model scenario."""

    name: str
    sst_k: float
    air_temperature_k: float
    wind_speed_m_s: float


class HalogenBoxModel:
    """Six-species conceptual CHBr3-Br-BrO-HO2-OH-HOBr box model."""

    state_names = ("CHBr3", "Br", "BrO", "HO2", "OH", "HOBr")

    # Solubility-form Henry coefficient, Hcp = Caq / p.
    kh_reference = 1.9e-2  # mol m-3 Pa-1 at 298.15 K
    kh_reference_temperature_k = 298.15
    kh_temperature_parameter_k = 4700.0

    dissolved_chbr3_mol_m3 = 5.0e-9
    chbr3_photolysis_s = 2.0e-6
    hobr_photolysis_s = 5.0e-3
    hox_production_molecule_cm3_s = 1.0e6
    chlorine_molecule_cm3 = 1.0e4
    ozone_molecule_cm3 = 7.0e11
    br_o3_rate_cm3_molecule_s = 1.2e-12

    # Retained prototype scaling assumptions.
    gas_transfer_prefactor = 0.01
    ocean_source_scale = 1.0e10
    ho2_background_loss_s = 5.0e-2
    oh_background_loss_s = 4.0
    bro_ho2_activation_parameter_k = 580.0
    hobr_termination_s = 1.0e-4
    steady_state_fractional_drift_per_hour = 1.0e-4

    def __init__(self) -> None:
        self.initial_state = np.array(
            [1.0e9, 0.0, 0.0, 1.0e8, 1.0e6, 0.0]
        )

    def henry_chbr3(self, temperature_k: float) -> float:
        """Return CHBr3 Hcp in mol m-3 Pa-1."""
        return self.kh_reference * np.exp(
            self.kh_temperature_parameter_k
            * (
                1.0 / temperature_k
                - 1.0 / self.kh_reference_temperature_k
            )
        )

    def kinetic_rates(self, air_temperature_k: float) -> dict[str, float]:
        """Return bimolecular rates in cm3 molecule-1 s-1."""
        return {
            "OH_CHBr3": 9.94e-13 * np.exp(-387.0 / air_temperature_k),
            "Cl_CHBr3": 0.43e-11 * np.exp(-809.0 / air_temperature_k),
            "BrO_HO2": 4.5e-12
            * np.exp(
                self.bro_ho2_activation_parameter_k / air_temperature_k
            ),
        }

    def ocean_source_rate(self, scenario: MarineState) -> float:
        """Return the prototype ocean source in molecule cm-3 s-1."""
        equilibrium_pressure_pa = (
            self.dissolved_chbr3_mol_m3 / self.henry_chbr3(scenario.sst_k)
        )
        transfer_proxy = (
            self.gas_transfer_prefactor * scenario.wind_speed_m_s**2
        )
        return (
            transfer_proxy
            * equilibrium_pressure_pa
            * self.ocean_source_scale
        )

    def tendency(
        self,
        _time_s: float,
        state: np.ndarray,
        scenario: MarineState,
    ) -> np.ndarray:
        """Return tendencies for CHBr3, Br, BrO, HO2, OH, and HOBr."""
        chbr3, br, bro, ho2, oh, hobr = state
        rates = self.kinetic_rates(scenario.air_temperature_k)
        k_term_HOBr = self.hobr_termination_s  # Approx. 2.78-hour lifetime

        chbr3_loss = (
            self.chbr3_photolysis_s
            + rates["OH_CHBr3"] * oh
            + rates["Cl_CHBr3"] * self.chlorine_molecule_cm3
        ) * chbr3

        br_to_bro = (
            self.br_o3_rate_cm3_molecule_s
            * self.ozone_molecule_cm3
            * br
        )
        hobr_formation = rates["BrO_HO2"] * bro * ho2
        hobr_photolysis = self.hobr_photolysis_s * hobr
        hobr_termination = k_term_HOBr * hobr
        br_production = (
            3.0 * self.chbr3_photolysis_s * chbr3 + hobr_photolysis
        )

        return np.array(
            [
                self.ocean_source_rate(scenario) - chbr3_loss,
                br_production - br_to_bro,
                br_to_bro - hobr_formation,
                (
                    self.hox_production_molecule_cm3_s
                    - hobr_formation
                    - self.ho2_background_loss_s * ho2
                ),
                (
                    self.hox_production_molecule_cm3_s
                    + hobr_photolysis
                    - self.oh_background_loss_s * oh
                ),
                hobr_formation - hobr_photolysis - hobr_termination,
            ]
        )

    def fractional_drift_per_hour(
        self,
        time_s: float,
        state: np.ndarray,
        scenario: MarineState,
    ) -> np.ndarray:
        """Return absolute fractional state drift over one hour."""
        concentration_scale = np.maximum(np.abs(state), 1.0)
        return (
            np.abs(self.tendency(time_s, state, scenario))
            / concentration_scale
            * 3600.0
        )

    def is_steady_state(
        self,
        time_s: float,
        state: np.ndarray,
        scenario: MarineState,
    ) -> bool:
        """Return whether every state meets the fractional-drift criterion."""
        return bool(
            np.all(
                self.fractional_drift_per_hour(time_s, state, scenario)
                <= self.steady_state_fractional_drift_per_hour
            )
        )

    def run(
        self,
        scenario: MarineState,
        duration_s: float = 86_400.0,
        output_interval_s: float = 300.0,
    ):
        """Integrate one scenario with SciPy's implicit Radau solver."""
        output_times = np.arange(
            0.0,
            duration_s + output_interval_s,
            output_interval_s,
        )
        output_times = output_times[output_times <= duration_s]

        solution = solve_ivp(
            self.tendency,
            (0.0, duration_s),
            self.initial_state,
            args=(scenario,),
            method="Radau",
            t_eval=output_times,
            rtol=1.0e-8,
            atol=1.0e-3,
        )
        if not solution.success:
            raise RuntimeError(solution.message)
        if not np.isfinite(solution.y).all() or np.min(solution.y) < -1.0e-3:
            raise RuntimeError("Solver produced an invalid concentration")
        return solution
