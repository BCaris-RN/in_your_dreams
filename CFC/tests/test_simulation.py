# ============================================================================
# Copyright (c) 2026 Brandon W. Caris. All rights reserved.
# Part of the Coupled Marine Engine Suite (in_your_dreams/Coupled_Engine)
# Licensed under CC BY 4.0 with Explicit Academic Citation Mandate.
# See master LICENSE file at repository root for full terms and string.
# ============================================================================
import numpy as np
import pandas as pd

from cfc_flux_model.forcing import Forcing
from cfc_flux_model.parameters import CFC11, ModelConfig
from cfc_flux_model.simulation import initial_state, run_simulation


def make_forcing(anomaly: float, months: int = 13) -> Forcing:
    index = pd.period_range("2020-01", periods=months, freq="M")
    frame = pd.DataFrame(
        {
            "sst_anomaly_c": np.full(months, anomaly),
            "cfc11_ppt": np.linspace(230.0, 228.0, months),
            "cfc12_ppt": np.linspace(510.0, 507.0, months),
            "sst_estimated": False,
            "cfc11_estimated": False,
            "cfc12_estimated": False,
        },
        index=index,
    )
    frame.index.name = "month"
    return Forcing(frame)


def test_coupled_inventory_closes_against_external_forcing() -> None:
    result = run_simulation(
        make_forcing(0.5),
        CFC11,
        ModelConfig(),
        scenario="test",
    )

    assert result.frame["mass_balance_residual_mol"].abs().max() < 0.01


def test_zero_anomaly_scenarios_are_numerically_identical() -> None:
    forcing = make_forcing(0.0)
    state = initial_state(forcing, CFC11, ModelConfig())
    observed = run_simulation(
        forcing,
        CFC11,
        ModelConfig(),
        scenario="observed",
        initial=state,
        use_sst_anomaly=True,
    )
    counterfactual = run_simulation(
        forcing,
        CFC11,
        ModelConfig(),
        scenario="counterfactual",
        initial=state,
        use_sst_anomaly=False,
    )

    numeric_columns = observed.frame.select_dtypes(include="number").columns
    np.testing.assert_allclose(
        observed.frame[numeric_columns],
        counterfactual.frame[numeric_columns],
        rtol=0.0,
        atol=0.0,
    )


def test_positive_sst_anomaly_shifts_flux_toward_atmosphere() -> None:
    zero_forcing = make_forcing(0.0)
    warm_forcing = make_forcing(1.0)
    state = initial_state(zero_forcing, CFC11, ModelConfig())
    warm = run_simulation(
        warm_forcing,
        CFC11,
        ModelConfig(),
        scenario="warm",
        initial=state,
        use_sst_anomaly=True,
    )
    zero = run_simulation(
        zero_forcing,
        CFC11,
        ModelConfig(),
        scenario="zero",
        initial=state,
        use_sst_anomaly=False,
    )

    warm_net_release = -warm.frame["cumulative_air_to_sea_mol"].iloc[-1]
    zero_net_release = -zero.frame["cumulative_air_to_sea_mol"].iloc[-1]
    assert warm_net_release > zero_net_release
