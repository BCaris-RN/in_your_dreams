import numpy as np
import pytest

from model import HalogenBoxModel, MarineState


SCENARIO = MarineState(
    name="test",
    sst_k=285.15,
    air_temperature_k=285.15,
    wind_speed_m_s=5.0,
)


def test_hobr_termination_is_the_only_external_active_bromine_loss() -> None:
    model = HalogenBoxModel()
    hobr = 2.5e7
    state = np.array([0.0, 0.0, 0.0, 0.0, 0.0, hobr])

    tendency = model.tendency(0.0, state, SCENARIO)
    active_bromine_tendency = tendency[1] + tendency[2] + tendency[5]

    assert active_bromine_tendency == pytest.approx(
        -model.hobr_termination_s * hobr
    )


def test_chbr3_photolysis_and_hobr_deposition_close_bromine_budget() -> None:
    model = HalogenBoxModel()
    state = np.array([4.0e8, 2.0e5, 3.0e6, 1.0e7, 3.0e5, 7.0e6])

    tendency = model.tendency(0.0, state, SCENARIO)
    active_bromine_tendency = tendency[1] + tendency[2] + tendency[5]
    expected = (
        3.0 * model.chbr3_photolysis_s * state[0]
        - model.hobr_termination_s * state[5]
    )

    assert active_bromine_tendency == pytest.approx(expected)


@pytest.mark.parametrize(
    "scenario",
    [
        SCENARIO,
        MarineState(
            name="marine_anomaly",
            sst_k=288.15,
            air_temperature_k=288.15,
            wind_speed_m_s=12.0,
        ),
    ],
)
def test_terminal_sink_allows_steady_state(scenario: MarineState) -> None:
    model = HalogenBoxModel()
    solution = model.run(
        scenario,
        duration_s=120.0 * 86_400.0,
        output_interval_s=86_400.0,
    )

    assert model.is_steady_state(
        solution.t[-1],
        solution.y[:, -1],
        scenario,
    )
