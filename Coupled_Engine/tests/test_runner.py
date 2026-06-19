import csv

import numpy as np

from engine import STATE_SIZE, UnifiedConfig, UnifiedResult
from engine.runner import (
    CLIMATE_ACCOUNTING_ERROR_COLUMNS,
    DAY_S,
    Experiment,
    build_forcing,
    climate_accounting_error_matrix,
    write_outputs,
)


def make_result() -> UnifiedResult:
    time_s = np.array([0.0, DAY_S, 2.0 * DAY_S])
    cfc_matrix = np.column_stack(
        (
            time_s,
            np.array([1.0, 1.0, 1.0]),
            np.array([0.78, 0.78, 0.78]),
            np.array([0.0, 0.1, 0.2]),
            np.array([0.0, 0.03, 0.06]),
        )
    )
    methane_matrix = np.column_stack(
        (
            time_s,
            np.array([7.7, 7.7, 7.7]),
            np.array([2.0, 2.0, 2.0]),
            np.array([3.0e5, 3.1e5, 3.2e5]),
            np.array([2.9e5, 2.9e5, 2.9e5]),
            np.array([-3.0, -4.0, -5.0]),
        )
    )
    state = np.zeros((STATE_SIZE, time_s.size))
    return UnifiedResult(
        time_s=time_s,
        state=state,
        no_stratification_state=state.copy(),
        no_ph_biology_state=state.copy(),
        cfc_inversion_discrepancy_matrix=cfc_matrix,
        methane_lifetime_distortion_matrix=methane_matrix,
    )


def test_climate_accounting_error_matrix_joins_on_time_step() -> None:
    matrix = climate_accounting_error_matrix(make_result())

    assert matrix.shape == (3, len(CLIMATE_ACCOUNTING_ERROR_COLUMNS))
    np.testing.assert_allclose(matrix[:, 0], [0.0, DAY_S, 2.0 * DAY_S])
    np.testing.assert_allclose(matrix[:, 2], [0.0, DAY_S, DAY_S])
    np.testing.assert_allclose(matrix[:, 6], [0.0, 0.1, 0.2])
    np.testing.assert_allclose(matrix[:, 10], [-3.0, -4.0, -5.0])


def test_build_forcing_reports_default_fallback(tmp_path) -> None:
    experiment = Experiment(duration_days=2.0)
    forcing, provenance = build_forcing(experiment, tmp_path)

    assert provenance.mode == "deterministic_validation_defaults"
    assert provenance.missing_inputs == ("ph", "chlorophyll", "salinity")
    np.testing.assert_allclose(forcing.ph, experiment.ph)
    np.testing.assert_allclose(
        forcing.chlorophyll_a_mg_m3,
        experiment.chlorophyll_a_mg_m3,
    )


def test_write_outputs_creates_paper_artifacts(tmp_path) -> None:
    output_dir = tmp_path / "results"
    experiment = Experiment(duration_days=2.0)
    _, provenance = build_forcing(experiment, tmp_path)

    write_outputs(
        output_dir,
        make_result(),
        UnifiedConfig(),
        experiment,
        provenance,
    )

    expected_files = {
        "cfc_inversion_discrepancy_matrix.csv",
        "methane_lifetime_distortion_matrix.csv",
        "climate_accounting_error_matrix.csv",
        "unified_state_trajectory.csv",
        "parameters.csv",
        "summary.md",
    }
    assert expected_files == {path.name for path in output_dir.iterdir()}

    with (output_dir / "climate_accounting_error_matrix.csv").open(
        newline="",
        encoding="utf-8",
    ) as handle:
        reader = csv.reader(handle)
        assert next(reader) == list(CLIMATE_ACCOUNTING_ERROR_COLUMNS)

    summary = (output_dir / "summary.md").read_text(encoding="utf-8")
    assert "Reported-timestep sum, CFC-11 discrepancy column: 0.3 Gg" in summary
    assert "Forcing mode: `deterministic_validation_defaults`" in summary
