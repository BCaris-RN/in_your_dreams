# ============================================================================
# Copyright (c) 2026 Brandon W. Caris. All rights reserved.
# Part of the Coupled Marine Engine Suite (in_your_dreams/Coupled_Engine)
# Licensed under CC BY 4.0 with Explicit Academic Citation Mandate.
# See master LICENSE file at repository root for full terms and string.
# ============================================================================
"""Command-line runner and artifact writer for the coupled engine."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Sequence

import numpy as np

from . import (
    CFC11,
    CFC12,
    STATE_NAMES,
    CoupledForcing,
    UnifiedConfig,
    UnifiedResult,
    integrate_unified_state_vector,
)

DAY_S = 86_400.0

CLIMATE_ACCOUNTING_ERROR_COLUMNS = (
    "time_s",
    "day",
    "delta_t_s",
    "delta_t_days",
    "gradient_richardson_number",
    "kw_stratification_multiplier",
    "cfc11_inversion_discrepancy_Gg",
    "cfc12_inversion_discrepancy_Gg",
    "oh_dynamic_molecule_cm3",
    "oh_no_ph_biology_molecule_cm3",
    "methane_lifetime_distortion_percent",
)


@dataclass(frozen=True)
class Experiment:
    """Deterministic Paper 3 validation scenario."""

    name: str = "deterministic_validation"
    duration_days: float = 60.0
    output_interval_days: float = 1.0
    sst_k: float = 288.15
    air_temperature_k: float = 288.15
    wind_speed_m_s: float = 8.0
    salinity_psu: float = 35.0
    gradient_richardson_number: float = 1.0
    ph: float = 7.7
    chlorophyll_a_mg_m3: float = 2.0
    cfc11_start_ppt: float = 230.0
    cfc11_end_ppt: float = 228.0
    cfc12_start_ppt: float = 510.0
    cfc12_end_ppt: float = 507.0


@dataclass(frozen=True)
class ForcingProvenance:
    """Data source summary for the generated forcing arrays."""

    mode: str
    loaded_files: tuple[str, ...]
    missing_inputs: tuple[str, ...]
    notes: tuple[str, ...]


@dataclass(frozen=True)
class RunArtifacts:
    """Paths and high-level diagnostics from one runner invocation."""

    output_dir: Path
    provenance: ForcingProvenance
    final_cfc11_discrepancy_gg: float
    final_cfc12_discrepancy_gg: float
    final_methane_lifetime_distortion_percent: float


def _time_grid(experiment: Experiment) -> np.ndarray:
    if experiment.duration_days <= 0.0:
        raise ValueError("duration_days must be positive")
    if experiment.output_interval_days <= 0.0:
        raise ValueError("output_interval_days must be positive")

    duration_s = experiment.duration_days * DAY_S
    interval_s = experiment.output_interval_days * DAY_S
    time_s = np.arange(0.0, duration_s + interval_s * 0.5, interval_s)
    if time_s[-1] < duration_s:
        time_s = np.append(time_s, duration_s)
    return time_s


def _first_csv(directory: Path) -> Path | None:
    if not directory.exists():
        return None
    return next(iter(sorted(directory.glob("*.csv"))), None)


def _read_csv_series(path: Path, value_columns: Sequence[str]) -> tuple[np.ndarray, np.ndarray]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or ()
        if "time_s" in fieldnames:
            time_column = "time_s"
            time_scale = 1.0
        elif "day" in fieldnames:
            time_column = "day"
            time_scale = DAY_S
        else:
            raise ValueError(f"{path} must contain either time_s or day")

        value_column = next((name for name in value_columns if name in fieldnames), None)
        if value_column is None:
            columns = ", ".join(value_columns)
            raise ValueError(f"{path} must contain one of: {columns}")

        times: list[float] = []
        values: list[float] = []
        for row in reader:
            times.append(float(row[time_column]) * time_scale)
            values.append(float(row[value_column]))

    if len(times) < 2:
        raise ValueError(f"{path} must contain at least two records")
    time_array = np.asarray(times, dtype=float)
    value_array = np.asarray(values, dtype=float)
    if np.any(np.diff(time_array) <= 0.0):
        raise ValueError(f"{path} time values must be strictly increasing")
    if not np.isfinite(value_array).all():
        raise ValueError(f"{path} contains non-finite values")
    return time_array, value_array


def _interpolate_csv(
    path: Path,
    value_columns: Sequence[str],
    target_time_s: np.ndarray,
) -> np.ndarray:
    source_time_s, source_values = _read_csv_series(path, value_columns)
    return np.interp(target_time_s, source_time_s, source_values)


def build_forcing(
    experiment: Experiment,
    data_dir: Path = Path("data"),
) -> tuple[CoupledForcing, ForcingProvenance]:
    """Build forcing arrays, using CSV inputs only when all required files exist."""
    time_s = _time_grid(experiment)
    required_paths = {
        "ph": _first_csv(data_dir / "ph"),
        "chlorophyll": _first_csv(data_dir / "chlorophyll"),
        "salinity": _first_csv(data_dir / "salinity"),
    }
    missing_inputs = tuple(
        name for name, path in required_paths.items() if path is None
    )

    if missing_inputs:
        provenance = ForcingProvenance(
            mode="deterministic_validation_defaults",
            loaded_files=(),
            missing_inputs=missing_inputs,
            notes=(
                "Using static validation forcing because one or more forcing "
                "CSV files are missing.",
            ),
        )
        ph = np.full(time_s.size, experiment.ph)
        chlorophyll = np.full(time_s.size, experiment.chlorophyll_a_mg_m3)
        salinity = np.full(time_s.size, experiment.salinity_psu)
        gradient_richardson_number = np.full(
            time_s.size,
            experiment.gradient_richardson_number,
        )
    else:
        assert required_paths["ph"] is not None
        assert required_paths["chlorophyll"] is not None
        assert required_paths["salinity"] is not None
        ph = _interpolate_csv(required_paths["ph"], ("ph",), time_s)
        chlorophyll = _interpolate_csv(
            required_paths["chlorophyll"],
            ("chlorophyll_a_mg_m3", "chlorophyll", "chl_a_mg_m3"),
            time_s,
        )
        salinity = _interpolate_csv(
            required_paths["salinity"],
            ("salinity_psu", "salinity"),
            time_s,
        )
        try:
            gradient_richardson_number = _interpolate_csv(
                required_paths["salinity"],
                (
                    "gradient_richardson_number",
                    "richardson_number",
                    "rig",
                ),
                time_s,
            )
            notes = ()
        except ValueError:
            gradient_richardson_number = np.full(
                time_s.size,
                experiment.gradient_richardson_number,
            )
            notes = (
                "Salinity CSV did not include a Richardson-number column; "
                "using the experiment default.",
            )

        provenance = ForcingProvenance(
            mode="csv_forcing",
            loaded_files=tuple(
                str(path)
                for path in required_paths.values()
                if path is not None
            ),
            missing_inputs=(),
            notes=notes,
        )

    forcing = CoupledForcing(
        time_s=time_s,
        sst_k=np.full(time_s.size, experiment.sst_k),
        air_temperature_k=np.full(time_s.size, experiment.air_temperature_k),
        wind_speed_m_s=np.full(time_s.size, experiment.wind_speed_m_s),
        salinity_psu=salinity,
        gradient_richardson_number=gradient_richardson_number,
        ph=ph,
        chlorophyll_a_mg_m3=chlorophyll,
        cfc11_ppt=np.linspace(
            experiment.cfc11_start_ppt,
            experiment.cfc11_end_ppt,
            time_s.size,
        ),
        cfc12_ppt=np.linspace(
            experiment.cfc12_start_ppt,
            experiment.cfc12_end_ppt,
            time_s.size,
        ),
    )
    return forcing, provenance


def climate_accounting_error_matrix(result: UnifiedResult) -> np.ndarray:
    """Join the CFC and methane diagnostics into the Paper 3 B matrix."""
    cfc = result.cfc_inversion_discrepancy_matrix
    methane = result.methane_lifetime_distortion_matrix
    if cfc.shape[0] != methane.shape[0]:
        raise ValueError("diagnostic matrices must have the same row count")
    if not np.allclose(cfc[:, 0], methane[:, 0], rtol=0.0, atol=1.0e-9):
        raise ValueError("diagnostic matrices must share identical time_s values")

    time_s = cfc[:, 0]
    delta_t_s = np.concatenate(([0.0], np.diff(time_s)))
    return np.column_stack(
        (
            time_s,
            time_s / DAY_S,
            delta_t_s,
            delta_t_s / DAY_S,
            cfc[:, 1],
            cfc[:, 2],
            cfc[:, 3],
            cfc[:, 4],
            methane[:, 3],
            methane[:, 4],
            methane[:, 5],
        )
    )


def _write_csv(path: Path, columns: Sequence[str], rows: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        writer.writerows(rows)


def write_state_trajectory(output_dir: Path, result: UnifiedResult) -> None:
    rows = np.column_stack((result.time_s, result.time_s / DAY_S, result.state.T))
    _write_csv(
        output_dir / "unified_state_trajectory.csv",
        ("time_s", "day", *STATE_NAMES),
        rows,
    )


def write_parameters(
    output_dir: Path,
    config: UnifiedConfig,
    experiment: Experiment,
    provenance: ForcingProvenance,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[tuple[str, str, object, str]] = []
    for field in fields(experiment):
        rows.append(("experiment", field.name, getattr(experiment, field.name), ""))
    for field in fields(config):
        rows.append(("config", field.name, getattr(config, field.name), ""))
    rows.extend(
        [
            ("derived", "mixed_layer_volume_m3", config.mixed_layer_volume_m3, ""),
            ("derived", "deep_ocean_volume_m3", config.deep_ocean_volume_m3, ""),
            ("forcing", "mode", provenance.mode, ""),
            (
                "forcing",
                "loaded_files",
                "; ".join(provenance.loaded_files),
                "",
            ),
            (
                "forcing",
                "missing_inputs",
                "; ".join(provenance.missing_inputs),
                "",
            ),
            ("forcing", "notes", " ".join(provenance.notes), ""),
            ("species", "CFC11_molecular_weight_g_mol", CFC11.molecular_weight_g_mol, ""),
            ("species", "CFC12_molecular_weight_g_mol", CFC12.molecular_weight_g_mol, ""),
        ]
    )

    with (output_dir / "parameters.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("section", "parameter", "value", "unit"))
        writer.writerows(rows)


def write_summary(
    output_dir: Path,
    experiment: Experiment,
    provenance: ForcingProvenance,
    b_matrix: np.ndarray,
) -> None:
    final = b_matrix[-1]
    cfc11_sum = float(np.sum(b_matrix[:, 6]))
    cfc12_sum = float(np.sum(b_matrix[:, 7]))
    lines = [
        "# Coupled Engine Validation Summary",
        "",
        f"Experiment: `{experiment.name}`",
        f"Duration: {experiment.duration_days:g} days",
        f"Output interval: {experiment.output_interval_days:g} days",
        f"Forcing mode: `{provenance.mode}`",
        "",
        "## Forcing Pre-Flight",
        "",
    ]
    if provenance.missing_inputs:
        lines.append(
            "Missing optional forcing inputs: "
            + ", ".join(provenance.missing_inputs)
        )
    if provenance.loaded_files:
        lines.append("Loaded forcing files:")
        lines.extend(f"- `{path}`" for path in provenance.loaded_files)
    if provenance.notes:
        lines.extend(provenance.notes)
    lines.extend(
        [
            "",
            "## Final Diagnostics",
            "",
            (
                "Final CFC-11 inversion discrepancy: "
                f"{final[6]:.8g} Gg"
            ),
            (
                "Final CFC-12 inversion discrepancy: "
                f"{final[7]:.8g} Gg"
            ),
            (
                "Final methane-lifetime distortion: "
                f"{final[10]:.8g}%"
            ),
            (
                "Reported-timestep sum, CFC-11 discrepancy column: "
                f"{cfc11_sum:.8g} Gg"
            ),
            (
                "Reported-timestep sum, CFC-12 discrepancy column: "
                f"{cfc12_sum:.8g} Gg"
            ),
            "",
            "## Interpretation Boundary",
            "",
            "The climate accounting error matrix is a reporting join over the",
            "existing CFC inversion discrepancy and methane-lifetime distortion",
            "matrices. It does not add new model physics beyond the 17-state",
            "implicit Radau integration.",
            "",
        ]
    )
    (output_dir / "summary.md").write_text("\n".join(lines), encoding="utf-8")


def write_outputs(
    output_dir: Path,
    result: UnifiedResult,
    config: UnifiedConfig,
    experiment: Experiment,
    provenance: ForcingProvenance,
) -> np.ndarray:
    """Write all Paper 3 runner artifacts and return the joined B matrix."""
    output_dir.mkdir(parents=True, exist_ok=True)
    b_matrix = climate_accounting_error_matrix(result)
    _write_csv(
        output_dir / "cfc_inversion_discrepancy_matrix.csv",
        result.cfc_inversion_columns,
        result.cfc_inversion_discrepancy_matrix,
    )
    _write_csv(
        output_dir / "methane_lifetime_distortion_matrix.csv",
        result.methane_lifetime_columns,
        result.methane_lifetime_distortion_matrix,
    )
    _write_csv(
        output_dir / "climate_accounting_error_matrix.csv",
        CLIMATE_ACCOUNTING_ERROR_COLUMNS,
        b_matrix,
    )
    write_state_trajectory(output_dir, result)
    write_parameters(output_dir, config, experiment, provenance)
    write_summary(output_dir, experiment, provenance, b_matrix)
    return b_matrix


def run_model(
    output_dir: Path = Path("results"),
    data_dir: Path = Path("data"),
    experiment: Experiment | None = None,
    config: UnifiedConfig | None = None,
) -> RunArtifacts:
    """Run the default coupled experiment and write publication artifacts."""
    experiment = experiment or Experiment()
    config = config or UnifiedConfig()
    forcing, provenance = build_forcing(experiment, data_dir)
    result = integrate_unified_state_vector(forcing, config=config)
    b_matrix = write_outputs(output_dir, result, config, experiment, provenance)
    final = b_matrix[-1]
    return RunArtifacts(
        output_dir=output_dir,
        provenance=provenance,
        final_cfc11_discrepancy_gg=float(final[6]),
        final_cfc12_discrepancy_gg=float(final[7]),
        final_methane_lifetime_distortion_percent=float(final[10]),
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the coupled biogeochemical Paper 3 validation model."
    )
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--days", type=float, default=60.0)
    parser.add_argument("--interval-days", type=float, default=1.0)
    parser.add_argument("--ph", type=float, default=7.7)
    parser.add_argument("--chlorophyll", type=float, default=2.0)
    parser.add_argument("--salinity", type=float, default=35.0)
    parser.add_argument("--gradient-richardson-number", type=float, default=1.0)
    args = parser.parse_args()

    experiment = Experiment(
        duration_days=args.days,
        output_interval_days=args.interval_days,
        ph=args.ph,
        chlorophyll_a_mg_m3=args.chlorophyll,
        salinity_psu=args.salinity,
        gradient_richardson_number=args.gradient_richardson_number,
    )
    artifacts = run_model(
        output_dir=args.output_dir,
        data_dir=args.data_dir,
        experiment=experiment,
    )
    if artifacts.provenance.missing_inputs:
        print(
            "Using deterministic validation forcing; missing optional inputs: "
            + ", ".join(artifacts.provenance.missing_inputs)
        )
    print(f"Results written to {artifacts.output_dir.resolve()}")
    print(
        "Final CFC discrepancies: "
        f"CFC-11={artifacts.final_cfc11_discrepancy_gg:.8g} Gg, "
        f"CFC-12={artifacts.final_cfc12_discrepancy_gg:.8g} Gg"
    )
    print(
        "Final methane-lifetime distortion: "
        f"{artifacts.final_methane_lifetime_distortion_percent:.8g}%"
    )


if __name__ == "__main__":
    main()
