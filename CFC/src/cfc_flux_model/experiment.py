# ============================================================================
# Copyright (c) 2026 Brandon W. Caris. All rights reserved.
# Part of the Coupled Marine Engine Suite (in_your_dreams/Coupled_Engine)
# Licensed under CC BY 4.0 with Explicit Academic Citation Mandate.
# See master LICENSE file at repository root for full terms and string.
# ============================================================================
"""Paired observed-SST and zero-anomaly experiments."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .forcing import Forcing
from .parameters import ModelConfig, SPECIES, Species
from .simulation import (
    SimulationResult,
    reset_diagnostics,
    run_simulation,
)


def _paired_species_run(
    forcing: Forcing,
    species: Species,
    config: ModelConfig,
    experiment_start: str,
    experiment_end: str,
    enthalpy_scale: float = 1.0,
) -> tuple[SimulationResult, SimulationResult]:
    spinup = forcing.subset(
        str(forcing.frame.index.min()),
        experiment_start,
        include_endpoint=False,
    )
    spinup_result = run_simulation(
        spinup,
        species,
        config,
        scenario="spinup",
        use_sst_anomaly=True,
        enthalpy_scale=enthalpy_scale,
    )
    branch_state = reset_diagnostics(spinup_result.final_state)
    experiment = forcing.subset(experiment_start, experiment_end)
    observed = run_simulation(
        experiment,
        species,
        config,
        scenario="observed_sst",
        initial=branch_state,
        use_sst_anomaly=True,
        enthalpy_scale=enthalpy_scale,
    )
    counterfactual = run_simulation(
        experiment,
        species,
        config,
        scenario="zero_anomaly",
        initial=branch_state,
        use_sst_anomaly=False,
        enthalpy_scale=enthalpy_scale,
    )
    return observed, counterfactual


def _mass_gigagrams(moles: float, species: Species) -> float:
    return moles * species.molecular_weight_g_mol / 1e9


def summarize_pair(
    observed: SimulationResult,
    counterfactual: SimulationResult,
) -> dict[str, float | str]:
    species = observed.species
    observed_net = -float(observed.frame["cumulative_air_to_sea_mol"].iloc[-1])
    counterfactual_net = -float(
        counterfactual.frame["cumulative_air_to_sea_mol"].iloc[-1]
    )
    observed_gross = float(observed.frame["cumulative_outgassing_mol"].iloc[-1])
    counterfactual_gross = float(
        counterfactual.frame["cumulative_outgassing_mol"].iloc[-1]
    )
    accelerated_net = observed_net - counterfactual_net
    accelerated_gross = observed_gross - counterfactual_gross
    observed_uptake = max(-observed_net, 0.0)
    counterfactual_uptake = max(-counterfactual_net, 0.0)
    uptake_reduction = counterfactual_uptake - observed_uptake
    uptake_reduction_percent = (
        uptake_reduction / counterfactual_uptake * 100.0
        if counterfactual_uptake > 0
        else np.nan
    )
    percent = (
        accelerated_gross / counterfactual_gross * 100.0
        if counterfactual_gross > 0
        else np.nan
    )
    return {
        "species": species.name,
        "reference_henry_mol_m3_pa": observed.reference_henry_mol_m3_pa,
        "effective_dissolution_enthalpy_kj_mol": (
            observed.dissolution_enthalpy_j_mol / 1000.0
        ),
        "observed_net_ocean_to_air_mol": observed_net,
        "counterfactual_net_ocean_to_air_mol": counterfactual_net,
        "temperature_driven_net_release_mol": accelerated_net,
        "temperature_driven_net_release_Gg": _mass_gigagrams(accelerated_net, species),
        "observed_ocean_uptake_Gg": _mass_gigagrams(observed_uptake, species),
        "counterfactual_ocean_uptake_Gg": _mass_gigagrams(
            counterfactual_uptake, species
        ),
        "temperature_driven_uptake_reduction_Gg": _mass_gigagrams(
            uptake_reduction, species
        ),
        "ocean_uptake_reduction_percent": uptake_reduction_percent,
        "observed_gross_outgassing_mol": observed_gross,
        "counterfactual_gross_outgassing_mol": counterfactual_gross,
        "temperature_driven_gross_outgassing_mol": accelerated_gross,
        "temperature_driven_gross_outgassing_Gg": _mass_gigagrams(
            accelerated_gross, species
        ),
        "gross_outgassing_acceleration_percent": percent,
        "observed_global_flux_regime": (
            "net source" if observed_net > 0 else "net sink"
        ),
        "max_abs_mass_balance_residual_mol": max(
            observed.frame["mass_balance_residual_mol"].abs().max(),
            counterfactual.frame["mass_balance_residual_mol"].abs().max(),
        ),
    }


def annual_fluxes(results: list[SimulationResult]) -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    for result in results:
        frame = result.frame.copy()
        cumulative_net = -frame["cumulative_air_to_sea_mol"]
        cumulative_gross = frame["cumulative_outgassing_mol"]
        interval = pd.DataFrame(
            {
                "month": pd.PeriodIndex(frame["month"].iloc[:-1], freq="M"),
                "net_ocean_to_air_mol": np.diff(cumulative_net),
                "gross_outgassing_mol": np.diff(cumulative_gross),
            }
        )
        interval["year"] = interval["month"].dt.year
        grouped = interval.groupby("year", as_index=False)[
            ["net_ocean_to_air_mol", "gross_outgassing_mol"]
        ].sum()
        for record in grouped.to_dict(orient="records"):
            rows.append(
                {
                    "species": result.species.name,
                    "scenario": result.scenario,
                    **record,
                    "net_ocean_to_air_Gg": _mass_gigagrams(
                        float(record["net_ocean_to_air_mol"]), result.species
                    ),
                    "gross_outgassing_Gg": _mass_gigagrams(
                        float(record["gross_outgassing_mol"]), result.species
                    ),
                }
            )
    return pd.DataFrame(rows)


def run_experiment(
    forcing: Forcing,
    config: ModelConfig | None = None,
    experiment_start: str = "2016-01",
    experiment_end: str = "2025-12",
) -> tuple[list[SimulationResult], pd.DataFrame]:
    config = config or ModelConfig()
    results: list[SimulationResult] = []
    summaries: list[dict[str, float | str]] = []
    for species in SPECIES:
        observed, counterfactual = _paired_species_run(
            forcing, species, config, experiment_start, experiment_end
        )
        results.extend([observed, counterfactual])
        summaries.append(summarize_pair(observed, counterfactual))
    return results, pd.DataFrame(summaries)


def run_sensitivity(
    forcing: Forcing,
    config: ModelConfig | None = None,
    experiment_start: str = "2016-01",
    experiment_end: str = "2025-12",
) -> pd.DataFrame:
    base = config or ModelConfig()
    cases = [
        ("baseline", base, 1.0),
        (
            "gas_transfer_-20pct",
            replace(
                base,
                gas_transfer_velocity_m_yr=base.gas_transfer_velocity_m_yr * 0.8,
            ),
            1.0,
        ),
        (
            "gas_transfer_+20pct",
            replace(
                base,
                gas_transfer_velocity_m_yr=base.gas_transfer_velocity_m_yr * 1.2,
            ),
            1.0,
        ),
        (
            "mixed_layer_depth_-20pct",
            replace(base, mixed_layer_depth_m=base.mixed_layer_depth_m * 0.8),
            1.0,
        ),
        (
            "mixed_layer_depth_+20pct",
            replace(base, mixed_layer_depth_m=base.mixed_layer_depth_m * 1.2),
            1.0,
        ),
        ("enthalpy_-20pct", base, 0.8),
        ("enthalpy_+20pct", base, 1.2),
    ]
    rows: list[dict[str, float | str]] = []
    for case_name, case_config, enthalpy_scale in cases:
        for species in SPECIES:
            observed, counterfactual = _paired_species_run(
                forcing,
                species,
                case_config,
                experiment_start,
                experiment_end,
                enthalpy_scale,
            )
            row = summarize_pair(observed, counterfactual)
            row["case"] = case_name
            rows.append(row)
    return pd.DataFrame(rows)


def write_outputs(
    output_dir: Path,
    results: list[SimulationResult],
    summary: pd.DataFrame,
    sensitivity: pd.DataFrame | None,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    monthly = pd.concat(
        [result.frame for result in results],
        ignore_index=True,
    )
    monthly.to_csv(output_dir / "monthly_box_states.csv", index=False)
    annual_fluxes(results).to_csv(output_dir / "annual_fluxes.csv", index=False)
    summary.to_csv(output_dir / "summary.csv", index=False)
    if sensitivity is not None:
        sensitivity.to_csv(output_dir / "sensitivity.csv", index=False)

    config = results[0].config
    parameter_rows = [
        ("ocean_area", config.ocean_area_m2, "m2", "representative global ocean"),
        (
            "mixed_layer_depth",
            config.mixed_layer_depth_m,
            "m",
            "baseline model assumption",
        ),
        (
            "deep_ocean_depth",
            config.deep_ocean_depth_m,
            "m",
            "representative deep reservoir",
        ),
        (
            "gas_transfer_velocity",
            config.gas_transfer_velocity_m_yr,
            "m yr-1",
            "baseline model assumption",
        ),
        (
            "vertical_exchange_timescale",
            config.vertical_exchange_timescale_yr,
            "yr",
            "baseline model assumption",
        ),
        (
            "atmospheric_nudging_timescale",
            config.atmospheric_nudging_timescale_yr,
            "yr",
            "keeps atmosphere tied to NOAA history",
        ),
        (
            "tropospheric_dry_air",
            config.tropospheric_dry_air_mol,
            "mol",
            "representative tropospheric inventory",
        ),
        (
            "reference_temperature",
            config.reference_temperature_k,
            "K",
            "Van 't Hoff reference",
        ),
        ("salinity", config.salinity_psu, "psu", "reference seawater"),
        (
            "seawater_density",
            config.seawater_density_kg_m3,
            "kg m-3",
            "reference seawater",
        ),
    ]
    for result in results[::2]:
        parameter_rows.extend(
            [
                (
                    f"{result.species.name}_reference_Hcp",
                    result.reference_henry_mol_m3_pa,
                    "mol m-3 Pa-1",
                    "Warner-Weiss at reference T and S",
                ),
                (
                    f"{result.species.name}_effective_dissolution_enthalpy",
                    result.dissolution_enthalpy_j_mol / 1000.0,
                    "kJ mol-1",
                    "local derivative of Warner-Weiss fit",
                ),
            ]
        )
    pd.DataFrame(
        parameter_rows,
        columns=["parameter", "value", "unit", "basis"],
    ).to_csv(output_dir / "parameters.csv", index=False)

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    for axis, species in zip(axes, SPECIES, strict=True):
        selected = monthly.loc[monthly["species"] == species.name].copy()
        for scenario, group in selected.groupby("scenario"):
            dates = pd.PeriodIndex(group["month"], freq="M").to_timestamp()
            cumulative = -group["cumulative_air_to_sea_mol"].to_numpy()
            axis.plot(
                dates,
                [_mass_gigagrams(value, species) for value in cumulative],
                label=scenario.replace("_", " "),
            )
        axis.set_title(species.name)
        axis.set_ylabel("Net ocean-to-air flux (Gg)")
        axis.grid(alpha=0.25)
        axis.legend()
    axes[-1].set_xlabel("Month")
    fig.suptitle("Cumulative CFC flux: observed SST vs zero anomaly")
    fig.tight_layout()
    fig.savefig(output_dir / "cumulative_flux.png", dpi=180)
    plt.close(fig)

    report_lines = [
        "# Temperature-Driven CFC Flux Shift",
        "",
        "Paired global box-model experiment for 2016-2025. Positive reported",
        "net release means that observed SST anomalies shifted CFC toward the",
        "atmosphere relative to the zero-anomaly counterfactual.",
        "",
    ]
    for row in summary.to_dict(orient="records"):
        report_lines.extend(
            [
                f"## {row['species']}",
                "",
                (
                    "Temperature-driven flux shift toward the atmosphere: "
                    f"{row['temperature_driven_net_release_Gg']:.4f} Gg"
                ),
                (
                    "Reduction in ocean uptake: "
                    f"{row['ocean_uptake_reduction_percent']:.2f}%"
                ),
                (
                    "Global flux regime with observed SST: "
                    f"{row['observed_global_flux_regime']}"
                ),
                "",
            ]
        )
    if sensitivity is not None:
        report_lines.extend(["## Sensitivity range", ""])
        grouped = sensitivity.groupby("species")[
            "temperature_driven_net_release_Gg"
        ].agg(["min", "max"])
        for species, row in grouped.iterrows():
            report_lines.append(f"{species}: {row['min']:.4f} to {row['max']:.4f} Gg")
        report_lines.append("")
    report_lines.extend(
        [
            "## Interpretation boundary",
            "",
            "This is a temperature-only attribution experiment in a global",
            "representative box model, not a spatial ocean inversion. Atmospheric",
            "values after the final NOAA observation are trend-extrapolated and",
            "flagged in `data/processed/monthly_forcing.csv`. A positive warming",
            "effect does not by itself mean the global ocean became a net source;",
            "it can mean that warming weakened an otherwise positive ocean sink.",
            "",
        ]
    )
    (output_dir / "RESULTS.md").write_text(
        "\n".join(report_lines),
        encoding="utf-8",
    )
