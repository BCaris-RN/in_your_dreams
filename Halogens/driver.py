# ============================================================================
# Copyright (c) 2026 Brandon W. Caris. All rights reserved.
# Part of the Coupled Marine Engine Suite (in_your_dreams/Coupled_Engine)
# Licensed under CC BY 4.0 with Explicit Academic Citation Mandate.
# See master LICENSE file at repository root for full terms and string.
# ============================================================================
"""Run cold-baseline and storm-coupled marine heatwave scenarios."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from model import HalogenBoxModel, MarineState

OUTPUT_DIR = Path(__file__).parent / "outputs"
SIMULATION_DURATION_DAYS = 120.0
OUTPUT_INTERVAL_S = 3600.0


def save_scenario(model, scenario, solution):
    header = "time_s,time_h," + ",".join(
        f"{name}_molecule_cm3" for name in model.state_names
    )
    table = np.column_stack((solution.t, solution.t / 3600.0, solution.y.T))
    np.savetxt(
        OUTPUT_DIR / f"{scenario.name}.csv",
        table,
        delimiter=",",
        header=header,
        comments="",
    )


def convergence_diagnostics(model, scenario, solution):
    drifts = np.column_stack(
        [
            model.fractional_drift_per_hour(time_s, state, scenario)
            for time_s, state in zip(solution.t, solution.y.T)
        ]
    )
    meets_criterion = np.all(
        drifts <= model.steady_state_fractional_drift_per_hour,
        axis=0,
    )
    remains_converged = np.flip(
        np.logical_and.accumulate(np.flip(meets_criterion))
    )
    converged_indices = np.flatnonzero(remains_converged)
    convergence_index = (
        int(converged_indices[0]) if converged_indices.size else None
    )
    return drifts, convergence_index


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    model = HalogenBoxModel()

    cold_baseline = MarineState(
        name="cold_baseline",
        sst_k=285.15,
        air_temperature_k=285.15,
        wind_speed_m_s=5.0,
    )
    storm_heatwave = MarineState(
        name="marine_anomaly",
        sst_k=288.15,
        air_temperature_k=288.15,
        wind_speed_m_s=12.0,
    )

    scenarios = (cold_baseline, storm_heatwave)
    duration_s = SIMULATION_DURATION_DAYS * 86_400.0
    runs = {
        scenario.name: model.run(
            scenario,
            duration_s=duration_s,
            output_interval_s=OUTPUT_INTERVAL_S,
        )
        for scenario in scenarios
    }
    diagnostics = {
        scenario.name: convergence_diagnostics(
            model,
            scenario,
            runs[scenario.name],
        )
        for scenario in scenarios
    }
    for scenario in scenarios:
        save_scenario(model, scenario, runs[scenario.name])

    figure, axes = plt.subplots(4, 1, figsize=(10, 12), sharex=True)
    colors = {"cold_baseline": "royalblue", "marine_anomaly": "crimson"}
    labels = {
        "cold_baseline": "Case A: 12 C, 5 m/s",
        "marine_anomaly": "Case B: 15 C, 12 m/s",
    }

    for scenario in scenarios:
        solution = runs[scenario.name]
        days = solution.t / 86_400.0
        axes[0].plot(
            days,
            solution.y[0],
            color=colors[scenario.name],
            label=labels[scenario.name],
        )
        axes[1].plot(
            days,
            solution.y[1] + solution.y[2] + solution.y[5],
            color=colors[scenario.name],
            label=labels[scenario.name],
        )
        axes[2].plot(
            days,
            solution.y[5],
            color=colors[scenario.name],
            label=labels[scenario.name],
        )
        axes[3].plot(
            days,
            solution.y[4],
            color=colors[scenario.name],
            label=labels[scenario.name],
        )

    axes[0].set_ylabel("CHBr3\n(molecule cm-3)")
    axes[1].set_ylabel("Active Br family\n(molecule cm-3)")
    axes[2].set_ylabel("HOBr\n(molecule cm-3)")
    axes[3].set_ylabel("OH\n(molecule cm-3)")
    axes[3].set_xlabel("Simulation time (days)")
    for axis in axes:
        axis.grid(alpha=0.25)
        axis.legend()
    figure.suptitle("Terminal-sink equilibration of the bromine family")
    figure.tight_layout()
    figure.savefig(OUTPUT_DIR / "radical_trends.png", dpi=180)
    plt.close(figure)

    baseline = runs[cold_baseline.name]
    heatwave = runs[storm_heatwave.name]
    baseline_endpoint_oh = baseline.y[4, -1]
    heatwave_endpoint_oh = heatwave.y[4, -1]
    baseline_endpoint_tendency = model.tendency(
        baseline.t[-1],
        baseline.y[:, -1],
        cold_baseline,
    )
    heatwave_endpoint_tendency = model.tendency(
        heatwave.t[-1],
        heatwave.y[:, -1],
        storm_heatwave,
    )
    baseline_drift, baseline_convergence_index = diagnostics[
        cold_baseline.name
    ]
    heatwave_drift, heatwave_convergence_index = diagnostics[
        storm_heatwave.name
    ]
    baseline_endpoint_drift = baseline_drift[:, -1]
    heatwave_endpoint_drift = heatwave_drift[:, -1]
    baseline_is_steady = model.is_steady_state(
        baseline.t[-1], baseline.y[:, -1], cold_baseline
    )
    heatwave_is_steady = model.is_steady_state(
        heatwave.t[-1], heatwave.y[:, -1], storm_heatwave
    )
    endpoint_oh_shift_percent = (
        (heatwave_endpoint_oh - baseline_endpoint_oh)
        / baseline_endpoint_oh
        * 100.0
    )
    source_shift_percent = (
        (
            model.ocean_source_rate(storm_heatwave)
            - model.ocean_source_rate(cold_baseline)
        )
        / model.ocean_source_rate(cold_baseline)
        * 100.0
    )

    summary = "\n".join(
        [
            "Halogens scenario comparison",
            "============================",
            (
                "Cold baseline ocean source: "
                f"{model.ocean_source_rate(cold_baseline):.6g} "
                "molecule cm-3 s-1"
            ),
            (
                "Marine anomaly ocean source: "
                f"{model.ocean_source_rate(storm_heatwave):.6g} "
                "molecule cm-3 s-1"
            ),
            f"Ocean-source shift: {source_shift_percent:.3f}%",
            (
                "Equilibrated OH, cold baseline: "
                f"{baseline_endpoint_oh:.12g} molecule cm-3"
            ),
            (
                "Equilibrated OH, marine anomaly: "
                f"{heatwave_endpoint_oh:.12g} molecule cm-3"
            ),
            f"Equilibrated OH shift: {endpoint_oh_shift_percent:.9g}%",
            (
                "Equilibrated active Br family, cold baseline: "
                f"{baseline.y[[1, 2, 5], -1].sum():.12g} "
                "molecule cm-3"
            ),
            (
                "Equilibrated active Br family, marine anomaly: "
                f"{heatwave.y[[1, 2, 5], -1].sum():.12g} "
                "molecule cm-3"
            ),
            (
                "HOBr deposition balance, cold baseline: "
                f"{model.hobr_termination_s * baseline.y[5, -1]:.12g} "
                "deposition vs. "
                f"{3.0 * model.chbr3_photolysis_s * baseline.y[0, -1]:.12g} "
                "photolytic Br production, molecule cm-3 s-1"
            ),
            (
                "HOBr deposition balance, marine anomaly: "
                f"{model.hobr_termination_s * heatwave.y[5, -1]:.12g} "
                "deposition vs. "
                f"{3.0 * model.chbr3_photolysis_s * heatwave.y[0, -1]:.12g} "
                "photolytic Br production, molecule cm-3 s-1"
            ),
            (
                "Convergence criterion: every species fractional drift <= "
                f"{model.steady_state_fractional_drift_per_hour:.3g} per hour"
            ),
            (
                "First sustained convergence, cold baseline: "
                f"{baseline.t[baseline_convergence_index] / 86_400.0:.6g} days"
                if baseline_convergence_index is not None
                else "First sustained convergence, cold baseline: not reached"
            ),
            (
                "First sustained convergence, marine anomaly: "
                f"{heatwave.t[heatwave_convergence_index] / 86_400.0:.6g} days"
                if heatwave_convergence_index is not None
                else "First sustained convergence, marine anomaly: not reached"
            ),
            (
                "Maximum endpoint fractional drift, cold baseline: "
                f"{baseline_endpoint_drift.max():.12g} per hour "
                "("
                f"{model.state_names[int(np.argmax(baseline_endpoint_drift))]}"
                ")"
            ),
            (
                "Maximum endpoint fractional drift, marine anomaly: "
                f"{heatwave_endpoint_drift.max():.12g} per hour "
                "("
                f"{model.state_names[int(np.argmax(heatwave_endpoint_drift))]}"
                ")"
            ),
            (
                "Endpoint active Br residual, cold baseline: "
                f"{baseline_endpoint_tendency[[1, 2, 5]].sum():.12g} "
                "molecule cm-3 s-1"
            ),
            (
                "Endpoint active Br residual, marine anomaly: "
                f"{heatwave_endpoint_tendency[[1, 2, 5]].sum():.12g} "
                "molecule cm-3 s-1"
            ),
            (
                "Steady-state status, cold baseline: "
                f"{'reached' if baseline_is_steady else 'not reached'}"
            ),
            (
                "Steady-state status, marine anomaly: "
                f"{'reached' if heatwave_is_steady else 'not reached'}"
            ),
            "",
            "Note: source conversion and background HOx loss frequencies are",
            "prototype scaling assumptions, not validated atmospheric budgets.",
            "",
        ]
    )
    (OUTPUT_DIR / "summary.txt").write_text(summary, encoding="utf-8")
    print(summary)


if __name__ == "__main__":
    main()
