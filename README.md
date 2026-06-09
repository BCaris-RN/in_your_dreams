# In Your Dreams

Python project for a non-equilibrium, multi-box kinetic model of CFC-11 and
CFC-12 air-sea exchange between the troposphere and a warming marine mixed
layer. It combines Warner-Weiss solubility equations, historical NOAA climate
forcing, paired counterfactual simulations, and parametric sensitivity
analyses.

## Scientific question

How much did observed sea-surface temperature anomalies shift oceanic CFC
flux toward the atmosphere during the last ten complete calendar years
(2016-01-01 through 2025-12-31), holding all non-temperature drivers fixed?

## Core experiment

Run two otherwise identical simulations:

1. **Observed-SST case:** temperature-dependent Henry coefficient forced by
   the observed SST anomaly time series.
2. **Counterfactual case:** the same model with SST anomalies set to zero.

The temperature-only effect is the difference in cumulative ocean-to-
atmosphere CFC flux between those runs. In the nominal global-box result,
warming weakens continued ocean uptake; it does not make the global ocean a
net CFC source during this decade.

## Planned model boxes

- Troposphere
- Marine surface mixed layer
- Subsurface/deep-ocean reservoir

Air-sea exchange will be kinetic rather than instantaneously equilibrated.
The Henry coefficient for each CFC will be adjusted with a Van 't Hoff
relationship. The implementation must state the Henry-law convention and flux
sign convention explicitly before literature parameters are entered.

## Layout

```text
data/raw/          External source data, kept unchanged
data/processed/    Reproducible model-ready inputs
notebooks/         Exploration and validation
src/cfc_flux_model Model implementation
tests/             Regression and scientific-invariant tests
PROJECT_BRIEF.md   Scope, equations, assumptions, and milestones
```

## Reproduce the experiment

```powershell
cd G:\in_your_dreams
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -e ".[dev]"
cfc-fetch-data
cfc-run-model
pytest
```

Generated artifacts are written to `results/`:

- `RESULTS.md`: concise interpretation
- `summary.csv`: nominal paired-experiment result
- `sensitivity.csv`: one-at-a-time +/-20% parameter tests
- `parameters.csv`: model parameter and unit table
- `annual_fluxes.csv`: annual scenario fluxes
- `monthly_box_states.csv`: complete monthly state and diagnostic history
- `cumulative_flux.png`: observed-SST and zero-anomaly comparison

## Nominal result

For 2016-2025 SST anomalies averaging 0.760 °C above the NOAA baseline:

- CFC-11 flux shifted 0.386 Gg toward the atmosphere, reducing modeled ocean
  uptake by 8.59%.
- CFC-12 flux shifted 0.181 Gg toward the atmosphere, reducing modeled ocean
  uptake by 7.39%.

Across the configured sensitivity tests, the shifts are 0.310-0.461 Gg for
CFC-11 and 0.145-0.216 Gg for CFC-12.

These are transparent global representative-box estimates, not a spatial
ocean inversion. Seven atmospheric forcing months in 2025 are trend-
extrapolated and flagged in `data/processed/monthly_forcing.csv`.
