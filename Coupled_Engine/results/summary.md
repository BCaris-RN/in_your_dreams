# Coupled Engine Validation Summary

Experiment: `deterministic_validation`
Duration: 60 days
Output interval: 1 days
Forcing mode: `deterministic_validation_defaults`

## Forcing Pre-Flight

Missing optional forcing inputs: ph, chlorophyll, salinity
Using static validation forcing because one or more forcing CSV files are missing.

## Final Diagnostics

Final CFC-11 inversion discrepancy: 0.0030091264 Gg
Final CFC-12 inversion discrepancy: 0.0019339123 Gg
Final methane-lifetime distortion: -0.12498843%
Reported-timestep sum, CFC-11 discrepancy column: 0.12962835 Gg
Reported-timestep sum, CFC-12 discrepancy column: 0.076846693 Gg

## Interpretation Boundary

The climate accounting error matrix is a reporting join over the
existing CFC inversion discrepancy and methane-lifetime distortion
matrices. It does not add new model physics beyond the 17-state
implicit Radau integration.
