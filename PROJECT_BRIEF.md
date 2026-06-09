# Project Brief

## Task

Write a Python script that builds a non-equilibrium, multi-box kinetic flux
model. Couple the tropospheric box directly to a warming marine mixed-layer
box using temperature-dependent Van 't Hoff corrections to Henry's Law
coefficients for CFC-11 and CFC-12.

## Goal

Quantify how much CFC outgassing was accelerated purely by sea-surface
temperature anomalies over 2016-2025.

## Proposed state variables

- Atmospheric burden or tropospheric mixing ratio for each species
- Dissolved mixed-layer inventory for each species
- Subsurface/deep-ocean inventory for each species

## Governing relationships

The model uses the solubility-form Henry coefficient,
$H_{cp} = C_{\mathrm{aq}}/p$, in units of
$\mathrm{mol\,m^{-3}\,Pa^{-1}}$.

The temperature dependence is represented by a local Van 't Hoff
approximation anchored to the Warner-Weiss seawater solubility at the
reference temperature:

$$
\ln\left(\frac{H_{cp}(T)}{H_{cp}(T_{\mathrm{ref}})}\right)
= -\frac{\Delta H_{\mathrm{sol}}}{R}
\left(\frac{1}{T}-\frac{1}{T_{\mathrm{ref}}}\right)
$$

The sign and value of $\Delta H_{\mathrm{sol}}$ use the same solubility
convention as $H_{cp}$.

The kinetic air-sea flux is:

$$
F_{\mathrm{air\to sea}}
= k_{\mathrm{gas}}A
\left[H_{cp}(T,S)p_{\mathrm{air}}-C_{\mathrm{mixed}}\right]
$$

Positive $F_{\mathrm{air\to sea}}$ is atmosphere-to-ocean uptake. Outgassing
is therefore the negative part of this flux.

Exchange between the mixed layer and deep-ocean reservoir is:

$$
F_{\mathrm{mixed\to deep}}
= \frac{V_{\mathrm{mixed}}}{\tau_{\mathrm{exchange}}}
\left(C_{\mathrm{mixed}}-C_{\mathrm{deep}}\right)
$$

Later refinements may include Schmidt-number scaling, wind forcing, sea ice,
salinity, and spatially resolved mixed-layer properties.

## Isolation strategy

The observed and counterfactual runs must use identical:

- Initial inventories
- Atmospheric boundary conditions
- Wind or gas-transfer velocity
- Mixed-layer depth and ocean exchange rates
- Salinity
- Numerical solver settings

Only the SST anomaly term may differ. Report:

- Annual and cumulative flux for each species and scenario
- Absolute temperature-driven increase in outgassing
- Percentage increase relative to the zero-anomaly counterfactual
- Sensitivity to Henry parameters, mixed-layer depth, and transfer velocity

## Inputs to source and document

- Monthly SST anomalies and their climatological baseline
- Reference Henry coefficients for CFC-11 and CFC-12
- Dissolution enthalpy or an equivalent empirical temperature relationship
- Atmospheric CFC histories
- Initial ocean inventories
- Mixed-layer depth, area, and vertical exchange assumptions
- Gas-transfer velocity parameterization

## Milestones

1. [x] Lock units, Henry convention, sign conventions, and literature parameters.
2. [x] Implement temperature correction and air-sea flux with unit tests.
3. [x] Implement the coupled box ODE system and mass-balance tests.
4. [x] Build the observed-SST and zero-anomaly forcing pipelines.
5. [x] Run the paired experiment and uncertainty analysis.
6. [x] Produce plots, a parameter table, and a reproducible results summary.

## Completed experiment

The model was run for January 2016 through December 2025 after a common
July 1978 through January 2016 spin-up. The observed and counterfactual runs
share the same initial state and all non-temperature forcing. The nominal
global result remains an ocean sink, but observed SST anomalies reduce uptake
by 8.59% for CFC-11 and 7.39% for CFC-12.
