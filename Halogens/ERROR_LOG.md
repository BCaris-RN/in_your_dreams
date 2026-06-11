# Model Error Log

## HAL-001: Artificial HOx loss frequencies clamp radical response

- **Status:** Resolved for the requested recalibration
- **Severity:** Critical scientific-model defect
- **Detected:** 2026-06-10
- **Affected files:** `model.py`, `driver.py`

### Summary

The prototype assigns uncalibrated pseudo-first-order background losses:

```python
ho2_background_loss_s = 1.0e5
oh_background_loss_s = 1.0e6
```

These values imply lifetimes of:

$$
\tau_{\mathrm{HO_2}} = 10^{-5}\ \mathrm{s}
$$

and

$$
\tau_{\mathrm{OH}} = 10^{-6}\ \mathrm{s}.
$$

The losses therefore dominate the radical budget and numerically pin the
steady-state concentrations near:

- $\mathrm{HO_2} \approx 10\ \mathrm{molecule\ cm^{-3}}$
- $\mathrm{OH} \approx 1\ \mathrm{molecule\ cm^{-3}}$

This prevents the modeled OH concentration from responding meaningfully to
changes in marine CHBr3 forcing.

### Reproduction

The exact diagnostic boundaries are:

- **Case A, cold baseline:** SST = 285.15 K, wind = 5.0 m/s
- **Case B, marine anomaly:** SST = 288.15 K, wind = 12.0 m/s

Observed model output:

- Case A ocean-source proxy: 320.656 molecule cm-3 s-1
- Case B ocean-source proxy: 2192.75 molecule cm-3 s-1
- Ocean-source increase: 583.832%
- Case A steady-state OH, `y[4]`: 1.00000016609 molecule cm-3
- Case B steady-state OH, `y[4]`: 1.00000017595 molecule cm-3
- OH increase: 2.85703036e-7%
- Final OH tendency in both cases: `dOH/dt = 0`

The previously reported 659% source increase belonged to an earlier,
incorrect boundary pair and must not be attributed to Cases A and B above.

### Scientific context

Measured total OH reactivity in the global marine boundary layer averages
about 1.9 s-1, corresponding to an OH lifetime near 0.53 s. This demonstrates
that the prototype value of 1.0e6 s-1 is many orders of magnitude too large
for a background marine-air loss frequency.

Water vapor should not be listed as a routine direct OH sink alongside CO and
CH4 in this background-reactivity term. The HO2 loss must also be represented
with explicit or independently calibrated chemistry; the OH reactivity value
cannot simply be copied into the HO2 equation.

### Required correction

1. Replace the OH placeholder with a documented, configurable total OH
   reactivity or explicit sink reactions.
2. Replace the HO2 placeholder with explicit reactions or a separately
   justified pseudo-first-order loss.
3. Rerun Cases A and B with identical initial states and report the radical
   steady states and convergence residuals.
4. Add regression tests that fail if configured radical lifetimes fall below
   a documented physical threshold.
5. Continue labeling the ocean-source conversion as conceptual until its
   units and boundary-layer volume conversion are derived.

### Remediation

On 2026-06-10, the configured background losses were changed to:

```python
ho2_background_loss_s = 5.0e-2
oh_background_loss_s = 4.0
```

These values correspond to configured lifetimes of 20 s for HO2 and 0.25 s
for OH. The exact one-day diagnostic now produces:

- Case A 24-hour OH endpoint: 312175.160748 molecule cm-3
- Case B 24-hour OH endpoint: 314910.058595 molecule cm-3
- 24-hour endpoint OH increase: 0.876077981%
- Mean OH increase after hour 1: 0.307279403%
- Case A final `dOH/dt`: 0.501469512 molecule cm-3 s-1
- Case B final `dOH/dt`: 0.558692310 molecule cm-3 s-1

The radical clamp is removed, but these endpoint values must not be described
as strict steady states because of HAL-002 below. The 0.05 s-1 HO2 loss remains
a user-specified proxy that requires independent scientific calibration.

### Reference

The marine OH-reactivity comparison is based on Thames et al. (2020),
"Missing OH Reactivity in the Global Marine Boundary Layer," which reports a
mean measured marine-boundary-layer OH reactivity of 1.9 s-1:

https://doi.org/10.5194/acp-20-4013-2020

## HAL-002: Bromine family has no terminal sink

- **Status:** Resolved on 2026-06-11
- **Severity:** High scientific-model defect
- **Detected:** 2026-06-10
- **Affected files:** `model.py`, `driver.py`

### Summary

Before remediation, the Br and BrO equations recycled bromine through the
HOBr pathway but did not remove bromine from the modeled family. Adding the
two tendencies gave:

$$
\frac{d([\mathrm{Br}] + [\mathrm{BrO}])}{dt}
= 3J_{\mathrm{CHBr_3}}[\mathrm{CHBr_3}],
$$

which remains positive while CHBr3 is present. The full chemical system
therefore cannot reach a strict steady state under continuous ocean forcing.

Longer integrations demonstrate the continuing drift:

- 1 day: OH endpoint increase = 0.8761%
- 7 days: OH endpoint increase = 6.2566%
- 30 days: OH endpoint increase = 8.3896%

At 30 days, BrO is still increasing in both scenarios. The driver now labels
the one-day values as 24-hour endpoints and reports `dOH/dt` rather than
claiming that they are strict steady-state concentrations.

### Required correction

1. Add documented bromine-family terminal losses such as deposition,
   heterogeneous uptake, or transport out of the modeled box.
2. Define a convergence criterion for every state variable.
3. Report a steady state only when all normalized tendencies satisfy that
   criterion.

### Remediation

The audited `model.py` now carries HOBr as an explicit sixth state. The
`BrO + HO2` reaction forms HOBr, HOBr photolysis recycles it to `Br + OH`,
and deposition/washout removes it at:

```python
k_term_HOBr = 1.0e-4  # s^-1; approximately 2.7-hour lifetime
```

The active gas-phase bromine budget is now:

$$
\frac{d([\mathrm{Br}] + [\mathrm{BrO}] + [\mathrm{HOBr}])}{dt}
= 3J_{\mathrm{CHBr_3}}[\mathrm{CHBr_3}]
- k_{\mathrm{term}}[\mathrm{HOBr}].
$$

The driver integrates Cases A and B for 120 days with hourly output. A state
is reported as converged only when every species has an absolute fractional
drift no greater than `1.0e-4` per hour and remains within that threshold for
the rest of the integration.

Verified results:

- Case A first sustained convergence: 34.9167 days
- Case B first sustained convergence: 12.75 days
- Case A equilibrated active bromine family: 8.16962164e7 molecule cm-3
- Case B equilibrated active bromine family: 7.66395939e8 molecule cm-3
- Case A equilibrated OH: 261621.045 molecule cm-3
- Case B equilibrated OH: 328766.097 molecule cm-3
- Equilibrated OH increase: 25.6650%
- Case A deposition/production balance: 929.683574/929.683573
  molecule cm-3 s-1
- Case B deposition/production balance: 6301.287784/6301.287784
  molecule cm-3 s-1

The one-day runs do not meet the convergence criterion. At 24 hours, BrO
still changes by approximately 2.33% per hour in Case A and 2.81% per hour in
Case B. The terminal sink permits a steady state, but does not produce an
instantaneous or one-day equilibrium.
