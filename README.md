# In Your Dreams

This repository contains **three Python equation-based applications** for
atmospheric and climate research:

1. **CFC air-sea flux model**: attributes changes in CFC-11 and CFC-12 ocean
   uptake to observed sea-surface temperature anomalies.
2. **Marine halogen box model**: follows CHBr3 oxidation and reactive bromine
   chemistry, including an explicit HOBr deposition and washout sink.
3. **Coupled biogeochemical engine**: combines both systems with pH,
   chlorophyll-a, and haline-stratification forcing for Paper 3.

The applications address different scientific questions and use separate
environments, assumptions, state vectors, and numerical integrations.

| Application | Scientific question | Numerical system |
| --- | --- | --- |
| [`CFC/`](CFC/) | How did 2016-2025 SST anomalies alter modeled CFC-11 and CFC-12 air-sea exchange? | Three-reservoir atmosphere, mixed-layer, and deep-ocean ODE model |
| [`Halogens/`](Halogens/) | How does continuous marine CHBr3 forcing affect reactive bromine and HOx, and when does deposition permit equilibrium? | Six-species stiff chemical box model |
| [`Coupled_Engine/`](Coupled_Engine/) | How do pH, biology, and salinity stratification jointly alter halogen chemistry and CFC inversion diagnostics? | Unified 17-state implicit Radau integration |

## Application 1: CFC air-sea flux

### Purpose

The CFC application is a temperature-attribution experiment. It compares an
observed-SST simulation with a counterfactual simulation in which SST
anomalies are set to zero. Atmospheric histories, initial inventories,
exchange parameters, salinity, and numerical settings are otherwise
identical.

The model uses the solubility-form Henry coefficient:

$$
H_{cp} = \frac{C_{\mathrm{aq}}}{p}.
$$

Its local temperature response is represented with a Van 't Hoff correction:

$$
\ln\left(\frac{H_{cp}(T)}{H_{cp}(T_{\mathrm{ref}})}\right)
= -\frac{\Delta H_{\mathrm{sol}}}{R}
\left(\frac{1}{T}-\frac{1}{T_{\mathrm{ref}}}\right).
$$

Air-sea and mixed-layer-to-deep-ocean exchange are:

$$
F_{\mathrm{air\to sea}}
= k_{\mathrm{gas}}A
\left[H_{cp}(T,S)p_{\mathrm{air}}-C_{\mathrm{mixed}}\right],
$$

$$
F_{\mathrm{mixed\to deep}}
= \frac{V_{\mathrm{mixed}}}{\tau_{\mathrm{exchange}}}
\left(C_{\mathrm{mixed}}-C_{\mathrm{deep}}\right).
$$

The physical inventories obey:

$$
\frac{dN_{\mathrm{air}}}{dt}
= -F_{\mathrm{air\to sea}} + F_{\mathrm{external}},
$$

$$
\frac{dN_{\mathrm{mixed}}}{dt}
= F_{\mathrm{air\to sea}}-F_{\mathrm{mixed\to deep}},
$$

$$
\frac{dN_{\mathrm{deep}}}{dt}
= F_{\mathrm{mixed\to deep}}.
$$

Positive air-to-sea flux is ocean uptake. A positive reported
temperature-driven release means warming shifted flux toward the atmosphere
relative to the zero-anomaly case; it does not necessarily mean the global
ocean became a net source.

See [`CFC/README.md`](CFC/README.md) for inputs, results, and reproduction.

## Application 2: Marine halogen chemistry

### Purpose

The Halogens application is a conceptual marine boundary-layer kinetic box
model. It integrates:

$$
\left[\mathrm{CHBr_3},\ \mathrm{Br},\ \mathrm{BrO},\
\mathrm{HO_2},\ \mathrm{OH},\ \mathrm{HOBr}\right].
$$

Define:

$$
L_{\mathrm{CHBr_3}}
= \left(J_{\mathrm{CHBr_3}}
+ k_{\mathrm{OH+CHBr_3}}[\mathrm{OH}]
+ k_{\mathrm{Cl+CHBr_3}}[\mathrm{Cl}]\right)[\mathrm{CHBr_3}],
$$

$$
R_{\mathrm{BrO}} = k_{\mathrm{Br+O_3}}[\mathrm{Br}][\mathrm{O_3}],
$$

$$
R_{\mathrm{HOBr}} = k_{\mathrm{BrO+HO_2}}[\mathrm{BrO}][\mathrm{HO_2}],
$$

$$
P_{\mathrm{HOBr}} = J_{\mathrm{HOBr}}[\mathrm{HOBr}],
\qquad
D_{\mathrm{HOBr}} = k_{\mathrm{term}}[\mathrm{HOBr}].
$$

The integrated equations are:

$$
\frac{d[\mathrm{CHBr_3}]}{dt}
= S_{\mathrm{ocean}}-L_{\mathrm{CHBr_3}},
$$

$$
\frac{d[\mathrm{Br}]}{dt}
= 3J_{\mathrm{CHBr_3}}[\mathrm{CHBr_3}]
+P_{\mathrm{HOBr}}-R_{\mathrm{BrO}},
$$

$$
\frac{d[\mathrm{BrO}]}{dt}
= R_{\mathrm{BrO}}-R_{\mathrm{HOBr}},
$$

$$
\frac{d[\mathrm{HO_2}]}{dt}
= P_{\mathrm{HOx}}-R_{\mathrm{HOBr}}-k_{\mathrm{HO_2}}[\mathrm{HO_2}],
$$

$$
\frac{d[\mathrm{OH}]}{dt}
= P_{\mathrm{HOx}}+P_{\mathrm{HOBr}}-k_{\mathrm{OH}}[\mathrm{OH}],
$$

$$
\frac{d[\mathrm{HOBr}]}{dt}
= R_{\mathrm{HOBr}}-P_{\mathrm{HOBr}}-D_{\mathrm{HOBr}}.
$$

The active gas-phase bromine family closes as:

$$
\frac{d([\mathrm{Br}]+[\mathrm{BrO}]+[\mathrm{HOBr}])}{dt}
= 3J_{\mathrm{CHBr_3}}[\mathrm{CHBr_3}]
-k_{\mathrm{term}}[\mathrm{HOBr}].
$$

The current experiment uses
$k_{\mathrm{term}}=1.0\times10^{-4}\ \mathrm{s^{-1}}$, corresponding to an
HOBr deposition and washout lifetime of approximately 2.78 hours.

See [`Halogens/README.md`](Halogens/README.md) for convergence results,
assumptions, and reproduction.

## Application 3: Coupled biogeochemical engine

The Paper 3 engine preserves the stiff halogen mechanism and Warner-Weiss CFC
solubility equations while adding three time-dependent indices:

- pH-scaled aqueous CHBr3 hydrolysis
- chlorophyll-a-scaled CHBr3 production
- exponential haline damping of the Wanninkhof transfer velocity

It reports matched-counterfactual matrices for stratification-driven CFC
inversion discrepancy in Gg and OH-driven methane-lifetime distortion.
See [`Coupled_Engine/README.md`](Coupled_Engine/README.md) for the state vector,
units, and output conventions.

## Reproducibility

Create a separate virtual environment in each application directory.

```powershell
cd CFC
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -e ".[dev]"
cfc-fetch-data
cfc-run-model
pytest
```

```powershell
cd Halogens
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -e ".[dev]"
halogens-run
pytest
```

## Interpretation boundary

These are transparent research prototypes, not operational atmospheric
forecast or inversion systems.

- The CFC model is a global representative-box attribution experiment. It
  does not resolve circulation, regional mixed-layer depth, wind, sea ice, or
  spatial tracer inventories.
- The Halogens model uses a conceptual ocean-source scaling and prescribed
  background HOx losses. Its abundances should not be interpreted as
  validated marine boundary-layer predictions.

Each application documents its scientific sources and known limitations in
its own directory.
