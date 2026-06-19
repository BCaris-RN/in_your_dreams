# Coupled Biogeochemical Engine

Paper 3 combines the stiff marine halogen mechanism from `Halogens/model.py`
with the time-dependent Warner-Weiss CFC solubility model from `CFC/src`.
The model is a conceptual box experiment, not a spatial ocean inversion.

## Unified state

`src/engine/__init__.py` integrates one 17-element state vector with SciPy's implicit
Radau solver:

- CHBr3, Br, BrO, HO2, OH, and HOBr in the marine boundary layer
- dissolved mixed-layer CHBr3
- atmospheric, mixed-layer, and deep-ocean CFC-11 inventories
- atmospheric, mixed-layer, and deep-ocean CFC-12 inventories
- cumulative air-to-sea and external CFC fluxes for both species

Time is in seconds. Gas chemistry uses `molecule cm-3`, aqueous CHBr3 uses
`mol m-3`, and CFC inventories use `mol`.

## Dynamic Paper 3 indices

1. **pH hydrolysis:** the bimolecular coefficient is
   `1.23e17 * exp(-107300 / (R*T))`. Liquid `[OH-]` is calculated from
   `Kw(T) / [H+]`, and their product is the first-order aqueous CHBr3 loss
   rate. Moving from pH 8.1 to 7.7 suppresses this pathway by 60.2%.
2. **Chlorophyll-a production:** the extra-tropical coastal emission flux is
   `E = 1.127e5 * 2.0 * 2.5 * Chl-a` in `molecule cm-2 s-1`. It is converted
   to `mol m-3 s-1` before entering the mixed-layer CHBr3 tendency.
3. **Haline stratification:** Wanninkhof `kw` is multiplied by the Gradient
   Richardson Number scaling `(1 + gamma * Rig) ** -0.35`. The multiplier is
   fixed at 1 for `Rig <= 0` and approaches zero for strong stratification.

## Outputs

`integrate_unified_state_vector` returns the full trajectory plus two
time-indexed matrices:

- `cfc_inversion_discrepancy_matrix`: CFC-11 and CFC-12 cumulative uptake
  differences between unstratified and stratified runs, converted to Gg.
- `methane_lifetime_distortion_matrix`: the percentage change in the local
  methane-OH lifetime proxy relative to a no-pH/no-biology chemistry run.

A negative methane-lifetime distortion means enhanced OH shortened the
chemical lifetime. The CFC discrepancy is positive when an inversion that
omits stratification would overestimate ocean uptake.

The command-line runner writes reproducible Paper 3 artifacts to `results/`:

- `cfc_inversion_discrepancy_matrix.csv`
- `methane_lifetime_distortion_matrix.csv`
- `climate_accounting_error_matrix.csv`
- `unified_state_trajectory.csv`
- `parameters.csv`
- `summary.md`

If `data/ph/*.csv`, `data/chlorophyll/*.csv`, or `data/salinity/*.csv` are
missing, the runner uses the deterministic 60-day validation forcing and
records that fallback in `summary.md` and `parameters.csv`.

## Run model

```powershell
cd G:\in_your_dreams\Coupled_Engine
.\.venv\Scripts\Activate.ps1
coupled-run-model
```

## Run tests

```powershell
cd G:\in_your_dreams\Coupled_Engine
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -e ".[dev]"
pytest
```
