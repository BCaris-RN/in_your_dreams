# Scientific Sources

- NIST Chemistry WebBook, bromoform Henry-law data. The selected entry gives
  `kH = 1.9 mol kg-1 bar-1` and a temperature parameter of `4700 K`; this is
  approximately `1.9e-2 mol m-3 Pa-1` after unit conversion.
  https://webbook.nist.gov/cgi/inchi?ID=C75252&Mask=10
- Orkin et al. (2013), measurements of OH reactions with bromoform. The
  reported CHBr3 rate expression is
  `9.94e-13 exp(-387/T) cm3 molecule-1 s-1`.
  https://doi.org/10.1021/jp3128753
- Kambanis et al. (1997), absolute rates for Cl with bromomethanes. The
  CHBr3 expression is `0.43e-11 exp(-809/T) cm3 molecule-1 s-1`.
  https://doi.org/10.1021/jp9719671
- IUPAC atmospheric kinetics evaluations. Recent mechanism implementations
  commonly use `4.5e-12 exp(500/T)` for BrO + HO2; the submitted prototype's
  `580/T` value remains configurable for reproducibility.
  https://iupac.org/project/2017-024-1-100/
