# Raw input data

Files in this directory are downloaded unchanged by:

```powershell
cfc-fetch-data
```

- `noaa_global_ocean_sst_anomaly.json`: NOAA NCEI Climate at a Glance,
  monthly global ocean temperature anomalies. The series begins January 1850
  and uses the NOAA global anomaly baseline.
- `HATS_global_F11.txt`: NOAA GML/HATS global monthly CFC-11 dry-air mole
  fractions in ppt.
- `HATS_global_F12.txt`: NOAA GML/HATS global monthly CFC-12 dry-air mole
  fractions in ppt.

The atmospheric files currently end during 2025. The forcing-preparation step
linearly extrapolates the latest 24-month trend through December 2025 and
marks every non-observed month with an `*_estimated` flag.
