"""Prepare monthly SST and atmospheric CFC forcing."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

SST_SERIES_START = pd.Period("1850-01", freq="M")


@dataclass(frozen=True)
class Forcing:
    frame: pd.DataFrame

    def subset(self, start: str, end: str, include_endpoint: bool = True) -> "Forcing":
        start_period = pd.Period(start, freq="M")
        end_period = pd.Period(end, freq="M")
        if include_endpoint:
            end_period += 1
        selected = self.frame.loc[
            (self.frame.index >= start_period) & (self.frame.index <= end_period)
        ].copy()
        if include_endpoint and end_period not in selected.index:
            endpoint = selected.iloc[-1].copy()
            endpoint.name = end_period
            selected = pd.concat([selected, endpoint.to_frame().T])
        return Forcing(selected)


def load_noaa_ocean_anomalies(path: Path) -> pd.Series:
    """Load NOAA CAG global-ocean monthly anomalies from its raw JSON array."""
    values = np.asarray(json.loads(path.read_text(encoding="utf-8")), dtype=float)
    index = pd.period_range(SST_SERIES_START, periods=len(values), freq="M")
    return pd.Series(values, index=index, name="sst_anomaly_c")


def load_noaa_hats_global(path: Path, name: str) -> pd.Series:
    """Load the global monthly mean from a NOAA HATS combined text file."""
    records: list[tuple[pd.Period, float]] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or not line[0].isdigit():
            continue
        fields = line.split()
        if len(fields) < 8:
            continue
        year, month = int(fields[0]), int(fields[1])
        global_mean = float(fields[6])
        if np.isfinite(global_mean):
            records.append((pd.Period(year=year, month=month, freq="M"), global_mean))
    if not records:
        raise ValueError(f"no monthly NOAA HATS records found in {path}")
    return pd.Series(
        [value for _, value in records],
        index=[period for period, _ in records],
        name=name,
        dtype=float,
    ).sort_index()


def _complete_monthly_series(
    series: pd.Series,
    start: pd.Period,
    end: pd.Period,
    trend_months: int = 24,
) -> tuple[pd.Series, pd.Series]:
    """Interpolate gaps and linearly extend the latest observed trend."""
    target_index = pd.period_range(start, end, freq="M")
    completed = series.reindex(target_index)
    observed = completed.notna()
    completed = completed.interpolate(method="linear", limit_area="inside")

    last_observed = series.index.max()
    future = target_index[target_index > last_observed]
    if len(future):
        history = series.loc[:last_observed].tail(trend_months)
        x = np.arange(len(history), dtype=float)
        slope, intercept = np.polyfit(x, history.to_numpy(dtype=float), 1)
        future_x = np.arange(len(history), len(history) + len(future), dtype=float)
        completed.loc[future] = intercept + slope * future_x

    first_observed = series.index.min()
    earlier = target_index[target_index < first_observed]
    if len(earlier):
        completed.loc[earlier] = series.iloc[0]

    if completed.isna().any():
        missing = ", ".join(str(period) for period in completed[completed.isna()].index)
        raise ValueError(f"forcing contains unresolved missing months: {missing}")
    return completed, ~observed


def prepare_forcing(
    raw_dir: Path,
    output_csv: Path | None = None,
    start: str = "1978-07",
    end: str = "2025-12",
) -> Forcing:
    """Combine source files into a continuous, model-ready monthly table."""
    start_period = pd.Period(start, freq="M")
    end_period = pd.Period(end, freq="M")
    sst = load_noaa_ocean_anomalies(raw_dir / "noaa_global_ocean_sst_anomaly.json")
    cfc11 = load_noaa_hats_global(raw_dir / "HATS_global_F11.txt", "cfc11_ppt")
    cfc12 = load_noaa_hats_global(raw_dir / "HATS_global_F12.txt", "cfc12_ppt")

    sst_completed, sst_estimated = _complete_monthly_series(
        sst, start_period, end_period
    )
    cfc11_completed, cfc11_estimated = _complete_monthly_series(
        cfc11, start_period, end_period
    )
    cfc12_completed, cfc12_estimated = _complete_monthly_series(
        cfc12, start_period, end_period
    )

    frame = pd.concat(
        [sst_completed, cfc11_completed, cfc12_completed],
        axis=1,
    )
    frame["sst_estimated"] = sst_estimated
    frame["cfc11_estimated"] = cfc11_estimated
    frame["cfc12_estimated"] = cfc12_estimated
    frame.index.name = "month"

    if output_csv is not None:
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        export = frame.copy()
        export.index = export.index.astype(str)
        export.to_csv(output_csv)
    return Forcing(frame)
