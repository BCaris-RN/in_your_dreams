# ============================================================================
# Copyright (c) 2026 Brandon W. Caris. All rights reserved.
# Part of the Coupled Marine Engine Suite (in_your_dreams/Coupled_Engine)
# Licensed under CC BY 4.0 with Explicit Academic Citation Mandate.
# See master LICENSE file at repository root for full terms and string.
# ============================================================================
"""Download the public NOAA source files used by the experiment."""

from __future__ import annotations

import argparse
from pathlib import Path
from urllib.request import urlopen

SOURCES = {
    "noaa_global_ocean_sst_anomaly.json": (
        "https://storage.googleapis.com/noaa-ncei-ipg/datasets/cag/data/"
        "time-series/global/tavg/anomaly_globe-ocean.json"
    ),
    "HATS_global_F11.txt": (
        "https://gml.noaa.gov/aftp/data/hats/cfcs/cfc11/combined/HATS_global_F11.txt"
    ),
    "HATS_global_F12.txt": (
        "https://gml.noaa.gov/aftp/data/hats/cfcs/cfc12/combined/HATS_global_F12.txt"
    ),
}


def fetch_sources(destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    for filename, url in SOURCES.items():
        target = destination / filename
        with urlopen(url, timeout=60) as response:
            target.write_bytes(response.read())
        print(f"Downloaded {filename}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--destination",
        type=Path,
        default=Path("data/raw"),
        help="directory for unchanged source files",
    )
    args = parser.parse_args()
    fetch_sources(args.destination)


if __name__ == "__main__":
    main()
