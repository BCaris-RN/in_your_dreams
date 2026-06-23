# ============================================================================
# Copyright (c) 2026 Brandon W. Caris. All rights reserved.
# Part of the Coupled Marine Engine Suite (in_your_dreams/Coupled_Engine)
# Licensed under CC BY 4.0 with Explicit Academic Citation Mandate.
# See master LICENSE file at repository root for full terms and string.
# ============================================================================
"""Run the complete CFC temperature-attribution experiment."""

from __future__ import annotations

import argparse
from pathlib import Path

from .experiment import run_experiment, run_sensitivity, write_outputs
from .forcing import prepare_forcing
from .parameters import ModelConfig


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the CFC-11/CFC-12 kinetic multi-box experiment."
    )
    parser.add_argument("--raw-dir", type=Path, default=Path("data/raw"))
    parser.add_argument(
        "--processed-file",
        type=Path,
        default=Path("data/processed/monthly_forcing.csv"),
    )
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    parser.add_argument("--start", default="2016-01")
    parser.add_argument("--end", default="2025-12")
    parser.add_argument("--skip-sensitivity", action="store_true")
    args = parser.parse_args()

    forcing = prepare_forcing(
        args.raw_dir,
        args.processed_file,
        start="1978-07",
        end=args.end,
    )
    config = ModelConfig()
    results, summary = run_experiment(
        forcing,
        config,
        experiment_start=args.start,
        experiment_end=args.end,
    )
    sensitivity = None
    if not args.skip_sensitivity:
        sensitivity = run_sensitivity(
            forcing,
            config,
            experiment_start=args.start,
            experiment_end=args.end,
        )
    write_outputs(args.output_dir, results, summary, sensitivity)
    print(summary.to_string(index=False))
    print(f"\nResults written to {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
