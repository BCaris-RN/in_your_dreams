# ============================================================================
# Copyright (c) 2026 Brandon W. Caris. All rights reserved.
# Part of the Coupled Marine Engine Suite (in_your_dreams/Coupled_Engine)
# Licensed under CC BY 4.0 with Explicit Academic Citation Mandate.
# See master LICENSE file at repository root for full terms and string.
# ============================================================================
from pathlib import Path

from cfc_flux_model.forcing import load_noaa_hats_global


def test_noaa_hats_parser_reads_global_mean(tmp_path: Path) -> None:
    source = tmp_path / "hats.txt"
    source.write_text(
        "\n".join(
            [
                "# global monthly data",
                "HATS_F11_YYYY HATS_F11_MM HATS_NH HATS_NH_sd "
                "HATS_SH HATS_SH_sd HATS_Global HATS_Global_sd",
                "2024 11 214.0 0.1 212.0 0.2 213.0 0.3",
                "2024 12 213.8 0.1 211.8 0.2 212.8 0.3",
            ]
        ),
        encoding="utf-8",
    )

    series = load_noaa_hats_global(source, "cfc11_ppt")

    assert list(series.index.astype(str)) == ["2024-11", "2024-12"]
    assert list(series) == [213.0, 212.8]
