"""Physical and numerical parameters for the box model."""

from dataclasses import dataclass


@dataclass(frozen=True)
class WarnerWeissCoefficients:
    """Coefficients for Warner and Weiss (1985), mol kg-1 atm-1."""

    a1: float
    a2: float
    a3: float
    b1: float
    b2: float
    b3: float


@dataclass(frozen=True)
class Species:
    name: str
    atmospheric_column: str
    molecular_weight_g_mol: float
    coefficients: WarnerWeissCoefficients


@dataclass(frozen=True)
class ModelConfig:
    """Representative global-ocean box geometry and exchange parameters."""

    ocean_area_m2: float = 3.61e14
    mixed_layer_depth_m: float = 50.0
    deep_ocean_depth_m: float = 3650.0
    gas_transfer_velocity_m_yr: float = 30.0
    vertical_exchange_timescale_yr: float = 10.0
    atmospheric_nudging_timescale_yr: float = 1.0 / 12.0
    tropospheric_dry_air_mol: float = 1.51e20
    surface_pressure_pa: float = 101_325.0
    reference_temperature_k: float = 288.15
    salinity_psu: float = 35.0
    seawater_density_kg_m3: float = 1025.0
    initial_deep_saturation_fraction: float = 0.15
    relative_tolerance: float = 1e-8
    absolute_tolerance: float = 1e-3

    @property
    def mixed_layer_volume_m3(self) -> float:
        return self.ocean_area_m2 * self.mixed_layer_depth_m

    @property
    def deep_ocean_volume_m3(self) -> float:
        return self.ocean_area_m2 * self.deep_ocean_depth_m


CFC11 = Species(
    name="CFC-11",
    atmospheric_column="cfc11_ppt",
    molecular_weight_g_mol=137.368,
    coefficients=WarnerWeissCoefficients(
        a1=-136.2685,
        a2=206.1150,
        a3=57.2805,
        b1=-0.148598,
        b2=0.095114,
        b3=-0.0163396,
    ),
)

CFC12 = Species(
    name="CFC-12",
    atmospheric_column="cfc12_ppt",
    molecular_weight_g_mol=120.913,
    coefficients=WarnerWeissCoefficients(
        a1=-124.4395,
        a2=185.4299,
        a3=51.6383,
        b1=-0.149779,
        b2=0.094668,
        b3=-0.0160043,
    ),
)

SPECIES = (CFC11, CFC12)
