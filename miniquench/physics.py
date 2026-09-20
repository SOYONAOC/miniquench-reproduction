from dataclasses import dataclass, asdict
from pathlib import Path
import hashlib
import json
import numpy as np
from astropy.constants import G as _G, M_sun, R_sun, k_B, h, c, m_p

ROOT = Path(__file__).resolve().parents[1]
G = _G.cgs.value
MSUN = M_sun.cgs.value
RSUN = R_sun.cgs.value
KB = k_B.cgs.value
HP = h.cgs.value
C = c.cgs.value
MP = m_p.cgs.value
MYR = 365.25 * 86400 * 1e6
KPC = 3.0856775814913673e21
FB = 0.0489 / 0.3111


@dataclass(frozen=True)
class Config:
    mass: float = 1e8
    z: float = 6.0
    fstar: float = 0.001
    cover: float = 1.0
    gas_fraction: float = 1.0
    kappa: float = 0.035
    alpha: float = 2.35
    mchar: float = 0.2
    beta: float = 1.6
    epsilon: float = 1.0
    model: str = "full_approx"
    gravity_factor: float = 2.0
    loading: str = "outward_only"
    boundary_work: bool = True
    growth: bool = True
    cooling: bool = True
    shell_thickness: float = 0.1
    stellar_radius_power: float = 0.5
    lya_velocity: bool = False
    metallicity: float = 1e-4
    sf_window: float = 100.0
    tmax: float = 700.0
    history_dt: float = 0.1
    max_step: float = 0.5
    output_dt: float = 0.5
    rtol: float = 2e-7
    atol: float = 1e-9

    def __post_init__(self):
        for key in (
            "mass",
            "kappa",
            "mchar",
            "beta",
            "shell_thickness",
            "sf_window",
            "tmax",
            "history_dt",
            "max_step",
            "output_dt",
            "rtol",
            "atol",
        ):
            if not np.isfinite(getattr(self, key)) or getattr(self, key) <= 0:
                raise ValueError(f"{key} must be finite and positive")
        if not 0 < self.cover <= 1 or not 0 < self.gas_fraction <= 1:
            raise ValueError("cover and gas_fraction must be in (0,1]")
        if self.fstar < 0 or self.epsilon < 0 or self.z < 0:
            raise ValueError("negative physical input")
        if self.model not in ("full_approx", "sn_only", "no_feedback"):
            raise ValueError("Unknown model; exact full author closure is unavailable")
        if self.loading not in ("outward_only", "until_return"):
            raise ValueError("unknown loading convention")
        if self.tmax < self.sf_window:
            raise ValueError("tmax < sf_window")


def hubble(z):
    return 67.66 * 1e5 / (1000 * KPC) * np.sqrt(0.3111 * (1 + z) ** 3 + 0.6889)


def growth_rate(cfg):
    return 0.03 / 1000 * (1 + cfg.z) ** 2.5 if cfg.growth else 0.0


def halo(t, cfg):
    """Total mass Msun, virial radius cm, circular/inflow velocity cm/s."""
    mass = cfg.mass * np.exp(growth_rate(cfg) * np.asarray(t))
    rv = (2 * G * mass * MSUN / (200 * 0.3111 * hubble(cfg.z) ** 2)) ** (1 / 3)
    return mass, rv, np.sqrt(G * mass * MSUN / rv)


def boundary(t, cfg):
    r = cfg.kappa * halo(t, cfg)[1]
    a = growth_rate(cfg) / 3
    return r, a * r / MYR, a * a * r / MYR**2


def reservoir(t, cfg):
    """Exact prelaunch solution: gas, formed stars (Msun), SFR (Msun/Myr)."""
    t = np.asarray(t)
    a = growth_rate(cfg)
    b = cfg.fstar / (0.0026 / hubble(cfg.z) / MYR)
    g0 = cfg.gas_fraction * FB * cfg.mass
    if a + b == 0:
        gas = np.full_like(t, g0, dtype=float)
    else:
        gas = g0 * np.exp(-b * t) + FB * cfg.mass * a / (a + b) * (
            np.exp(a * t) - np.exp(-b * t)
        )
    formed = g0 + FB * cfg.mass * np.expm1(a * t) - gas
    return gas, formed, b * gas


def enclosed_mass(r, t, cfg):
    mass, rv, _ = halo(t, cfg)
    return cfg.gravity_factor * mass * r / rv


def pressure_rhs(p, r, v, luminosity, z, cooling=True):
    """CGS dp/ds, equation (16); no hidden floors."""
    tc = 120 * MYR * ((1 + z) / 10) ** -4
    return luminosity / (2 * np.pi * r**3) - 5 * p * v / r - (p / tc if cooling else 0)


def code_fingerprint():
    digest = hashlib.sha256()
    paths = sorted((ROOT / "miniquench").glob("*.py"))
    paths += sorted((ROOT / "third_party/Lyman-alpha-feedback/M_F_fit").glob("*.py"))
    paths += [ROOT / "data/schure_cooling.csv"]
    for path in paths:
        digest.update(str(path.relative_to(ROOT)).encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def cache_key(cfg):
    raw = json.dumps(
        {"config": asdict(cfg), "code": code_fingerprint()}, sort_keys=True
    )
    return hashlib.sha256(raw.encode()).hexdigest()
