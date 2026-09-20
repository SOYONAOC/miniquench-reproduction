"""All radiative terms with explicit, non-author-confirmed shell closure."""

import sys
import numpy as np
from .physics import ROOT, MSUN, MP, C

sys.path.insert(0, str(ROOT / "third_party/Lyman-alpha-feedback/M_F_fit"))
from M_F_fit import M_F_fit

TABLE = np.genfromtxt(ROOT / "data/schure_cooling.csv", delimiter=",", names=True)
USE = TABLE["logT"] >= 4
LOGT = TABLE["logT"][USE]
LOGL = TABLE["logLambdahd"][USE]


def thermal_equilibrium(heating, nH, volume):
    if heating <= nH * nH * volume * 10 ** LOGL[0]:
        return 1e4
    target = np.log10(heating / (nH * nH * volume))
    candidates = np.flatnonzero((LOGL[:-1] < target) & (LOGL[1:] >= target))
    if not len(candidates):
        raise ValueError(
            "No stable thermal equilibrium inside Schure table: closure fails"
        )
    i = candidates[0]
    return 10 ** np.interp(target, LOGL[i : i + 2], LOGT[i : i + 2])


def radiation_forces(r, v, mseg, lum, cfg):
    if cfg.model != "full_approx":
        return 0.0, 0.0, float("nan"), float("nan")
    _, luv, llya = lum
    area = 4 * np.pi * r * r * cfg.cover
    sigma = mseg * MSUN / area
    dust = sigma * cfg.metallicity / 162
    uv = (-np.expm1(-1e3 * dust) + 10**0.7 * dust) * luv / C
    N = 0.753 * sigma / MP
    n = N / (cfg.shell_thickness * r)
    volume = area * cfg.shell_thickness * r
    temp = thermal_equilibrium(0.5 * (luv + llya) * cfg.cover, n, volume)
    mult = M_F_fit(
        N,
        temp,
        cfg.metallicity,
        nHI_average=n,
        z=cfg.z,
        vmax=v / 1e5 if cfg.lya_velocity else 0.0,
    )
    return float(uv), float(mult * llya / C), float(temp), float(mult)
