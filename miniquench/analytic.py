"""Section 3 two-force estimates, deliberately separate from shell trajectories."""

import numpy as np
from scipy.integrate import cumulative_trapezoid
from .physics import Config, MSUN, G, halo, reservoir
from .feedback import History


def analytic_mass(time, fstar=0.01):
    t = np.asarray(time)
    factor = np.where(
        t <= 55, (t / 55) ** 2.31, (np.maximum(t - 55, 0) / 36 + 1) ** 1.5
    )
    return 4.4e11 * (fstar / 0.01) ** 1.5 * factor


def numerical_mass(time, z=6, mchar=0.2, growth=True, history_dt=0.05):
    cfg = Config(
        mass=1e8,
        fstar=0.01,
        z=z,
        mchar=mchar,
        model="sn_only",
        growth=growth,
        history_dt=history_dt,
        tmax=100,
        sf_window=100,
    )
    hist = History(cfg)
    n = cumulative_trapezoid(hist.values[0], hist.times, initial=0)
    nsn = np.interp(time, hist.times, n)
    mh, rv, _ = halo(time, cfg)
    gas = reservoir(time, cfg)[0]
    # Fp=2 N E/R0; Fg=2 G Mh Mgas/(rvir R0): radius cancels.
    ratio = nsn * 1e51 * rv / (G * mh * MSUN * gas * MSUN)
    return cfg.mass * ratio**1.5
