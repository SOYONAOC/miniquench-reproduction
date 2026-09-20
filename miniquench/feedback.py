"""Unit-mass age kernels and causal, immutable SFH convolution."""

from functools import lru_cache
import numpy as np
from scipy.integrate import cumulative_trapezoid, quad
from .physics import C, HP, KB, RSUN, ROOT, code_fingerprint


class IMF:
    def __init__(self, alpha=2.35, mchar=0.2, beta=1.6):
        self.alpha, self.mchar, self.beta = alpha, mchar, beta
        self.norm = quad(lambda m: m * self.raw(m), 0.08, 100, epsabs=1e-11)[0]

    def raw(self, m):
        return np.asarray(m) ** (-self.alpha) * np.exp(
            -((self.mchar / np.asarray(m)) ** self.beta)
        )

    def __call__(self, m):
        return self.raw(m) / self.norm

    def sn_number(self):
        return quad(self, 8, 100, epsabs=1e-12)[0]


def lifetime(m):
    return 1e4 * np.asarray(m) ** -2.5  # Myr


def sn_kernel(age, imf):
    age = np.asarray(age, dtype=float)
    result = np.zeros_like(age)
    use = (age >= lifetime(100)) & (age <= lifetime(8))
    m = (age[use] / 1e4) ** -0.4
    result[use] = imf(m) * 0.4 * m / age[use]
    return result


class Kernels:
    def __init__(self, alpha, mchar, beta, radius_power):
        self.imf = IMF(alpha, mchar, beta)
        self.mass = np.unique(np.r_[np.geomspace(0.08, 100, 8193), 8.0])
        m = self.mass
        phi = self.imf(m)
        self.count = cumulative_trapezoid(phi, m, initial=0)
        temp = 5772 * m**0.625
        radius = RSUN * m**radius_power

        # Integrate the blackbody spectrum in dimensionless photon energy.
        def spectral(T, photon):
            lo = HP * C / (912e-8 * KB * T) if photon else HP * C / (4000e-8 * KB * T)
            hi = 700 if photon else min(700, HP * C / (100e-8 * KB * T))
            if lo >= hi:
                return 0.0
            power = 2 if photon else 3
            return quad(
                lambda x: x**power * np.exp(-x) / (-np.expm1(-x)), lo, hi, epsabs=1e-10
            )[0]

        # Tabulate spectra on a smaller smooth mass grid, then interpolate.
        sm = np.geomspace(0.08, 100, 512)
        st = 5772 * sm**0.625
        uv = np.array([spectral(T, False) for T in st])
        q = np.array([spectral(T, True) for T in st])
        luv = (
            8
            * np.pi**2
            * radius**2
            * (KB * temp) ** 4
            / (HP**3 * C**2)
            * np.interp(m, sm, uv)
        )
        qion = (
            8
            * np.pi**2
            * radius**2
            * (KB * temp) ** 3
            / (HP**3 * C**2)
            * np.interp(m, sm, q)
        )
        lya = 0.9 * (2 / 3) * HP * C / (1215.67e-8) * qion
        self.integrals = []
        for light in (luv, lya):
            self.integrals.append(
                (
                    cumulative_trapezoid(phi * light, m, initial=0),
                    cumulative_trapezoid(phi * light * lifetime(m), m, initial=0),
                )
            )

    def cumulative(self, age, channel):
        """Integral of age kernel from 0 to age: SN/Msun or erg/s/Msun*Myr."""
        age = np.asarray(age, dtype=float)
        positive = np.maximum(age, 0)
        m = np.full_like(positive, 100.0)
        use = positive > 0
        m[use] = np.clip((positive[use] / 1e4) ** -0.4, 0.08, 100)
        if channel == 0:
            lower = np.maximum(m, 8)
            result = self.count[-1] - np.interp(lower, self.mass, self.count)
        else:
            live, dead = self.integrals[channel - 1]
            result = (
                positive * np.interp(m, self.mass, live)
                + dead[-1]
                - np.interp(m, self.mass, dead)
            )
        return np.where(age > 0, result, 0.0)

    def instantaneous(self, age, channel):
        age = np.asarray(age, dtype=float)
        if channel == 0:
            return sn_kernel(age, self.imf)
        m = np.full_like(age, 100.0)
        use = age > 0
        m[use] = np.clip((age[use] / 1e4) ** -0.4, 0.08, 100)
        return np.interp(m, self.mass, self.integrals[channel - 1][0])


@lru_cache(maxsize=32)
def kernels(alpha, mchar, beta, radius_power):
    """Persist precomputed kernels, bound to implementation and physical parameters."""
    import hashlib
    import json

    key = hashlib.sha256(
        json.dumps([alpha, mchar, beta, radius_power, code_fingerprint()]).encode()
    ).hexdigest()
    path = ROOT / "data/cache/kernels" / f"{key}.npz"
    if path.exists():
        with np.load(path, allow_pickle=False) as d:
            k = Kernels.__new__(Kernels)
            k.imf = IMF(alpha, mchar, beta)
            k.mass = d["mass"]
            k.count = d["count"]
            k.integrals = [(d["uv0"], d["uv1"]), (d["lya0"], d["lya1"])]
        return k
    k = Kernels(alpha, mchar, beta, radius_power)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        path,
        mass=k.mass,
        count=k.count,
        uv0=k.integrals[0][0],
        uv1=k.integrals[0][1],
        lya0=k.integrals[1][0],
        lya1=k.integrals[1][1],
    )
    return k


class History:
    """Causal exponential SFH convolution, with a continuous optional cutoff.

    The exact reservoir solution is a sum of two exponentials. Preintegrate
    exp(-rate*age)*kernel(age) on the independently controlled age grid.
    No future birth mass or mutable solver history is used.
    """

    def __init__(self, cfg, stop=None):
        from .physics import growth_rate, FB, hubble, MYR

        self.cfg = cfg
        self.stop = stop
        dt = cfg.history_dt
        self.times = np.arange(int(np.ceil(cfg.tmax / dt)) + 1) * dt
        self.k = kernels(cfg.alpha, cfg.mchar, cfg.beta, cfg.stellar_radius_power)
        a = growth_rate(cfg)
        b = cfg.fstar / (0.0026 / hubble(cfg.z) / MYR)
        q = FB * cfg.mass * a / (a + b) if a + b else 0.0
        self.rates = np.array([a, -b])
        self.amplitudes = np.array([b * q, b * (cfg.gas_fraction * FB * cfg.mass - q)])
        mid = (self.times[:-1] + self.times[1:]) / 2
        self.weighted = []
        for rate in self.rates:
            channels = []
            for ch in range(3):
                weights = np.diff(self.k.cumulative(self.times, ch)) * np.exp(
                    -rate * mid
                )
                channels.append(np.r_[0, np.cumsum(weights)])
            self.weighted.append(channels)
        self.values = self.evaluate(self.times)

    def evaluate(self, t):
        t = np.asarray(t, dtype=float)
        lower = np.zeros_like(t) if self.stop is None else np.maximum(t - self.stop, 0)
        out = np.zeros((3,) + t.shape)
        for amplitude, rate, channels in zip(
            self.amplitudes, self.rates, self.weighted
        ):
            for ch, arr in enumerate(channels):
                out[ch] += (
                    amplitude
                    * np.exp(rate * t)
                    * (
                        np.interp(t, self.times, arr)
                        - np.interp(lower, self.times, arr)
                    )
                )
        # Exactly known finite SN support after a cutoff, not a numerical floor.
        if self.stop is not None:
            out[0] = np.where(lower >= lifetime(8), 0.0, out[0])
        return out

    def __call__(self, t):
        if t < 0 or t > self.times[-1]:
            raise ValueError("history time outside table")
        return self.evaluate(t)
