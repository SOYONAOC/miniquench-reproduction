"""Constrained formation -> launch -> outward/fallback -> moving-boundary return."""

from dataclasses import asdict
import time
import json
import numpy as np
from scipy.integrate import solve_ivp
from .physics import (
    Config,
    ROOT,
    G,
    MSUN,
    MYR,
    FB,
    halo,
    boundary,
    reservoir,
    growth_rate,
    pressure_rhs,
    cache_key,
    code_fingerprint,
)
from .feedback import History
from .radiation import radiation_forces


def shell_budget(t, r, v, p, mseg, cfg, lum, loading):
    mh, rv, vin = halo(t, cfg)
    mdot = FB * growth_rate(cfg) * mh if loading else 0.0  # Msun/Myr
    grav = -cfg.gravity_factor * G * mh * MSUN * mseg * MSUN / (rv * r)
    acc = -mdot * MSUN / MYR * (vin + v)
    uv, lya, temp, mult = radiation_forces(r, v, mseg, lum, cfg)
    pressure = 4 * np.pi * r * r * p
    if cfg.model == "no_feedback":
        uv = lya = pressure = 0.0
    forces = np.array(
        [grav, acc, pressure * cfg.cover, uv * cfg.cover, lya * cfg.cover]
    )
    return forces, mdot, temp, mult


def run(cfg: Config):
    start = time.perf_counter()
    hist = History(cfg)
    rscale = boundary(0, cfg)[0]
    mscale = FB * cfg.mass
    pscale = 1e51 * (cfg.mass / 1e8) / (2 * np.pi * rscale**3)
    vscale = rscale / MYR
    nfev = 0
    rhs_calls = 0

    def luminosity(t, history):
        lum = history(t)
        if cfg.model == "no_feedback":
            lum[:] = 0
        return lum

    def formation_rhs(t, y):
        r, v, _ = boundary(t, cfg)
        lum = luminosity(t, hist)
        dp = pressure_rhs(
            y[0] * pscale,
            r,
            v if cfg.boundary_work else 0,
            cfg.epsilon * lum[0] * 1e51 / MYR,
            cfg.z,
            cfg.cooling,
        )
        return [dp * MYR / pscale]

    def launch(t, y):
        r, v, a = boundary(t, cfg)
        gas = reservoir(t, cfg)[0]
        forces, _, _, _ = shell_budget(
            t, r, v, y[0] * pscale, cfg.cover * gas, cfg, luminosity(t, hist), True
        )
        return (np.sum(forces) - cfg.cover * gas * MSUN * a) / (
            mscale * MSUN * rscale / MYR**2
        )

    launch.terminal = True
    launch.direction = 1
    solA = solve_ivp(
        formation_rhs,
        (0, cfg.sf_window),
        [0.0],
        events=launch,
        rtol=cfg.rtol,
        atol=cfg.atol,
        max_step=cfg.max_step,
        dense_output=True,
    )
    nfev += solA.nfev
    if not solA.success:
        raise ArithmeticError(solA.message)
    tsf = float(solA.t[-1]) if len(solA.t_events[0]) else None
    pieces = []
    turn = None
    returned = None
    relative_turn = None
    status = "not_launched_within_window"
    if tsf is not None:
        histB = History(cfg, stop=tsf)
        gas, stars, _ = reservoir(tsf, cfg)
        r, v, _ = boundary(tsf, cfg)
        y0 = [
            r / rscale,
            v / vscale,
            float(solA.y[0, -1]),
            cfg.cover * gas / mscale,
            0.0,
        ]

        def rhs(t, y):
            nonlocal rhs_calls
            rhs_calls += 1
            if rhs_calls > 50000:
                raise ArithmeticError(
                    "50000 RHS evaluations exceeded; investigate discontinuous loading / unresolved motion"
                )
            r, u, p, mseg, central = y
            if r <= 0 or mseg <= 0:
                raise ArithmeticError("nonpositive physical state")
            v = u * vscale
            lum = luminosity(t, histB)
            loading = v > 0 or cfg.loading == "until_return"
            forces, mdot, _, _ = shell_budget(
                t, r * rscale, v, p * pscale, mseg * mscale, cfg, lum, loading
            )
            dp = pressure_rhs(
                p * pscale,
                r * rscale,
                v,
                cfg.epsilon * lum[0] * 1e51 / MYR,
                cfg.z,
                cfg.cooling,
            )
            return [
                u,
                np.sum(forces) / (mseg * mscale * MSUN) * MYR / vscale,
                dp * MYR / pscale,
                mdot / mscale,
                (FB * growth_rate(cfg) * halo(t, cfg)[0] - mdot) / mscale,
            ]

        def ratio_turn(t, y):
            # d(r/R0)/dt=0, arm return only after an actual relative maximum.
            if t <= tsf + 1e-6:
                return 1e-13
            r0, v0, _ = boundary(t, cfg)
            return y[1] * vscale - y[0] * rscale * v0 / r0

        ratio_turn.terminal = True
        ratio_turn.direction = -1

        def physical_turn(t, y):
            return y[1]

        physical_turn.direction = -1
        physical_turn.terminal = False

        def return_event(t, y):
            return y[0] - boundary(t, cfg)[0] / rscale

        return_event.direction = -1
        return_event.terminal = True
        s1 = solve_ivp(
            rhs,
            (tsf, cfg.tmax),
            y0,
            events=(ratio_turn, physical_turn),
            rtol=cfg.rtol,
            atol=cfg.atol,
            max_step=cfg.max_step,
            dense_output=True,
        )
        pieces.append(s1)
        nfev += s1.nfev
        if not s1.success:
            raise ArithmeticError(s1.message)
        if len(s1.t_events[1]):
            turn = float(s1.t_events[1][0])
        status = "launched_not_returned_by_tmax"
        if len(s1.t_events[0]):
            relative_turn = float(s1.t_events[0][0])
            if s1.y[0, -1] <= boundary(relative_turn, cfg)[0] / rscale:
                raise ArithmeticError("launch did not separate from moving boundary")
            s2 = solve_ivp(
                rhs,
                (s1.t[-1], cfg.tmax),
                s1.y[:, -1],
                events=(return_event, physical_turn),
                rtol=cfg.rtol,
                atol=cfg.atol,
                max_step=cfg.max_step,
                dense_output=True,
            )
            pieces.append(s2)
            nfev += s2.nfev
            if not s2.success:
                raise ArithmeticError(s2.message)
            if len(s2.t_events[1]):
                turn = float(s2.t_events[1][0])
            if len(s2.t_events[0]):
                returned = float(s2.t_events[0][0])
                status = "completed_cycle"
    end = solA.t[-1] if not pieces else pieces[-1].t[-1]
    events = [x for x in (tsf, turn, relative_turn, returned, end) if x is not None]
    times = np.unique(np.r_[np.arange(0, end, cfg.output_dt), events])
    rows = []
    for t in times:
        if tsf is None or t < tsf:
            r, v, _ = boundary(t, cfg)
            p = float(solA.sol(t)[0]) * pscale
            gas, stars, sfr = reservoir(t, cfg)
            mseg = cfg.cover * gas
            center = 0.0
            other = (1 - cfg.cover) * gas
            lum = luminosity(t, hist)
            loading = True
        else:
            s = pieces[0] if len(pieces) == 1 or t <= pieces[0].t[-1] else pieces[1]
            y = s.sol(t)
            r, v, p, mseg, center = y * np.array(
                [rscale, vscale, pscale, mscale, mscale]
            )
            gas0, stars, _ = reservoir(tsf, cfg)
            other = (1 - cfg.cover) * gas0
            gas = mseg + other + center
            sfr = 0.0
            lum = luminosity(t, histB)
            loading = v > 0 or cfg.loading == "until_return"
        forces, mdot, temp, mult = shell_budget(t, r, v, p, mseg, cfg, lum, loading)
        rows.append(
            [
                t,
                r / boundary(t, cfg)[0],
                r,
                v / 1e5,
                p,
                sfr / 1e6,
                gas,
                stars,
                mseg,
                center,
                other,
                halo(t, cfg)[0],
                *forces,
                lum[0],
                lum[1],
                lum[2],
                temp,
                mult,
                mdot,
            ]
        )
    names = [
        "time_myr",
        "r_over_R0",
        "radius_cm",
        "velocity_km_s",
        "pressure_erg_cm3",
        "sfr_msun_yr",
        "gas_msun",
        "formed_msun",
        "incident_shell_msun",
        "central_gas_msun",
        "other_shell_msun",
        "halo_msun",
        "gravity_dyne",
        "accretion_dyne",
        "pressure_dyne",
        "uv_dyne",
        "lya_dyne",
        "sn_per_myr",
        "uv_erg_s",
        "lya_erg_s",
        "shell_temperature_K",
        "force_multiplier",
        "shell_loading_msun_myr",
    ]
    curve = np.array(rows)
    maxratio = float(np.max(curve[:, 1]))
    if relative_turn is not None:
        maxratio = float(pieces[0].y[0, -1] * rscale / boundary(relative_turn, cfg)[0])
    if tsf is not None and maxratio - 1 <= 10 * (cfg.rtol + cfg.atol):
        raise ArithmeticError(
            "Grazing launch: displacement is below integration resolution; no resolved quenching cycle"
        )
    result = dict(
        status=status,
        t_sf_myr=tsf,
        t_turn_myr=turn,
        t_relative_max_myr=relative_turn,
        t_return_myr=returned,
        t_quench_myr=None if returned is None else returned - tsf,
        t_quench_lower_bound_myr=float(end - tsf)
        if tsf is not None and returned is None
        else None,
        max_r_over_R0=maxratio,
        wall_seconds=time.perf_counter() - start,
        nfev=nfev,
        config=asdict(cfg),
        code_fingerprint=code_fingerprint(),
    )
    return result, dict(zip(names, curve.T))


def run_cached(cfg):
    key = cache_key(cfg)
    folder = ROOT / "data/cache/trajectories"
    folder.mkdir(parents=True, exist_ok=True)
    meta = folder / f"{key}.json"
    data = folder / f"{key}.npz"
    if meta.exists() and data.exists():
        result = json.loads(meta.read_text())
        if result["status"] == "numerical_failure":
            raise ArithmeticError(result["error"])
        return result, dict(np.load(data, allow_pickle=False)), True
    started = time.perf_counter()
    try:
        result, curve = run(cfg)
    except (ValueError, ArithmeticError, FloatingPointError) as error:
        meta.write_text(
            json.dumps(
                {
                    "status": "numerical_failure",
                    "error": str(error),
                    "config": asdict(cfg),
                    "wall_seconds": time.perf_counter() - started,
                    "code_fingerprint": code_fingerprint(),
                },
                indent=2,
            )
        )
        raise
    np.savez_compressed(data, **curve)
    meta.write_text(json.dumps(result, indent=2, allow_nan=False))
    return result, curve, False
