import numpy as np
import pytest
from dataclasses import replace
from scipy.integrate import quad, solve_ivp
from miniquench.physics import (
    Config,
    MYR,
    MSUN,
    FB,
    boundary,
    halo,
    reservoir,
    pressure_rhs,
    cache_key,
)
from miniquench.feedback import IMF, lifetime, sn_kernel, kernels, History
from miniquench.dynamics import run, shell_budget
from miniquench.radiation import thermal_equilibrium, radiation_forces


@pytest.mark.parametrize(
    "alpha,mc", [(2.35, 0.2), (2.35, 1), (2.35, 5), (1.35, 10), (-2.35, 0.2)]
)
def test_imf_mass(alpha, mc):
    imf = IMF(alpha, mc)
    assert quad(lambda m: m * imf(m), 0.08, 100)[0] == pytest.approx(1, rel=1e-9)


@pytest.mark.parametrize("mc", [0.2, 5])
def test_sn_integrated_count(mc):
    imf = IMF(mchar=mc)
    n = quad(
        lambda t: float(sn_kernel(t, imf)), lifetime(100), lifetime(8), epsabs=1e-12
    )[0]
    assert n == pytest.approx(imf.sn_number(), rel=1e-8)
    k = kernels(2.35, mc, 1.6, 0.5)
    assert float(k.cumulative(100, 0)) == pytest.approx(n, rel=2e-6)
    assert sn_kernel(0.01, imf) == 0 and sn_kernel(60, imf) == 0


@pytest.mark.parametrize("channel", [0, 1, 2])
def test_instant_constant_and_shutdown(channel):
    k = kernels(2.35, 0.2, 1.6, 0.5)
    # Independent time quadrature of the instantaneous kernel versus cumulative.
    val = quad(
        lambda t: float(k.instantaneous(t, channel)),
        0,
        30,
        points=[0.1],
        epsabs=1e-9,
        epsrel=1e-5,
        limit=200,
    )[0]
    assert val == pytest.approx(float(k.cumulative(30, channel)), rel=5e-5)

    # Exact top-hat SFH: feedback remains after shutdown; delayed SN die by 10+tau8.
    def response(t):
        return k.cumulative(t, channel) - k.cumulative(t - 10, channel)

    assert response(20) > 0
    if channel == 0:
        assert response(70) == pytest.approx(0, abs=1e-12)
        assert float(k.cumulative(100, 0)) == pytest.approx(k.imf.sn_number(), rel=2e-6)
    else:
        assert response(70) > 0


def test_history_causality_and_mass():
    cfg = Config(history_dt=0.02, tmax=100, sf_window=20)
    h = History(cfg, stop=10.013)
    full = History(cfg)
    np.testing.assert_allclose(h(5), full(5), rtol=2e-10)
    # Cross-check the convolution at non-grid stop against independent birth-time quadrature.
    k = h.k
    for ch in [0, 1, 2]:
        expected = quad(
            lambda birth: (
                float(reservoir(birth, cfg)[2]) * float(k.instantaneous(20 - birth, ch))
            ),
            0,
            10.013,
            epsrel=2e-5,
        )[0]
        assert h(20)[ch] == pytest.approx(expected, rel=5e-5)
    assert h(70)[0] == pytest.approx(0, abs=1e-8)


def test_adiabatic_pressure():
    r0 = 1e20
    v = 1e6
    p0 = 1e-9

    def rhs(t, y):
        return [pressure_rhs(y[0], r0 + v * t, v, 0, 6, False)]

    sol = solve_ivp(rhs, (0, 3e13), [p0], rtol=1e-10, atol=1e-23)
    invariant = sol.y[0] * (r0 + v * sol.t) ** 5
    np.testing.assert_allclose(invariant / invariant[0], 1, rtol=1e-9)


def test_compton():
    tc = 120 * MYR * (0.7) ** -4
    sol = solve_ivp(
        lambda t, y: [pressure_rhs(y[0], 1e20, 0, 0, 6)],
        (0, tc),
        [1e-9],
        rtol=1e-10,
        atol=1e-23,
    )
    assert sol.y[0, -1] / 1e-9 == pytest.approx(np.exp(-1), rel=1e-9)


def test_reservoir_conservation():
    cfg = Config()
    t = np.linspace(0, 100, 31)
    gas, stars, sfr = reservoir(t, cfg)
    np.testing.assert_allclose(gas + stars, FB * halo(t, cfg)[0], rtol=1e-14)
    assert np.all(sfr > 0)


def test_cover_and_momentum():
    cfg = Config(model="sn_only")
    r = boundary(10, cfg)[0]
    m = FB * cfg.mass
    v = 1e6
    f, mdot, _, _ = shell_budget(10, r, v, 1e-9, m, cfg, np.zeros(3), True)
    # d(Mv)/dt=external gravity+pressure - Mdot*v_in.
    lhs = np.sum(f) + mdot * MSUN / MYR * v
    rhs = f[0] + f[2] - mdot * MSUN / MYR * halo(10, cfg)[2]
    assert lhs == pytest.approx(rhs, rel=1e-13)
    fc, mc, _, _ = shell_budget(
        10, r, v, 1e-9, m * 0.1, replace(cfg, cover=0.1), np.zeros(3), True
    )
    assert mc == mdot
    assert fc[0] == pytest.approx(0.1 * f[0])
    assert fc[2] == pytest.approx(0.1 * f[2])
    assert fc[1] == f[1]
    # At fcover=1 segment equals whole shell, exactly the isotropic equation.
    assert f[2] == pytest.approx(4 * np.pi * r * r * 1e-9)


def test_no_feedback_no_launch():
    r, c = run(Config(model="no_feedback", history_dt=0.5, max_step=1, tmax=100))
    assert r["status"] == "not_launched_within_window"
    assert r["t_sf_myr"] is None and r["t_quench_myr"] is None
    np.testing.assert_allclose(c["r_over_R0"], 1, atol=1e-14)
    assert c["sfr_msun_yr"][-1] > 0  # observation window is not a forced shutdown


def test_completed_return_and_mass():
    cfg = Config(model="sn_only", history_dt=0.1, max_step=1, output_dt=1)
    r, c = run(cfg)
    assert r["status"] == "completed_cycle"
    assert r["t_sf_myr"] < r["t_turn_myr"] < r["t_return_myr"]
    assert r["t_quench_myr"] == pytest.approx(r["t_return_myr"] - r["t_sf_myr"])
    assert c["r_over_R0"][-1] == pytest.approx(1, abs=1e-8)
    np.testing.assert_allclose(
        c["gas_msun"] + c["formed_msun"], FB * c["halo_msun"], rtol=2e-7
    )
    after = c["time_myr"] >= r["t_sf_myr"]
    assert np.all(c["sfr_msun_yr"][after] == 0)
    assert np.any(c["sn_per_myr"][after] > 0)
    falling = c["velocity_km_s"] < 0
    assert np.all(c["shell_loading_msun_myr"][falling] == 0)


def test_right_censoring():
    r, _ = run(
        Config(model="sn_only", tmax=30, sf_window=20, history_dt=0.1, max_step=0.5)
    )
    assert r["status"] == "launched_not_returned_by_tmax"
    assert r["t_quench_myr"] is None
    assert r["t_quench_lower_bound_myr"] == pytest.approx(30 - r["t_sf_myr"])


def test_invalid_config_and_cache():
    with pytest.raises(ValueError):
        Config(cover=0)
    with pytest.raises(ValueError):
        Config(model="full")
    assert cache_key(Config()) != cache_key(Config(cover=0.1))
    assert cache_key(Config()) != cache_key(Config(loading="until_return"))


def test_radiation_closure():
    assert thermal_equilibrium(0, 1, 1e60) == 1e4
    cfg = Config()
    f = radiation_forces(
        boundary(0, cfg)[0], 0, FB * cfg.mass, np.array([0, 1e40, 1e40]), cfg
    )
    assert f[0] > 0 and f[1] > 0 and f[2] >= 1e4 and f[3] > 0
    with pytest.raises(ValueError):
        thermal_equilibrium(1e100, 1, 1)


def test_shutdown_is_continuous_at_non_grid_time():
    cfg = Config(history_dt=0.2)
    stop = 5.6234
    before = History(cfg)
    after = History(cfg, stop=stop)
    np.testing.assert_allclose(before(stop), after(stop), rtol=1e-14)
    assert np.all(after(stop + 0.01) > 0)


def test_isothermal_mass_and_analytic_join():
    from miniquench.physics import enclosed_mass
    from miniquench.analytic import analytic_mass

    cfg = Config()
    mass, rv, _ = halo(0, cfg)
    assert enclosed_mass(rv, 0, cfg) == pytest.approx(2 * mass)
    assert enclosed_mass(rv / 10, 0, cfg) == pytest.approx(0.2 * mass)
    assert analytic_mass(55) == pytest.approx(4.4e11)
    assert analytic_mass(55 + 1e-7) == pytest.approx(analytic_mass(55 - 1e-7), rel=1e-7)


def test_filament_cycle_mass_and_return():
    cfg = Config(
        cover=0.1,
        model="sn_only",
        loading="until_return",
        mass=1e9,
        history_dt=0.1,
        max_step=1,
        output_dt=1,
    )
    r, c = run(cfg)
    assert r["status"] == "completed_cycle"
    np.testing.assert_allclose(
        c["gas_msun"] + c["formed_msun"], FB * c["halo_msun"], rtol=2e-7
    )
    assert c["incident_shell_msun"][-1] > c["incident_shell_msun"][0]


def test_unresolved_radiative_launch_is_not_success():
    cfg = Config(mass=1e9, history_dt=0.05, loading="until_return")
    with pytest.raises(ArithmeticError, match="launch|Grazing"):
        run(cfg)
