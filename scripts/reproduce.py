"""Reproducible CPU scans; immutable configuration + code-bound checkpoint cache."""

import os

for key in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[key] = "1"
import argparse
import json
import time
from dataclasses import asdict, replace
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from miniquench.physics import Config, ROOT, FB, code_fingerprint
from miniquench.dynamics import run_cached
from miniquench.analytic import analytic_mass, numerical_mass
from miniquench.feedback import IMF

plt.style.use(ROOT / "data/apj.mplstyle")
plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "DejaVu Serif",
        "font.serif": ["DejaVu Serif"],
        "font.size": 10,
        "figure.dpi": 140,
    }
)


def worker(payload):
    label, mass, config = payload
    cfg = Config(**{**config, "mass": float(mass)})
    started = time.perf_counter()
    try:
        result, curve, cached = run_cached(cfg)
    except ArithmeticError as error:
        result = dict(
            status="numerical_failure",
            error=str(error),
            t_sf_myr=None,
            t_quench_myr=None,
            t_quench_lower_bound_myr=None,
            max_r_over_R0=None,
            wall_seconds=time.perf_counter() - started,
            config=asdict(cfg),
        )
        curve = None
        cached = False
    return label, result, curve, cached


MODEL_TAG = "full_approx"


def savefig(fig, name):
    if name.startswith(("figure2", "figure5", "figure8")):
        name += "_" + MODEL_TAG
    fig.savefig(ROOT / "outputs/figures" / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(ROOT / "outputs/figures" / f"{name}.png", bbox_inches="tight", dpi=170)
    plt.close(fig)


def reference(ax, figure, label, color=None):
    path = ROOT / "data/reference" / f"fig{figure}_{label}.csv"
    a = np.genfromtxt(path, delimiter=",", names=True)
    ax.plot(
        a["x"],
        a["y"],
        ":",
        color=color,
        lw=1.3,
        label="paper vector" if label == "m8" else None,
    )


def imf_plot():
    x = np.geomspace(0.08, 100, 1000)
    fig, ax = plt.subplots(figsize=(7, 4))
    columns = [x]
    names = ["mass_msun"]
    for alpha in (1.35, 2.35):
        for mc, color in ((0.2, "tab:purple"), (5, "tab:blue"), (10, "tab:green")):
            y = IMF(alpha, mc)(x)
            ax.loglog(
                x,
                y,
                "--" if alpha == 2.35 else "-",
                color=color,
                label=f"alpha={alpha}, mc={mc}",
            )
            columns.append(y)
            names.append(f"alpha{alpha}_mc{mc}")
            reference(ax, 1, f"alpha{alpha}_mc{mc:g}", color)
    ax.loglog(
        x, IMF(-2.35, 0.2)(x), ":", color="gray", label="literal Fig.5 alpha=-2.35"
    )
    ax.axvline(8, color="gray", lw=0.7)
    ax.set(
        xlabel=r"Stellar mass [$M_\odot$]",
        ylabel=r"$dN/dm$ per formed $M_\odot$",
        ylim=(1e-5, 10),
        xlim=(0.08, 100),
        title="Fig. 1: Eq. (7) and explicit caption-conflict check",
    )
    ax.legend(fontsize=8, ncol=2)
    np.savetxt(
        ROOT / "data/curves/figure1.csv",
        np.array(columns).T,
        delimiter=",",
        header=",".join(names),
        comments="",
    )
    savefig(fig, "figure1_imf")


def fig4(dt):
    t = np.linspace(10, 100, 181)
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.7), sharey=True)
    cols = [t, analytic_mass(t)]
    names = ["time_myr", "analytic_z6_msun"]
    for ax, z in zip(axes, [3, 6, 9]):
        for mc, color, label in [
            (0.2, "tab:red", "chabrier"),
            (5, "tab:blue", "topheavy"),
        ]:
            for growth, ls in [(True, "-"), (False, "--")]:
                m = numerical_mass(t, z, mc, growth, dt)
                ax.semilogy(
                    t,
                    m,
                    ls,
                    color=color,
                    label=f"{label}, " + ("growth" if growth else "fixed halo"),
                )
                cols.append(m)
                names.append(f"z{z}_{label}_growth{int(growth)}")
            reference(ax, 4, f"z{z}_{label}", color)
        if z == 6:
            ax.semilogy(t, analytic_mass(t), "-", color="k", lw=1, label="Eq. (33)")
        ax.set(
            title=f"z={z}", xlabel=r"$t_{SF}$ [Myr]", xlim=(10, 100), ylim=(5e9, 3e13)
        )
    axes[0].set_ylabel(r"Maximum initial halo mass [$M_\odot$]")
    axes[0].legend(fontsize=7)
    fig.suptitle("Fig. 4: two-force model; dotted = paper; no radiative closure")
    fig.tight_layout()
    np.savetxt(
        ROOT / "data/curves/figure4.csv",
        np.array(cols).T,
        delimiter=",",
        header=",".join(names),
        comments="",
    )
    savefig(fig, "figure4_mass_limit")


def families(figures, base):
    groups = {}
    if 2 in figures:
        groups["fig2"] = base
    if 5 in figures:
        for f, label in [(0.001, "f001"), (0.01, "f01"), (0.1, "f1")]:
            groups["fig5_" + label] = replace(base, fstar=f)
        for mc, label in [(1, "mc1"), (5, "mc5")]:
            groups["fig5_" + label] = replace(base, fstar=0.01, mchar=mc)
    if 8 in figures:
        for f, label in [(0.001, "f001"), (0.01, "f01")]:
            for cover, cl in [
                (1, "c1"),
                (0.1, "c01"),
                (0.05, "c005"),
                (0.005, "c0005"),
            ]:
                groups[f"fig8_{label}_{cl}"] = replace(base, fstar=f, cover=cover)
            groups[f"fig8_{label}_lowgas"] = replace(base, fstar=f, gas_fraction=0.1)
    return groups


def getxy(rows):
    rows = sorted(rows, key=lambda r: r["config"]["mass"])
    x = np.array([r["config"]["mass"] for r in rows])
    y = np.array(
        [
            r["t_quench_myr"] if r["status"] == "completed_cycle" else np.nan
            for r in rows
        ]
    )
    return x, y, rows


def plot_scan(ax, rows, label, color, ls="-"):
    x, y, rs = getxy(rows)
    ax.semilogx(x, y, ls, color=color, label=label)
    # Censoring is a lower bound, never an assigned zero-duration quench.
    for xi, r in zip(x, rs):
        if r["status"] == "launched_not_returned_by_tmax":
            ax.plot(xi, r["t_quench_lower_bound_myr"], "^", color=color, ms=3)
        elif r["status"] == "not_launched_within_window":
            ax.plot(xi, 0, "x", color=color, ms=3)
        elif r["status"] == "numerical_failure":
            ax.plot(xi, 0, "|", color="red", ms=7)


def scan_upper(results):
    values = [200.0]
    for rows in results.values():
        for r in rows:
            for k in ("t_quench_myr", "t_quench_lower_bound_myr"):
                if r.get(k) is not None:
                    values.append(r[k])
    return max(values) * 1.08


def plot_results(results, curves, figures, model):
    if 2 in figures:
        fig, axes = plt.subplots(2, 3, figsize=(12, 7))
        for mass, color in [(1e8, "tab:blue"), (1e9, "tab:red")]:
            if ("fig2", mass) not in curves:
                axes[0, 0].text(
                    0.05,
                    0.8,
                    f"Mh={mass:.0e}: numerical_failure",
                    transform=axes[0, 0].transAxes,
                    color=color,
                )
                reference(axes[0, 0], 2, "m8" if mass == 1e8 else "m9", color)
                continue
            c = curves[("fig2", mass)]
            r = next(r for r in results["fig2"] if r["config"]["mass"] == mass)
            t = c["time_myr"]
            lab = (
                f"Mh={mass:.0e}, quench={r['t_quench_myr']:.1f} Myr"
                if r["t_quench_myr"] is not None
                else r["status"]
            )
            for ax, name in zip(
                axes.flat,
                [
                    "r_over_R0",
                    "velocity_km_s",
                    "pressure_erg_cm3",
                    "sfr_msun_yr",
                    "gas_msun",
                    "formed_msun",
                ],
            ):
                ax.plot(t, c[name], color=color, label=lab)
                ax.set(
                    xlabel="Time [Myr]",
                    ylabel={
                        "r_over_R0": r"$r/R_0$",
                        "velocity_km_s": "Velocity [km/s]",
                        "pressure_erg_cm3": r"Pressure [erg cm$^{-3}$]",
                        "sfr_msun_yr": r"SFR [$M_\odot$/yr]",
                        "gas_msun": r"Gas mass [$M_\odot$]",
                        "formed_msun": r"Formed mass [$M_\odot$]",
                    }[name],
                )
            reference(axes[0, 0], 2, "m8" if mass == 1e8 else "m9", color)
            for ax in axes.flat:
                if r["t_sf_myr"] is not None:
                    ax.axvline(r["t_sf_myr"], color=color, ls="--", lw=0.6)
        axes[0, 0].legend(fontsize=7)
        fig.suptitle(f"Fig. 2: {model}; dotted trajectories = paper")
        fig.tight_layout()
        savefig(fig, "figure2_trajectories")
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        for ax, mass in zip(axes, [1e8, 1e9]):
            if ("fig2", mass) not in curves:
                ax.text(
                    0.1,
                    0.5,
                    "numerical_failure: no valid trajectory",
                    transform=ax.transAxes,
                )
                continue
            c = curves[("fig2", mass)]
            for name in [
                "gravity_dyne",
                "accretion_dyne",
                "pressure_dyne",
                "uv_dyne",
                "lya_dyne",
            ]:
                ax.plot(c["time_myr"], c[name], label=name)
            ax.set_yscale("symlog", linthresh=1e30)
            ax.set(
                xlabel="Time [Myr]",
                ylabel="Signed force [dyne]",
                title=f"Mh={mass:.0e}",
            )
        axes[0].legend(fontsize=7)
        fig.tight_layout()
        savefig(fig, "figure2_forces")
    if 5 in figures:
        fig, axes = plt.subplots(2, 2, figsize=(10, 8))
        sets = [
            [("f001", "tab:blue"), ("f01", "tab:purple"), ("f1", "tab:red")],
            [("f01", "tab:purple"), ("mc1", "darkorange"), ("mc5", "seagreen")],
        ]
        for row, items in enumerate(sets):
            for label, color in items:
                rows = results["fig5_" + label]
                plot_scan(axes[row, 0], rows, label, color)
                reference(
                    axes[row, 0],
                    5,
                    label if row == 0 or label != "f01" else "mc02",
                    color,
                )
                if label == "f1":
                    continue
                x, y, _ = getxy(rows)
                axes[row, 1].semilogx(x * FB * 0.01, y, color=color)
                axes[row, 1].semilogx(x * FB * 0.1, y, color=color, label=label)
                good = np.isfinite(y)
                axes[row, 1].fill(
                    np.r_[x[good] * FB * 0.01, (x[good] * FB * 0.1)[::-1]],
                    np.r_[y[good], y[good][::-1]],
                    color=color,
                    alpha=0.09,
                )
            for ax in axes[row]:
                ax.set(ylabel="Quenching duration [Myr]", ylim=(0, scan_upper(results)))
                ax.legend(fontsize=8)
            axes[row, 0].set(xlabel=r"Initial halo mass [$M_\odot$]", xlim=(1e8, 1e13))
            axes[row, 1].set(
                xlabel=r"Mapped total stellar mass [$M_\odot$]", xlim=(1e6, 1e11)
            )
        fig.suptitle(
            f"Fig. 5: {model}; dotted: paper; x: no launch; red tick: failure; triangle: lower bound"
        )
        fig.tight_layout()
        savefig(fig, "figure5_efficiency_imf")
    if 8 in figures:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
        for ax, f in zip(axes, ["f001", "f01"]):
            for label, color, ls in [
                ("c1", "tab:blue", "-"),
                ("c01", "tab:purple", "-"),
                ("c005", "tab:orange", "-"),
                ("c0005", "tab:green", "-"),
                ("lowgas", "tab:blue", "--"),
            ]:
                plot_scan(ax, results[f"fig8_{f}_{label}"], label, color, ls)
                reference(ax, 8, f"{f}_{label}", color)
            ax.set(
                xlabel=r"Initial halo mass [$M_\odot$]",
                ylabel="Quenching duration [Myr]",
                title=r"$f_\star=$" + ("0.001" if f == "f001" else "0.01"),
                ylim=(0, scan_upper(results)),
                xlim=(1e8, 1e13),
            )
        axes[0].legend(fontsize=8)
        fig.suptitle(
            f"Fig. 8: {model}; dotted: paper; x: no launch; red tick: failure; triangle: lower bound"
        )
        fig.tight_layout()
        savefig(fig, "figure8_cover_gas")


def main():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--figures",
        nargs="+",
        type=int,
        default=[1, 2, 4, 5, 8],
        choices=[1, 2, 4, 5, 8],
    )
    p.add_argument("--profile", choices=["quick", "full"], default="quick")
    p.add_argument("--model", choices=["full_approx", "sn_only"], default="full_approx")
    p.add_argument(
        "--loading", choices=["outward_only", "until_return"], default="until_return"
    )
    p.add_argument("--workers", type=int, default=4)
    args = p.parse_args()
    global MODEL_TAG
    MODEL_TAG = args.model if args.profile == "full" else args.model + "_quick"
    if not 1 <= args.workers <= 4:
        p.error("workers must be 1..4")
    start = time.perf_counter()
    base = Config(
        model=args.model,
        loading=args.loading,
        history_dt=0.2 if args.profile == "quick" else 0.05,
        max_step=1 if args.profile == "quick" else 0.5,
        output_dt=1 if args.profile == "quick" else 0.5,
    )
    groups = families(args.figures, base)
    masses = np.logspace(8, 13, 9 if args.profile == "quick" else 31)
    for figure in args.figures:
        config_data = {
            "figure": figure,
            "profile": args.profile,
            "baseline": asdict(base),
            "mass_grid": ([1e8, 1e9] if figure == 2 else masses.tolist()),
            "refinement_rounds": 4
            if args.profile == "full" and figure in (5, 8)
            else 0,
            "families": {
                k: asdict(v) for k, v in groups.items() if k.startswith(f"fig{figure}")
            },
        }
        if figure in (1, 4):
            config_data.pop("baseline")
            config_data.pop("mass_grid")
            config_data["physical_inputs"] = (
                {
                    "alpha": [1.35, 2.35],
                    "mchar": [0.2, 5, 10],
                    "beta": 1.6,
                    "mass_min": 0.08,
                    "mass_max": 100,
                    "literal_caption_alpha": -2.35,
                }
                if figure == 1
                else {
                    "fstar": 0.01,
                    "z": [3, 6, 9],
                    "mchar": [0.2, 5],
                    "alpha": 2.35,
                    "beta": 1.6,
                    "time_min_myr": 10,
                    "time_max_myr": 100,
                    "time_points": 181,
                    "forces": ["SN", "gravity"],
                    "growth": [True, False],
                    "cooling": False,
                    "radiation": False,
                    "history_dt_myr": base.history_dt,
                }
            )
        (
            ROOT / "configs" / f"figure{figure}_{args.profile}_{args.model}.json"
        ).write_text(json.dumps(config_data, indent=2))
    if 1 in args.figures:
        imf_plot()
    if 4 in args.figures:
        fig4(base.history_dt)
    payloads = [
        (label, mass, asdict(cfg))
        for label, cfg in groups.items()
        for mass in ([1e8, 1e9] if label == "fig2" else masses)
    ]
    results = {label: [] for label in groups}
    curves = {}
    hits = 0

    def consume(items):
        nonlocal hits
        for label, result, curve, cached in items:
            results[label].append(result)
            hits += int(cached)
            mass = result["config"]["mass"]
            if label == "fig2" and curve is not None:
                curves[label, mass] = curve
                path = (
                    ROOT
                    / "data/curves"
                    / f"{label}_{mass:.0e}_{args.profile}_{args.model}.csv"
                )
                np.savetxt(
                    path,
                    np.array(list(curve.values())).T,
                    delimiter=",",
                    header=",".join(curve),
                    comments="",
                )
            print(
                f"{label} M={mass:.3g}: {result['status']}; SF={result['t_sf_myr']}; quench={result['t_quench_myr']}; {result['wall_seconds']:.3f}s; cached={cached}",
                flush=True,
            )

    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        consume(pool.map(worker, payloads))
        # Refine launch-status boundaries, without changing physical parameters.
        if args.profile == "full":
            for level in range(4):
                extra = []
                for label, rs in results.items():
                    if label == "fig2":
                        continue
                    _, _, rs = getxy(rs)
                    for lo, hi in zip(rs[:-1], rs[1:]):
                        if (
                            lo["status"] == "numerical_failure"
                            or hi["status"] == "numerical_failure"
                        ):
                            continue
                        if (lo["t_sf_myr"] is None) != (hi["t_sf_myr"] is None):
                            m = np.sqrt(lo["config"]["mass"] * hi["config"]["mass"])
                            extra.append((label, m, asdict(groups[label])))
                if extra:
                    consume(pool.map(worker, extra))
    plot_results(results, curves, args.figures, args.model)
    report = {
        "profile": args.profile,
        "model": args.model,
        "loading": args.loading,
        "wall_seconds": time.perf_counter() - start,
        "workers": args.workers,
        "cache_hits": hits,
        "code_fingerprint": code_fingerprint(),
        "results": results,
    }
    out = (
        ROOT
        / "data/curves"
        / (
            f"summary_{args.profile}_{args.model}_"
            + "".join(map(str, args.figures))
            + ".json"
        )
    )
    out.write_text(json.dumps(report, indent=2, allow_nan=False))
    # Flat machine-readable summary, including null durations as blank cells.
    import csv

    with (ROOT / "data/curves" / f"scan_{args.profile}_{args.model}.csv").open(
        "w"
    ) as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "family",
                "mass_msun",
                "status",
                "t_sf_myr",
                "t_quench_myr",
                "quench_lower_bound_myr",
                "max_r_over_R0",
                "wall_seconds",
            ]
        )
        for label, rs in results.items():
            for r in sorted(rs, key=lambda r: r["config"]["mass"]):
                writer.writerow(
                    [
                        label,
                        r["config"]["mass"],
                        r["status"],
                        r["t_sf_myr"],
                        r["t_quench_myr"],
                        r["t_quench_lower_bound_myr"],
                        r["max_r_over_R0"],
                        r["wall_seconds"],
                    ]
                )
    print(
        f"Total elapsed {report['wall_seconds']:.3f}s; cache hits {hits}; {out}",
        flush=True,
    )


if __name__ == "__main__":
    main()
