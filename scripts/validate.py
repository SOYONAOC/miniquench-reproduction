"""Independent numerical convergence, finite convention branches and paper errors."""

import os

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
import json
import time
import csv
from dataclasses import replace
import numpy as np
from miniquench.physics import ROOT, Config, code_fingerprint
from miniquench.dynamics import run_cached


def main():
    start = time.perf_counter()
    rows = []
    summaries = []
    for mass in (1e8, 1e9):
        base = Config(
            mass=mass,
            model="full_approx" if mass == 1e8 else "sn_only",
            loading="until_return",
            history_dt=0.1,
            max_step=0.5,
            rtol=2e-7,
            atol=1e-9,
        )
        specs = [
            ("baseline", {}),
            ("history_half", {"history_dt": 0.05}),
            ("history_quarter", {"history_dt": 0.025}),
            ("tolerance_tight", {"rtol": 2e-9, "atol": 1e-11}),
            ("step_half", {"max_step": 0.25}),
            (
                "combined_fine",
                {"history_dt": 0.025, "max_step": 0.25, "rtol": 2e-9, "atol": 1e-11},
            ),
        ]
        group = []
        for label, change in specs:
            r, _, hit = run_cached(replace(base, **change))
            group.append((label, r))
            print(
                "convergence",
                mass,
                label,
                r["status"],
                r["t_sf_myr"],
                r["t_quench_myr"],
                flush=True,
            )
        ref = group[-1][1]
        for label, r in group:
            row = {
                "mass_msun": mass,
                "model": r["config"]["model"],
                "case": label,
                "status": r["status"],
                "wall_seconds": r["wall_seconds"],
            }
            for key in ("t_sf_myr", "t_quench_myr", "max_r_over_R0"):
                row[key] = r[key]
                row[key + "_relative_error"] = (
                    None
                    if r[key] is None or ref[key] is None
                    else abs(r[key] - ref[key]) / abs(ref[key])
                )
                row[key + "_absolute_error"] = (
                    None
                    if r[key] is None or ref[key] is None
                    else abs(r[key] - ref[key])
                )
            rows.append(row)
        tested = [r for r in rows if r["mass_msun"] == mass]
        passes = all(
            r[k + "_relative_error"] is not None and r[k + "_relative_error"] <= 0.01
            for r in tested
            for k in ("t_sf_myr", "t_quench_myr", "max_r_over_R0")
        )
        summaries.append(
            {
                "mass_msun": mass,
                "model": base.model,
                "all_variations_within_1_percent": passes,
            }
        )
    with (ROOT / "data/convergence.csv").open("w") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    branches = []
    for mass in (1e8, 1e9):
        base = Config(mass=mass, loading="until_return", history_dt=0.05)
        for label, change in [
            ("baseline", {}),
            ("sn_only", {"model": "sn_only"}),
            ("outward_only", {"loading": "outward_only"}),
            ("isothermal_factor1", {"gravity_factor": 1}),
            ("no_boundary_PdV", {"boundary_work": False}),
            ("thickness001", {"shell_thickness": 0.01}),
            ("thickness03", {"shell_thickness": 0.3}),
            ("lya_velocity", {"lya_velocity": True}),
            ("radius_power07", {"stellar_radius_power": 0.7}),
        ]:
            try:
                r, _, _ = run_cached(replace(base, **change))
            except ArithmeticError as error:
                r = {
                    "status": "numerical_failure",
                    "error": str(error),
                    "config": {"mass": mass},
                    "t_quench_myr": None,
                }
            branches.append({"branch": label, **r})
            print("branch", mass, label, r["t_quench_myr"], flush=True)
    (ROOT / "data/sensitivity.json").write_text(json.dumps(branches, indent=2))
    # Launch threshold convergence (log mass brackets), independent of quench censoring.
    thresholds = []
    for dt in (0.2, 0.1, 0.05):
        cfg = Config(
            fstar=0.001,
            model="sn_only",
            loading="until_return",
            history_dt=dt,
            tmax=700,
            max_step=0.5,
            output_dt=2,
        )
        lo, hi = 8.0, 13.0
        r, _, _ = run_cached(replace(cfg, mass=10**lo))
        s, _, _ = run_cached(replace(cfg, mass=10**hi))
        if r["t_sf_myr"] is None or s["t_sf_myr"] is not None:
            raise AssertionError("threshold not bracketed")
        for _ in range(12):
            mid = (lo + hi) / 2
            r, _, _ = run_cached(replace(cfg, mass=10**mid))
            if r["t_sf_myr"] is not None:
                lo = mid
            else:
                hi = mid
        thresholds.append(
            {
                "history_dt_myr": dt,
                "log10_mass_low": lo,
                "log10_mass_high": hi,
                "mass_low_msun": 10**lo,
                "mass_high_msun": 10**hi,
            }
        )
    (ROOT / "data/threshold_convergence.json").write_text(
        json.dumps(thresholds, indent=2)
    )
    # Compare trajectories only over their common physical support, independent references.
    comparisons = []
    for mass, label in [(1e8, "m8"), (1e9, "m9")]:
        ref = np.genfromtxt(
            ROOT / "data/reference" / f"fig2_{label}.csv", delimiter=",", names=True
        )
        r, c, _ = run_cached(
            Config(
                mass=mass,
                model="full_approx" if mass == 1e8 else "sn_only",
                loading="until_return",
                history_dt=0.05,
            )
        )
        good = ref["x"] <= c["time_myr"][-1]
        diff = np.interp(ref["x"][good], c["time_myr"], c["r_over_R0"]) - ref["y"][good]
        comparisons.append(
            {
                "figure": 2,
                "mass_msun": mass,
                "model": r["config"]["model"],
                "quench_myr": r["t_quench_myr"],
                "paper_bracket_rough_myr": 175 if mass == 1e8 else 75,
                "radius_ratio_rms_difference": float(np.sqrt(np.mean(diff**2))),
                "paper_max_r_over_R0": float(np.max(ref["y"])),
                "computed_max_r_over_R0": r["max_r_over_R0"],
                "reference": "published vector plot; NOT author solver data",
                "readout_uncertainty": "~2 Myr, ~0.1 r/R0; bracket labels themselves are rounded",
            }
        )
    summary = {
        "wall_seconds": time.perf_counter() - start,
        "code_fingerprint": code_fingerprint(),
        "convergence": summaries,
        "thresholds": thresholds,
        "paper_comparison": comparisons,
    }
    (ROOT / "data/validation.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)
    if not all(s["all_variations_within_1_percent"] for s in summaries):
        raise SystemExit("Numerical 1% target NOT met; see convergence.csv")


if __name__ == "__main__":
    main()
