"""Finite, source-motivated Fig.4 approximations; no free fitted coefficients."""

import json
import numpy as np
from scipy.integrate import cumulative_trapezoid
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from miniquench.physics import ROOT, Config, FB, MSUN, G, halo, reservoir
from miniquench.feedback import kernels
from miniquench.analytic import numerical_mass, analytic_mass

plt.style.use(ROOT / "data/apj.mplstyle")
plt.rcParams.update(
    {
        "font.family": "DejaVu Serif",
        "text.usetex": False,
        "axes.labelsize": 12,
        "font.size": 10,
    }
)
t = np.linspace(10, 100, 181)
age = np.arange(0, 100.025, 0.025)
fig, axes = plt.subplots(1, 3, figsize=(12, 3.5), sharey=True)
columns = [t]
names = ["time_myr"]
errors = []
for ax, z in zip(axes, [3, 6, 9]):
    for mc, color, label in [(0.2, "tab:red", "chabrier"), (5, "tab:blue", "topheavy")]:
        cfg = Config(fstar=0.01, z=z, mchar=mc)
        k = kernels(cfg.alpha, mc, cfg.beta, cfg.stellar_radius_power)
        n = reservoir(0, cfg)[2] * cumulative_trapezoid(
            k.cumulative(age, 0), age, initial=0
        )
        mh, rv, _ = halo(age, cfg)
        limit = cfg.mass * (n * 1e51 * rv / (G * mh * MSUN * FB * mh * MSUN)) ** 1.5
        simple = np.interp(t, age, limit)
        original = numerical_mass(t, z, mc, True, 0.05)
        ax.semilogy(t, simple, color=color, label=label + "; constant SFR")
        ax.semilogy(t, original, "--", color=color, label=label + "; reservoir SFH")
        ref = np.genfromtxt(
            ROOT / "data/reference" / f"fig4_z{z}_{label}.csv",
            delimiter=",",
            names=True,
        )
        ax.semilogy(ref["x"], ref["y"], ":", color=color)
        error = np.interp(ref["x"], t, simple) / ref["y"] - 1
        errors.append(
            {
                "z": z,
                "imf": label,
                "prescription": "constant_initial_SFR_and_Mshell_fbMh",
                "median_relative_error": float(np.median(error)),
                "max_absolute_relative_error": float(np.max(abs(error))),
            }
        )
        columns.append(simple)
        names.append(f"z{z}_{label}")
    if z == 6:
        ax.semilogy(t, analytic_mass(t), color="k", lw=0.8, label="Eq.33")
    ax.set(xlabel=r"$t_{SF}$ [Myr]", title=f"z={z}", ylim=(6e9, 2e13), xlim=(10, 100))
axes[0].set_ylabel(r"$M_{h,max}$ [$M_\odot$]")
axes[0].legend(fontsize=7)
fig.suptitle(
    "Fig. 4 conventions: solid = section-3 approximation; dashed = reservoir; dotted = paper",
    fontsize=11,
)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(
        ROOT / "outputs/figures" / f"figure4_conventions.{ext}",
        bbox_inches="tight",
        dpi=180,
    )
np.savetxt(
    ROOT / "data/curves/figure4_constant_sfr.csv",
    np.array(columns).T,
    delimiter=",",
    header=",".join(names),
    comments="",
)
(ROOT / "data/figure4_conventions.json").write_text(json.dumps(errors, indent=2))
print(json.dumps(errors, indent=2))
