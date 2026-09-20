"""Readable review-only views of already computed curves; no new simulation."""

import json
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from miniquench.physics import ROOT

plt.style.use(ROOT / "data/apj.mplstyle")
plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "DejaVu Serif",
        "axes.labelsize": 13,
        "font.size": 11,
    }
)
for model in ["sn_only", "full_approx"]:
    fig, ax = plt.subplots(figsize=(8, 4.3))
    for mass, color, tag in [(1e8, "tab:blue", "m8"), (1e9, "tab:red", "m9")]:
        ref = np.genfromtxt(
            ROOT / "data/reference" / f"fig2_{tag}.csv", delimiter=",", names=True
        )
        ax.plot(ref["x"], ref["y"], ":", color=color, label=f"Paper: {mass:.0e} Msun")
        path = ROOT / "data/curves" / f"fig2_{mass:.0e}_full_{model}.csv"
        if path.exists():
            d = np.genfromtxt(path, delimiter=",", names=True)
            ax.plot(
                d["time_myr"],
                d["r_over_R0"],
                color=color,
                label=f"{model}: {mass:.0e} Msun",
            )
        else:
            ax.text(
                0.55,
                0.7,
                "1e9 Msun: unresolved launch",
                color=color,
                transform=ax.transAxes,
            )
    ax.set(xlabel="Time [Myr]", ylabel=r"$r/R_0(t)$", ylim=(0.8, 9.5))
    ax.legend(fontsize=9, ncol=2)
    fig.tight_layout()
    fig.savefig(
        ROOT / "slides/assets" / f"review_fig2_{model}.pdf", bbox_inches="tight"
    )
    fig.savefig(
        ROOT / "outputs/figures" / f"review_fig2_{model}.png",
        bbox_inches="tight",
        dpi=160,
    )
    plt.close(fig)
d = json.loads((ROOT / "data/curves/summary_full_sn_only_258.json").read_text())
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
for ax, items in zip(
    axes,
    [
        [
            ("f001", r"$f_\star=0.001$", "tab:blue"),
            ("f01", r"$f_\star=0.01$", "tab:purple"),
            ("f1", r"$f_\star=0.1$", "tab:red"),
        ],
        [
            ("f01", r"$m_{char}=0.2$", "tab:purple"),
            ("mc1", r"$m_{char}=1$", "darkorange"),
            ("mc5", r"$m_{char}=5$", "seagreen"),
        ],
    ],
):
    for key, label, color in items:
        rows = sorted(d["results"]["fig5_" + key], key=lambda r: r["config"]["mass"])
        x = [r["config"]["mass"] for r in rows]
        y = [
            r["t_quench_myr"] if r["status"] == "completed_cycle" else np.nan
            for r in rows
        ]
        ax.semilogx(x, y, color=color, label=label)
        ref = np.genfromtxt(
            ROOT / "data/reference" / f"fig5_{key}.csv", delimiter=",", names=True
        )
        ax.semilogx(ref["x"], ref["y"], ":", color=color)
        for r in rows:
            if r["status"] == "launched_not_returned_by_tmax":
                ax.plot(
                    r["config"]["mass"], r["t_quench_lower_bound_myr"], "^", color=color
                )
    ax.set(
        xlabel=r"Initial halo mass [$M_\odot$]",
        ylabel="Quenching duration [Myr]",
        ylim=(0, 730),
        xlim=(1e8, 1e13),
    )
    ax.legend(fontsize=10)
axes[0].set_title("Efficiency")
axes[1].set_title("IMF at fstar=0.01")
fig.suptitle(
    "SN-only: solid = computed, dotted = paper; gaps are not successful cycles"
)
fig.tight_layout()
fig.savefig(ROOT / "slides/assets/review_fig5_sn_only.pdf", bbox_inches="tight")
fig.savefig(
    ROOT / "outputs/figures/review_fig5_sn_only.png", bbox_inches="tight", dpi=160
)
