"""Compare independent published vector paths; never fit or tune model inputs."""

import json
import numpy as np
from miniquench.physics import ROOT
from miniquench.feedback import IMF

records = []
for alpha in (1.35, 2.35):
    for mc in (0.2, 5, 10):
        ref = np.genfromtxt(
            ROOT / "data/reference" / f"fig1_alpha{alpha}_mc{mc:g}.csv",
            delimiter=",",
            names=True,
        )
        pred = IMF(alpha, mc)(ref["x"])
        err = pred / ref["y"] - 1
        records.append(
            {
                "figure": 1,
                "curve": f"alpha{alpha}_mc{mc}",
                "n": len(err),
                "median_relative_error": float(np.median(err)),
                "max_absolute_relative_error": float(np.max(np.abs(err))),
            }
        )
a = np.genfromtxt(ROOT / "data/curves/figure4.csv", delimiter=",", names=True)
for z in (3, 6, 9):
    for imf in ("chabrier", "topheavy"):
        ref = np.genfromtxt(
            ROOT / "data/reference" / f"fig4_z{z}_{imf}.csv", delimiter=",", names=True
        )
        pred = np.interp(ref["x"], a["time_myr"], a[f"z{z}_{imf}_growth1"])
        err = pred / ref["y"] - 1
        records.append(
            {
                "figure": 4,
                "curve": f"z{z}_{imf}",
                "n": len(err),
                "median_relative_error": float(np.median(err)),
                "max_absolute_relative_error": float(np.max(np.abs(err))),
            }
        )
for path in sorted((ROOT / "data/curves").glob("summary_full*.json")):
    d = json.loads(path.read_text())
    for family, rows in d["results"].items():
        if not family.startswith(("fig5", "fig8")):
            continue
        fig = int(family[3])
        label = family[5:]
        ref = np.genfromtxt(
            ROOT / "data/reference" / f"fig{fig}_{label}.csv", delimiter=",", names=True
        )
        diff = []
        frac = []
        for r in rows:
            m = r["config"]["mass"]
            if (
                r["status"] != "completed_cycle"
                or not ref["x"].min() <= m <= ref["x"].max()
            ):
                continue
            y = float(np.interp(np.log10(m), np.log10(ref["x"]), ref["y"]))
            if y <= 2:
                continue  # below adopted 2 Myr readout precision, omit relative-error metric only
            delta = r["t_quench_myr"] - y
            diff.append(delta)
            frac.append(delta / y)
        records.append(
            {
                "figure": fig,
                "model": d["model"],
                "loading": d["loading"],
                "curve": label,
                "n_valid_paired": len(diff),
                "n_total": len(rows),
                "median_difference_myr": float(np.median(diff)) if diff else None,
                "median_relative_error": float(np.median(frac)) if frac else None,
                "note": "Only completed, overlapping, positive reference points; missing/failed runs are NOT agreement.",
            }
        )
(ROOT / "data/paper_comparison.json").write_text(json.dumps(records, indent=2))
for r in records:
    if r["figure"] in (1, 4):
        print(r)
