"""Render saved results without recomputing trajectories."""

import json
import importlib.util
import numpy as np
from miniquench.physics import ROOT

spec = importlib.util.spec_from_file_location("plotting", ROOT / "scripts/reproduce.py")
p = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p)
for model, suffix, figures in [
    ("full_approx", "12458", [2, 5, 8]),
    ("sn_only", "258", [2, 5, 8]),
]:
    d = json.loads(
        (ROOT / "data/curves" / f"summary_full_{model}_{suffix}.json").read_text()
    )
    curves = {}
    for mass in (1e8, 1e9):
        source = ROOT / "data/curves" / f"fig2_{mass:.0e}_full_{model}.csv"
        result = next(r for r in d["results"]["fig2"] if r["config"]["mass"] == mass)
        if result["status"] == "numerical_failure":
            continue
        a = np.genfromtxt(source, delimiter=",", names=True)
        curves["fig2", mass] = {k: a[k] for k in a.dtype.names}
    p.MODEL_TAG = model
    p.plot_results(d["results"], curves, figures, model)
