"""Extract plotted vector vertices, not underlying author simulation arrays."""

from pathlib import Path
import json
import numpy as np
import pymupdf

ROOT = Path(__file__).resolve().parents[1]
# Axes calibration visually verified against original PDF tick positions.
SPECS = {
    "IMF_plot": [
        (1, f"alpha{a}_mc{m}", i, 192, 193, (np.log10(0.08), 2), (-5, 1), True, True)
        for a, m, i in [
            (1.35, 0.2, 171),
            (1.35, 5, 172),
            (1.35, 10, 173),
            (2.35, 0.2, 174),
            (2.35, 5, 175),
            (2.35, 10, 176),
        ]
    ],
    "tquench_exs": [
        (2, "m8", 91, 105, 106, (0, 230), (0.7, 11), False, False),
        (2, "m9", 92, 105, 106, (0, 230), (0.7, 11), False, False),
    ],
    "mh_max_tmax": [
        (4, "z3_chabrier", 95, 97, 98, (10, 100), (10, 13), False, True),
        (4, "z3_topheavy", 96, 97, 98, (10, 100), (10, 13), False, True),
        (4, "z6_chabrier", 198, 201, 202, (10, 100), (10, 13), False, True),
        (4, "z6_topheavy", 199, 201, 202, (10, 100), (10, 13), False, True),
        (4, "z6_analytic", 200, 201, 202, (10, 100), (10, 13), False, True),
        (4, "z9_chabrier", 301, 303, 304, (10, 100), (10, 13), False, True),
        (4, "z9_topheavy", 302, 303, 304, (10, 100), (10, 13), False, True),
    ],
    "tquench_IMFs": [
        (5, "f001", 133, 138, 139, (8, 13), (1, 200), True, False),
        (5, "f01", 134, 138, 139, (8, 13), (1, 200), True, False),
        (5, "f1", 135, 138, 139, (8, 13), (1, 200), True, False),
        (5, "mc02", 491, 484, 485, (8, 13), (1, 200), True, False),
        (5, "mc1", 482, 484, 485, (8, 13), (1, 200), True, False),
        (5, "mc5", 483, 484, 485, (8, 13), (1, 200), True, False),
    ],
    "tquench_diff_fcovers": [
        (8, "f001_lowgas", 1, 154, 155, (8, 13), (-5, 275), True, False),
        (8, "f001_c1", 150, 154, 155, (8, 13), (-5, 275), True, False),
        (8, "f001_c01", 151, 154, 155, (8, 13), (-5, 275), True, False),
        (8, "f001_c005", 152, 154, 155, (8, 13), (-5, 275), True, False),
        (8, "f001_c0005", 153, 154, 155, (8, 13), (-5, 275), True, False),
        (8, "f01_lowgas", 163, 316, 317, (8, 13), (-5, 275), True, False),
        (8, "f01_c1", 312, 316, 317, (8, 13), (-5, 275), True, False),
        (8, "f01_c01", 313, 316, 317, (8, 13), (-5, 275), True, False),
        (8, "f01_c005", 314, 316, 317, (8, 13), (-5, 275), True, False),
        (8, "f01_c0005", 315, 316, 317, (8, 13), (-5, 275), True, False),
    ],
}
manifest = []
for stem, specs in SPECS.items():
    path = ROOT / "external_data/paper/source/plots" / f"{stem}.pdf"
    page = pymupdf.open(path)[0]
    drawings = page.get_drawings()
    for fig, label, index, left, right, xlim, ylim, logx, logy in specs:
        d = drawings[index]
        left_axis = drawings[left]["rect"]
        r = drawings[right]["rect"]
        points = []
        for item in d["items"]:
            if item[0] != "l":
                raise ValueError("Expected straight-line plotted polyline")
            if not points:
                points.append(tuple(item[1]))
            points.append(tuple(item[2]))
        p = np.array(points)
        x = xlim[0] + (p[:, 0] - left_axis.x0) / (r.x0 - left_axis.x0) * (
            xlim[1] - xlim[0]
        )
        y = ylim[0] + (left_axis.y1 - p[:, 1]) / (left_axis.y1 - left_axis.y0) * (
            ylim[1] - ylim[0]
        )
        if logx:
            x = 10**x
        if logy:
            y = 10**y
        good = (
            (p[:, 0] >= left_axis.x0 - 0.01)
            & (p[:, 0] <= r.x0 + 0.01)
            & (p[:, 1] >= left_axis.y0)
            & (p[:, 1] <= left_axis.y1)
        )
        out = ROOT / "data/reference" / f"fig{fig}_{label}.csv"
        np.savetxt(
            out, np.c_[x[good], y[good]], delimiter=",", header="x,y", comments=""
        )
        manifest.append(
            dict(
                file=str(out.relative_to(ROOT)),
                source=str(path.relative_to(ROOT)),
                drawing=index,
                axis_pdf_points=[left_axis.x0, r.x0, left_axis.y0, left_axis.y1],
                xlimits=xlim,
                ylimits=ylim,
                logx=logx,
                logy=logy,
                kind="published_vector_plot_vertices",
                coordinate_uncertainty_pdf_points=0.1,
                notes="Clipped vertices discarded; piecewise-linear drawing, not original solver samples. Adopt conservative 2 Myr timing and 0.03 dex mass reading uncertainty.",
            )
        )
(ROOT / "data/reference/manifest.json").write_text(json.dumps(manifest, indent=2))
print("Extracted", len(manifest), "independent reference curves")
