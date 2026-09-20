"""Extract a published cooling table and source checksums; no model data generation."""

from pathlib import Path
import re
import json
import hashlib
import subprocess

ROOT = Path(__file__).resolve().parents[1]
source = ROOT / "external_data/schure/source/091025Coolingv4.tex"
text = source.read_text()
rows = []
for line in text.splitlines():
    match = re.match(
        r"^\s*(\d\.\d+)\s*&\s*(-\d+\.\d+)\s*&\s*(-\d+\.\d+)\s*&\s*([\d.E+-]+)\s*\\\\",
        line,
    )
    if match:
        rows.append([float(x) for x in match.groups()])
if len(rows) != 110:
    raise ValueError(f"Expected 110 Schure table 2 entries, got {len(rows)}")
import numpy as np

np.savetxt(
    ROOT / "data/schure_cooling.csv",
    rows,
    delimiter=",",
    header="logT,logLambdaN,logLambdahd,ne_nH",
    comments="",
)
entries = []
for folder, url, version in [
    ("paper", "https://arxiv.org/", "2609.19265v1"),
    (
        "nebrin",
        "https://arxiv.org/src/2409.19288",
        "downloaded latest source; publication 2025",
    ),
    (
        "schure",
        "https://arxiv.org/src/0909.5204",
        "downloaded latest source; publication 2009",
    ),
]:
    for p in sorted((ROOT / "external_data" / folder).rglob("*")):
        if p.is_file():
            entries.append(
                dict(
                    path=str(p.relative_to(ROOT)),
                    source=(
                        (
                            "https://arxiv.org/pdf/2609.19265v1"
                            if p.name in ("paper.pdf", "paper.txt")
                            else "https://arxiv.org/html/2609.19265v1"
                            if p.name == "paper.html"
                            else "https://arxiv.org/src/2609.19265v1"
                        )
                        if folder == "paper"
                        else url
                    ),
                    version=version,
                    sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                    bytes=p.stat().st_size,
                )
            )
repo = ROOT / "third_party/Lyman-alpha-feedback"
commit = subprocess.check_output(
    ["git", "-C", str(repo), "rev-parse", "HEAD"], text=True
).strip()
for p in sorted((repo / "M_F_fit").glob("*.py")):
    entries.append(
        dict(
            path=str(p.relative_to(ROOT)),
            source="https://github.com/olofnebrin/Lyman-alpha-feedback",
            version=commit,
            sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
            bytes=p.stat().st_size,
        )
    )
(ROOT / "data/manifest.json").write_text(
    json.dumps(
        {
            "retrieved_utc": "2026-09-20",
            "paper_version": "2609.19265v1",
            "transport": "direct HTTPS, certificate verification enabled; initial proxy TLS attempt failed without producing data",
            "author_model_code": "No linked code or machine-readable curves found in target paper/source archive",
            "entries": entries,
        },
        indent=2,
    )
)
print(
    f"Extracted {len(rows)} cooling rows; hashed {len(entries)} sources; Nebrin {commit}"
)
