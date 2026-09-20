"""Check stored evidence, code fingerprints, source hashes and real cache reuse."""

import hashlib
import json
import time
from collections import Counter
from miniquench.physics import ROOT, Config, code_fingerprint
from miniquench.dynamics import run_cached

fingerprint = code_fingerprint()
checks = {}
for p in sorted((ROOT / "data/curves").glob("summary_full*.json")):
    data = json.loads(p.read_text())
    assert data["code_fingerprint"] == fingerprint, (p, "stale code fingerprint")
    rows = sum(data["results"].values(), [])
    checks[p.name] = {
        "total": len(rows),
        "statuses": dict(Counter(r["status"] for r in rows)),
    }
    for r in rows:
        if r["status"] == "completed_cycle":
            assert r["t_quench_myr"] > 0
        else:
            assert r["t_quench_myr"] is None
        if r["status"] == "launched_not_returned_by_tmax":
            assert r["t_quench_lower_bound_myr"] > 0
manifest = json.loads((ROOT / "data/manifest.json").read_text())
for e in manifest["entries"]:
    assert hashlib.sha256((ROOT / e["path"]).read_bytes()).hexdigest() == e["sha256"], (
        e["path"]
    )
validation = json.loads((ROOT / "data/validation.json").read_text())
assert validation["code_fingerprint"] == fingerprint
assert all(r["all_variations_within_1_percent"] for r in validation["convergence"])
start = time.perf_counter()
r, _, hit = run_cached(Config(mass=1e8, loading="until_return", history_dt=0.05))
assert hit
checks["cache_resume"] = {
    "cache_hit": hit,
    "wall_seconds": time.perf_counter() - start,
    "status": r["status"],
}
checks["source_hashes_verified"] = len(manifest["entries"])
checks["code_fingerprint"] = fingerprint
(ROOT / "data/artifact_checks.json").write_text(json.dumps(checks, indent=2))
print(json.dumps(checks, indent=2))
