# Independent mini-quenching reproduction

This repository is independent of AuroraLF. Work from this repository root.
Read REVIEW_START_HERE.md and reproduction_report.md before modifying the model.
Use the fixed paper arXiv:2609.19265v1 as the primary physical specification.
Use the repository's own .venv explicitly and retain dependency provenance.

Preserve source data, prior outputs, failures and uncommitted changes. Keep future
results separate when changing physics. Do not label numerical convergence as
successful reproduction. Do not tune parameters to match reference curves or
silently replace missing inputs. Distinguish author definitions from inferred
closures, diagnostic branches and implementation errors.

For model changes run the affected tests, the two Figure 2 single trajectories,
and appropriate convergence checks before a larger scan. Use CPU only, at most
four workers, no nested parallelism. Do not run unrelated AuroraLF workflows.
Keep a Chinese report and left-original/right-computed review figures. Follow
explicit user authorization for commits, pushes and repository visibility.
