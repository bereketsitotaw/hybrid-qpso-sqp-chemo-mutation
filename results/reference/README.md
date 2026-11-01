# Reference outputs

Pinned artifacts for reproducibility checks (paper cost weights).

## Pinned single-drug solve

File: `single_solution.mat` (MATLAB 5.0, MACI64, created 2026-05-29)

| Metric | Value |
|--------|-------|
| Scenario | `single` (Fig. 1) |
| T | 10 |
| exitflag | 1 |
| cost | 270530.1263 |
| final average tumor | ≈ 2.157611315 |
| drug-1 exposure ∫u₁ dt | ≈ 25.18831062 |
| Weights | α=5, β=400, γ₁=4000 |

Also encoded in `configs/reference_metrics.m`.

## Expected generated artifacts (local)

Run `main('all')` to populate:

- `results/generated/*_solution.mat` (and legacy `results/*.mat`)
- `results/figures/case_*.png`, `figures/pmp_*.png`, Figs. 4–6

Do not overwrite files in this `reference/` directory with experimental-weight runs.
