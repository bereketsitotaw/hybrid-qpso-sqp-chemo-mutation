# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Added
- `CHANGELOG.md`
- `configs/` — `paper_cost_weights.m`, `experimental_cost_weights.m`, `reproducibility_settings.m`, `reference_metrics.m`
- `tests/` smoke suite (`main('test')`)
- `docs/REPRODUCIBILITY.md`
- `.github/workflows/repo-smoke.yml` layout/weight guards
- `results/reference/single_solution.mat` pinned single-drug artifact
- `results/generated/` output directory
- `main('smoke')` / `main('test')` entry points

### Changed
- Restored CDC 2025 Section V.A cost weights as the enforced default (`α=5`, `β=400`, `γ₁=4000`)
- Expanded `.gitignore` for MATLAB / OS artifacts and generated results
- README: method, weights, tests, reproducibility, layout
- Solver summary now prints active weight profile
- Generated solves also written under `results/generated/`

### Fixed
- Local uncommitted drift that had set `α=1`, `β=0.1`, `γ₁=20` (objective-changing; not paper)
- Redundant control clipping expressions after QPSO smoothing in `solveOCP.m`
