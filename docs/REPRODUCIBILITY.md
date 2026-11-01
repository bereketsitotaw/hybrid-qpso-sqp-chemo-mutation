# Reproducibility

Companion notes for *Optimal Control for Cancer Chemotherapy Using Hybrid Quantum Particle Swarm Optimization* (IEEE CDC 2025).

**DOI:** [10.1109/CDC57313.2025.11312610](https://doi.org/10.1109/CDC57313.2025.11312610)

## Environment

- MATLAB R2020a+ (tested R2024a / R2024b)
- Optimization Toolbox (`fmincon`, SQP)
- Parallel Computing Toolbox optional (`parfor` in QPSO; `UseParallel` in `fmincon`)

Pinned settings: `configs/reproducibility_settings.m`.

## Cost weights (critical)

Publication defaults are CDC 2025 Section V.A:

| Symbol | Value |
|--------|-------|
| α (terminal) | 5 |
| β (running tumor) | 400 |
| γ₁ (drug 1) | 4000 |
| γ₂ (drug 2) | 4000 |

Implemented in `configs/paper_cost_weights.m` and enforced when `weight_profile='paper'` (default).

Do **not** use `weight_profile='experimental'` for paper reproduction. That profile is a **non-paper exploratory** archive of a local trial (`configs/experimental_cost_weights.m`) — **not validated**, not publication defaults, and kept only for forensics / what-if experiments.

## Determinism

- `chemo_pipeline` calls `rng(42,'twister')` before solves.
- QPSO still uses stochastic updates; identical hardware/toolbox versions are required for bit-stability.
- Parallel pool size can change wall-clock and (rarely) floating-point reduction order. For stricter checks, set `problem.options.nlp.UseParallel = false`.

## Reproduce paper cases

```matlab
cd hybrid-qpso-sqp-chemo-mutation
main('single')            % Fig. 1
main('double')            % combination therapies
main('figures')           % Figs. 4–6 from cached/generated results
```

Outputs:

- `results/generated/*_solution.mat` (preferred)
- `results/*_solution.mat` (legacy mirror)
- `results/figures/*.png`

## Reference pin

`results/reference/single_solution.mat` is a pinned single-drug solve under paper weights:

- `T = 10`, `exitflag = 1`
- `cost ≈ 270530.13`
- final average tumor ≈ 2.158
- drug-1 exposure ≈ 25.19

Metrics are also listed in `configs/reference_metrics.m`. Compare new runs against this pin; tolerances are in `reproducibility_settings.m`.

## Smoke / tests (not paper-quality)

```matlab
main('test')    % structural + weight + reference checks
main('smoke')   % tiny QPSO/SQP solve for plumbing only
```

## CI

GitHub Actions workflow `.github/workflows/repo-smoke.yml` checks layout and that paper weight constants remain present in configs. Full MATLAB solves are local (license required).
