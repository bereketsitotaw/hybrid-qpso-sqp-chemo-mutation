# Hybrid QPSO–SQP Cancer Chemotherapy

[![DOI](https://img.shields.io/badge/DOI-10.1109%2FCDC57313.2025.11312610-blue)](https://doi.org/10.1109/CDC57313.2025.11312610)

MATLAB code for *Optimal Control for Cancer Chemotherapy Using Hybrid Quantum Particle Swarm Optimization* (Kidane, Motayed, Wang).

**Published in:** 2025 IEEE 64th Conference on Decision and Control (CDC), Rio de Janeiro, Brazil, 9–12 December 2025.  
**DOI:** [10.1109/CDC57313.2025.11312610](https://doi.org/10.1109/CDC57313.2025.11312610)

## Quick start

```matlab
cd hybrid-qpso-sqp-chemo-mutation
main                  % all scenarios
main('single')        % Fig. 1 — single drug
main('double')        % three combination therapies
main('figures')       % Figs. 4–6
main('figures', 5)    % Fig. 5 only
main('test')          % smoke tests
main('smoke')         % tiny non-paper plumbing solve
```

Results: `results/generated/*.mat`, `results/figures/*.png` (legacy mirror also written to `results/*.mat`).

## Method (overview)

Hybrid **QPSO** global search warm-starts a trapezoidal collocation NLP refined by **`fmincon` (SQP)**. Optional PMP post-processing checks Hamiltonian stationarity. Tumor trait dynamics use Michaelis–Menten cytotoxicity with optional mutation (Section II / V of the paper).

## Cost weights (Section V.A)

Publication defaults (do not change for paper reproduction):

- α = 5, β = 400, γ₁ = 4000, γ₂ = 4000  

Defined in `configs/paper_cost_weights.m`. Override only via `define_chemotherapy_problem(..., 'weight_profile', ...)`.

## Single vs double medication

```matlab
addpath('scripts'); addpath('configs');
problem = define_chemotherapy_problem('single');
problem = define_chemotherapy_problem('cosine_gaussian');
[solution, summary] = solveOCP(problem);
```

| `main(...)` | Paper |
|-------------|-------|
| `'single'` | Fig. 1 |
| `'cosine_gaussian'` | Fig. 2 |
| `'cosine_sine'` | Fig. 3 |
| `'exp_gaussian'` | Fig. 4 |
| `'figures', 5` | Fig. 5 |
| `'figures', 6` | Fig. 6 |

## Repository layout

```
main.m
configs/          Paper weights, reproducibility, reference metrics
scripts/          Model + hybrid solver (lean MATLAB core)
tests/            Smoke tests (main('test'))
docs/             REPRODUCIBILITY.md
results/
  reference/      Pinned single-drug .mat + notes
  generated/      Local solves (gitignored)
  figures/        Generated plots (gitignored)
.github/workflows/repo-smoke.yml
```

## Requirements

MATLAB R2020a+, Optimization Toolbox (`fmincon`). Parallel Computing Toolbox optional.

## Reproducing results

See [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md). Compare new single-drug runs to `results/reference/single_solution.mat` / `configs/reference_metrics.m`.

## Citation

If you use this code, please cite:

> B. S. Kidane, M. S. H. Motayed and S. Wang, "Optimal Control for Cancer Chemotherapy Using Hybrid Quantum Particle Swarm Optimization," in *2025 IEEE 64th Conference on Decision and Control (CDC)*, Rio de Janeiro, Brazil, 2025, doi: [10.1109/CDC57313.2025.11312610](https://doi.org/10.1109/CDC57313.2025.11312610).

```bibtex
@inproceedings{kidane2025chemotherapy,
  author    = {Kidane, Bereket Sitotaw and Motayed, Md Samiul Haque and Wang, Shuo},
  title     = {Optimal Control for Cancer Chemotherapy Using Hybrid Quantum Particle Swarm Optimization},
  booktitle = {2025 IEEE 64th Conference on Decision and Control (CDC)},
  address   = {Rio de Janeiro, Brazil},
  year      = {2025},
  month     = dec,
  publisher = {IEEE},
  doi       = {10.1109/CDC57313.2025.11312610},
  isbn      = {978-8-3315-2627-6},
  issn      = {2576-2370}
}
```

**Publication details**

| Field | Value |
|-------|--------|
| Conference | 2025 IEEE 64th CDC |
| Dates | 9–12 December 2025 |
| Location | Rio de Janeiro, Brazil |
| DOI | [10.1109/CDC57313.2025.11312610](https://doi.org/10.1109/CDC57313.2025.11312610) |
| Electronic ISBN | 978-8-3315-2627-6 |
| Print ISBN (PoD) | 978-8-3315-2628-3 |
| Electronic ISSN | 2576-2370 |
| Print ISSN (PoD) | 0743-1546 |
| IEEE Xplore | Added 12 January 2026 |

## Contact / Issues

Open issues at the GitHub repository linked in [CITATION.cff](CITATION.cff).

## License

MIT — see [LICENSE](LICENSE).
