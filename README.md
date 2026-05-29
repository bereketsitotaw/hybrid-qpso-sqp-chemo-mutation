# Hybrid QPSO–SQP Cancer Chemotherapy Code

[![MATLAB](https://img.shields.io/badge/MATLAB-R2024a-orange)](requirements.txt)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

MATLAB code for the hybrid **Quantum PSO + SQP** method in:

> *Optimal Control for Cancer Chemotherapy Using Hybrid Quantum Particle Swarm Optimization*  
> B. S. Kidane, M. S. H. Motayed, S. Wang (UT Arlington)

**Authors:** Bereket Sitotaw Kidane, Md Samiul Haque Motayed, Shuo Wang

## Quick start

```matlab
cd hybrid-qpso-sqp-chemotherapy
main                  % all scenarios (single + three double-drug cases)
main('single')        % Fig. 1 — single medication
main('double')        % combination therapy only (Sec. V.B)
reproduce_paper_figures('all')   % Figs. 1–6
reproduce_paper_figures([5 6])  % comparison figures only
```

Outputs: `results/*_solution.mat`, `results/figures/*.png`

## Single vs double medication

One API handles both modes:

```matlab
problem = define_chemotherapy_problem('single');           % nu = 1
problem = define_chemotherapy_problem('cosine_gaussian');  % nu = 2
problem = define_chemotherapy_problem('exp_gaussian', 'phi_scale', 2);
[solution, summary] = solveOCP(problem);
```

| Command | Paper | Drugs |
|---------|-------|-------|
| `main('single')` | Fig. 1 | 1 |
| `main('cosine_gaussian')` | Fig. 2 / Sec. V.B.1 | 2 |
| `main('cosine_sine')` | Fig. 3 / Sec. V.B.2 | 2 |
| `main('exp_gaussian')` | Fig. 4 / Sec. V.B.3 | 2 |
| `reproduce_paper_figures(5)` | Fig. 5 — final trait distributions | — |
| `reproduce_paper_figures(6)` | Fig. 6 — population evolution | — |

## Requirements

- MATLAB R2020a+ (R2024a tested)
- Optimization Toolbox (`fmincon`, SQP)

## Repository layout

```
main.m                          Entry point
reproduce_paper_figures.m         Figs. 4–6 (and optional 1–3)
init_project.m                    Path setup
define_chemotherapy_problem.m     Trait OCP (single/double drug)
solveOCP.m                        Hybrid QPSO + SQP solver
configs/                          Scenario parameters
src/utils/                        Run pipeline, plots, figure scripts
results/                          Generated locally (gitignored)
```

## Related work

Chemo–immunotherapy companion code (3-state model, Cases A/B/C):  
[github.com/bereketsitotaw/hybrid-qpso-sqp-chemo-immuno](https://github.com/bereketsitotaw/hybrid-qpso-sqp-chemo-immuno) · [doi:10.3934/dcdsb.2026072](https://doi.org/10.3934/dcdsb.2026072)

## Citation

```bibtex
@article{kidane2026chemotherapy,
  title   = {Optimal Control for Cancer Chemotherapy Using Hybrid Quantum Particle Swarm Optimization},
  author  = {Kidane, Bereket Sitotaw and Motayed, Md Samiul Haque and Wang, Shuo},
  year    = {2026}
}
```

## License

MIT — see [LICENSE](LICENSE).
