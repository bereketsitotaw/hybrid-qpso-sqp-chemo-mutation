# Hybrid QPSO–SQP Cancer Chemotherapy

[![DOI](https://img.shields.io/badge/DOI-10.1109%2FCDC57313.2025.11312610-blue)](https://doi.org/10.1109/CDC57313.2025.11312610)

MATLAB code for *Optimal Control for Cancer Chemotherapy Using Hybrid Quantum Particle Swarm Optimization* (Kidane, Motayed, Wang).

**Published in:** 2025 IEEE 64th Conference on Decision and Control (CDC), Rio de Janeiro, Brazil, 9–12 December 2025.  
**DOI:** [10.1109/CDC57313.2025.11312610](https://doi.org/10.1109/CDC57313.2025.11312610)

## Quick start

```matlab
cd hybrid-qpso-sqp-chemo-mutation
main                  % all scenarios
main('single')          % Fig. 1 — single drug
main('double')          % three combination therapies
main('figures')         % Figs. 4–6
main('figures', 5)      % Fig. 5 only
```

Results: `results/*.mat`, `results/figures/*.png`

## Single vs double medication

```matlab
addpath('scripts');
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

## Layout (9 MATLAB files)

```
main.m
scripts/
  define_chemotherapy_problem.m   Model + scenarios
  solveOCP.m                      Hybrid QPSO + SQP
  ocp_*.m                         Solver core (5 files)
  chemo_pipeline.m                Run, plot, figures
```

## Requirements

MATLAB R2020a+, Optimization Toolbox (`fmincon`).

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

## License

MIT — see [LICENSE](LICENSE).
