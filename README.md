# Hybrid QPSO–SQP Cancer Chemotherapy

MATLAB code for *Optimal Control for Cancer Chemotherapy Using Hybrid Quantum Particle Swarm Optimization* (Kidane, Motayed, Wang).

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

## License

MIT — see [LICENSE](LICENSE).
