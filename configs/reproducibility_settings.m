function repro = reproducibility_settings()
%REPRODUCIBILITY_SETTINGS  Pinned solver / RNG settings for this repository.
%
%   Pipeline seed is set in chemo_pipeline via rng(42,'twister').
%   Full paper runs use the QPSO sizes below; smoke tests may override down.

    repro = struct();

    repro.matlab_version_tested = 'R2024a / R2024b';
    repro.doi = '10.1109/CDC57313.2025.11312610';

    repro.rng_seed = 42;
    repro.rng_generator = 'twister';

    % Defaults embedded in define_chemotherapy_problem (paper-scale search)
    repro.qpso_swarm_size = 150;
    repro.qpso_max_iter = 250;
    repro.grid_N = 100;
    repro.T_fixed = 10.0;
    repro.trait_nodes = 21;

    % fmincon / SQP (solveOCP defaults)
    repro.sqp_algorithm = 'sqp';
    repro.sqp_max_iterations = 5000;
    repro.sqp_max_function_evaluations = 1e6;
    repro.sqp_constraint_tolerance = 1e-6;
    repro.sqp_optimality_tolerance = 1e-6;
    repro.sqp_use_parallel_default = true;  % requires Parallel Computing Toolbox

    % Smoke overrides (tests / CI-friendly; not paper-quality optima)
    repro.smoke_qpso_swarm_size = 8;
    repro.smoke_qpso_max_iter = 5;
    repro.smoke_grid_N = 20;
    repro.smoke_enable_pmp = false;
    repro.smoke_use_parallel = false;

    % Reference comparison tolerances (single-drug pinned .mat)
    repro.ref_cost_rel_tol = 0.05;
    repro.ref_final_avg_abs_tol = 0.5;
    repro.ref_dose1_abs_tol = 2.0;
end
