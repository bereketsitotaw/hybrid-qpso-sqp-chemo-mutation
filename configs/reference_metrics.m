function metrics = reference_metrics()
%REFERENCE_METRICS  Pinned single-drug reference from results/reference/.
%
%   Source file: results/reference/single_solution.mat
%   Generated: 2026-05-29 (MACI64) under paper cost weights (alpha=5, beta=400,
%   gamma1=4000), exitflag=1. Copied from a clean duplicate checkout during the
%   deep reproducibility pass. Re-validate after any scientific change.

    metrics = struct();
    metrics.scenario = 'single';
    metrics.mat_file = fullfile('results', 'reference', 'single_solution.mat');
    metrics.T_final = 10;
    metrics.cost = 270530.1263;
    metrics.exitflag = 1;
    metrics.final_avg_tumor = 2.157611315;
    metrics.dose1 = 25.18831062;
    metrics.state_size = [21, 101];
    metrics.control_size = [1, 101];
    metrics.weights = struct('alpha', 5, 'beta', 400, 'gamma1', 4000);
    metrics.notes = 'Paper Section V.A weights; Fig. 1 single-drug case.';
end
