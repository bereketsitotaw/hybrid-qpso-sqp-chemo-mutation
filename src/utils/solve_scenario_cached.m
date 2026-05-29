function [solution, summary, problem] = solve_scenario_cached(scenario_id, results_dir, varargin)
%SOLVE_SCENARIO_CACHED Load saved solution or run solveOCP if missing.

    if nargin < 2 || isempty(results_dir)
        repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
        results_dir = fullfile(repo_root, 'results');
    end
    if ~exist(results_dir, 'dir')
        mkdir(results_dir);
    end

    mat_path = fullfile(results_dir, [scenario_id '_solution.mat']);
    if exist(mat_path, 'file')
        data = load(mat_path);
        solution = data.solution;
        summary = data.summary;
        problem = define_chemotherapy_problem(scenario_id, varargin{:});
        return;
    end

    repro = reproducibility_settings();
    rng(repro.rng_seed, repro.rng_generator);

    problem = define_chemotherapy_problem(scenario_id, varargin{:});
    problem.grid.N = repro.N;
    problem.options.enable_pmp = false;

    fprintf('Solving scenario "%s" (no cache)...\n', scenario_id);
    t_start = tic;
    [solution, summary] = solveOCP(problem);
    summary.total_time = toc(t_start);
    save_chemotherapy_solution(solution, summary, problem, mat_path);
end
