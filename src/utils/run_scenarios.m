function summaries = run_scenarios(scenario_id, varargin)
%RUN_SCENARIOS Solve and plot chemotherapy OCP case studies.
%
%   summaries = run_scenarios('all')
%   summaries = run_scenarios('double')
%   summaries = run_scenarios('cosine_gaussian', 'medication', 'single')

    repro = reproducibility_settings();
    rng(repro.rng_seed, repro.rng_generator);

    repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    results_dir = fullfile(repo_root, 'results');
    figures_dir = fullfile(results_dir, 'figures');
    if ~exist(figures_dir, 'dir')
        mkdir(figures_dir);
    end

    catalog = struct( ...
        'single',           @scenario_single, ...
        'cosine_gaussian',  @scenario_cosine_gaussian, ...
        'cosine_sine',      @scenario_cosine_sine, ...
        'exp_gaussian',     @scenario_exp_gaussian);

    run_list = resolve_run_list(scenario_id, catalog);

    summaries = struct();
    fprintf('Hybrid QPSO-SQP Chemotherapy | rng(%d) | MATLAB %s\n\n', ...
        repro.rng_seed, version);

    for k = 1:numel(run_list)
        tag = run_list{k};
        cfg = catalog.(tag)();
        cfg = apply_overrides(cfg, varargin{:});

        problem = define_chemotherapy_problem(cfg.scenario, ...
            'medication', cfg.medication, ...
            'mutation', cfg.mutation);
        problem.grid.N = repro.N;

        fprintf('--- %s ---\n', cfg.label);
        t_start = tic;
        [solution, summary] = solveOCP(problem);
        summary.total_time = toc(t_start);
        summary.scenario_id = tag;
        summary.label = cfg.label;

        if isfield(problem, 'generate_summary') && isa(problem.generate_summary, 'function_handle')
            problem.generate_summary(solution, summary, problem);
        end

        summaries.(tag) = summary;
        save_chemotherapy_solution(solution, summary, problem, ...
            fullfile(results_dir, [tag '_solution.mat']));
        plot_chemotherapy_results(solution, summary, problem, figures_dir);
    end

    fprintf('\nDone. Results in: %s\n', results_dir);
end

%--------------------------------------------------------------------------
function run_list = resolve_run_list(scenario_id, catalog)
    tags = fieldnames(catalog);
    key = lower(char(scenario_id));
    switch key
        case 'all'
            run_list = tags';
        case 'double'
            run_list = {'cosine_gaussian', 'cosine_sine', 'exp_gaussian'};
        otherwise
            if ~isfield(catalog, key)
                error('Unknown scenario "%s". Use: all, double, %s', ...
                    key, strjoin(tags, ', '));
            end
            run_list = {key};
    end
end

%--------------------------------------------------------------------------
function cfg = apply_overrides(cfg, varargin)
    if isempty(varargin)
        return;
    end
    for i = 1:2:numel(varargin)
        name = lower(char(varargin{i}));
        value = varargin{i + 1};
        switch name
            case 'medication'
                cfg.medication = lower(char(value));
            case 'mutation'
                cfg.mutation = logical(value);
            otherwise
                error('Unknown override "%s".', name);
        end
    end
end
