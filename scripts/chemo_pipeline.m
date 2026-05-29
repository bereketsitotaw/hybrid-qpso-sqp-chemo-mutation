function summaries = chemo_pipeline(repo_root, action, arg, varargin)
%CHEMO_PIPELINE Run scenarios, save results, and plot paper figures.
%
%   chemo_pipeline(repo_root, 'run', 'all')
%   chemo_pipeline(repo_root, 'figures', [5 6])

    rng(42, 'twister');
    N_grid = 100;

    results_dir = fullfile(repo_root, 'results');
    figures_dir = fullfile(results_dir, 'figures');
    if ~exist(figures_dir, 'dir'), mkdir(figures_dir); end

    switch lower(char(action))
        case 'run'
            summaries = run_all(repo_root, results_dir, figures_dir, arg, N_grid, varargin{:});
        case 'figures'
            if nargin < 3 || isempty(arg), arg = [4, 5, 6]; end
            summaries = [];
            make_figures(repo_root, results_dir, figures_dir, arg, N_grid);
        otherwise
            error('Unknown action "%s". Use ''run'' or ''figures''.', action);
    end
end

%==========================================================================
function summaries = run_all(repo_root, results_dir, figures_dir, scenario_id, N_grid, varargin)
    catalog = scenario_catalog();
    run_list = resolve_list(scenario_id, catalog);

    summaries = struct();
    fprintf('Hybrid QPSO-SQP Chemotherapy | rng(42)\n\n');

    for k = 1:numel(run_list)
        tag = run_list{k};
        cfg = catalog.(tag);
        cfg = apply_opts(cfg, varargin{:});

        problem = define_chemotherapy_problem(cfg.scenario, ...
            'medication', cfg.medication, 'mutation', cfg.mutation);
        problem.grid.N = N_grid;

        fprintf('--- %s ---\n', cfg.label);
        t0 = tic;
        [solution, summary] = solveOCP(problem);
        summary.total_time = toc(t0);
        summary.label = cfg.label;

        if isfield(problem, 'generate_summary')
            problem.generate_summary(solution, summary, problem);
        end

        summaries.(tag) = summary;
        save(fullfile(results_dir, [tag '_solution.mat']), 'solution', 'summary');
        plot_case(solution, problem, figures_dir);
    end
    fprintf('\nDone. Results: %s\n', results_dir);
end

%==========================================================================
function make_figures(~, results_dir, figures_dir, which, N_grid)
    if ischar(which) && strcmpi(which, 'all')
        run_all(fileparts(results_dir), results_dir, figures_dir, 'all', N_grid);
        which = [4, 5, 6];
    end
    for f = which(:)'
        switch f
            case 4, fig4(results_dir, figures_dir, N_grid);
            case 5, fig5(results_dir, figures_dir);
            case 6, fig6(results_dir, figures_dir);
            otherwise, warning('Skip unknown figure %d.', f);
        end
    end
    fprintf('Figures saved to %s\n', figures_dir);
end

%==========================================================================
function catalog = scenario_catalog()
    catalog.single = struct('scenario', 'single', 'medication', 'single', ...
        'mutation', false, 'label', 'Single-drug PMP benchmark (Fig. 1)');
    catalog.cosine_gaussian = struct('scenario', 'cosine_gaussian', 'medication', 'double', ...
        'mutation', true, 'label', 'Cosine-Gaussian (Sec. V.B.1)');
    catalog.cosine_sine = struct('scenario', 'cosine_sine', 'medication', 'double', ...
        'mutation', true, 'label', 'Cosine-Sine (Sec. V.B.2)');
    catalog.exp_gaussian = struct('scenario', 'exp_gaussian', 'medication', 'double', ...
        'mutation', true, 'label', 'Exponential-Gaussian (Sec. V.B.3)');
end

function run_list = resolve_list(scenario_id, catalog)
    tags = fieldnames(catalog);
    key = lower(char(scenario_id));
    switch key
        case 'all',  run_list = tags';
        case 'double', run_list = {'cosine_gaussian', 'cosine_sine', 'exp_gaussian'};
        otherwise
            if ~isfield(catalog, key)
                error('Unknown scenario "%s". Use: all, double, %s', key, strjoin(tags, ', '));
            end
            run_list = {key};
    end
end

function cfg = apply_opts(cfg, varargin)
    for i = 1:2:numel(varargin)
        switch lower(char(varargin{i}))
            case 'medication', cfg.medication = lower(char(varargin{i+1}));
            case 'mutation',   cfg.mutation = logical(varargin{i+1});
        end
    end
end

%==========================================================================
function [sol, sumy, prob] = load_or_solve(tag, results_dir, N_grid, varargin)
    path = fullfile(results_dir, [tag '_solution.mat']);
    prob = define_chemotherapy_problem(tag, varargin{:});
    if exist(path, 'file')
        d = load(path); sol = d.solution; sumy = d.summary;
        return;
    end
    prob.grid.N = N_grid;
    prob.options.enable_pmp = false;
    fprintf('Solving "%s"...\n', tag);
    [sol, sumy] = solveOCP(prob);
    solution = sol;
    summary = sumy;
    save(path, 'solution', 'summary');
end

%==========================================================================
function plot_case(solution, problem, figures_dir)
    tv = solution.time_vec; X = solution.states; U = solution.controls;
    p = problem.param; T = solution.T_final; lw = 2;
    meta = problem.meta;

    fig = figure('Name', meta.scenario, 'Position', [100, 100, 1200, 800]);
    sgtitle(sprintf('%s  (T=%.1f, J=%.1f)', problem.name, T, solution.cost), 'FontWeight', 'bold');

    subplot(2,2,1); hold on;
    plot(tv, U(1,:), 'b-', 'LineWidth', lw);
    if size(U,1) >= 2, plot(tv, U(2,:), 'r--', 'LineWidth', lw); end
    hold off; title('Controls'); xlabel('Time'); grid on; xlim([0 T]);

    subplot(2,2,2);
    plot(tv, p.e' * X, 'k-', 'LineWidth', lw);
    title('Average tumor'); xlabel('Time'); grid on; xlim([0 T]);

    subplot(2,2,3); hold on;
    plot(p.x, X(:,1), 'k--', 'LineWidth', lw);
    plot(p.x, X(:,end), 'm-', 'LineWidth', lw);
    hold off; title('Trait distribution'); xlabel('x'); grid on;

    subplot(2,2,4); hold on;
    plot(p.x, diag(p.Phi1), 'b-', 'LineWidth', lw);
    if any(diag(p.Phi2)), plot(p.x, diag(p.Phi2), 'r--', 'LineWidth', lw); end
    hold off; title('Drug efficacy'); xlabel('x'); grid on;

    saveas(fig, fullfile(figures_dir, ['case_' meta.scenario '.png']));
end

function fig4(results_dir, figures_dir, N_grid)
    [sol, ~, prob] = load_or_solve('exp_gaussian', results_dir, N_grid);
    plot_case(sol, prob, figures_dir);
    saveas(gcf, fullfile(figures_dir, 'figure_4.png'));

    fig = figure('Position', [150, 150, 700, 500]); hold on;
    for s = [1, 2, 3]
        pr = define_chemotherapy_problem('exp_gaussian', 'phi_scale', s);
        pr.grid.N = N_grid; pr.options.enable_pmp = false;
        pr.options.qpso.swarmSize = 80; pr.options.qpso.maxIter = 120;
        [so, ~] = solveOCP(pr);
        plot(pr.param.x, so.states(:,end), 'LineWidth', 2, 'DisplayName', sprintf('\\Phi=%d', s));
    end
    hold off; xlabel('Trait'); ylabel('Final population'); legend; grid on;
    saveas(fig, fullfile(figures_dir, 'figure_4_phi_sensitivity.png'));
end

function fig5(results_dir, figures_dir)
    tags = {'cosine_gaussian', 'cosine_sine', 'exp_gaussian'};
    labs = {'Cosine-Gaussian', 'Cosine-Sine', 'Exp-Gaussian'};
    clr = [0.0,0.45,0.74; 0.85,0.33,0.10; 0.93,0.69,0.13];
    fig = figure('Position', [120, 120, 800, 500]); hold on;
    for k = 1:3
        [so, ~, pr] = load_or_solve(tags{k}, results_dir, 100);
        plot(pr.param.x, so.states(:,end), 'Color', clr(k,:), 'LineWidth', 2.5, 'DisplayName', labs{k});
    end
    hold off; xlabel('Trait'); ylabel('Final population'); title('Fig. 5'); legend; grid on;
    saveas(fig, fullfile(figures_dir, 'figure_5.png'));
end

function fig6(results_dir, figures_dir)
    tags = {'cosine_gaussian', 'cosine_sine', 'exp_gaussian'};
    tit = {'Cosine-Gaussian', 'Cosine-Sine', 'Exponential-Gaussian'};
    fig = figure('Position', [100, 100, 900, 900]);
    for k = 1:3
        [so, ~, pr] = load_or_solve(tags{k}, results_dir, 100);
        subplot(3,1,k);
        plot(so.time_vec, pr.param.e' * so.states, 'k-', 'LineWidth', 2);
        title(tit{k}); xlabel('Time'); ylabel('Avg. tumor'); grid on;
    end
    sgtitle('Fig. 6 — population evolution', 'FontWeight', 'bold');
    saveas(fig, fullfile(figures_dir, 'figure_6.png'));
end
