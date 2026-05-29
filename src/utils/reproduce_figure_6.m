function reproduce_figure_6(save_dir)
%REPRODUCE_FIGURE_6 Fig. 6 — tumor population evolution under each regimen.

    if nargin < 1
        repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
        save_dir = fullfile(repo_root, 'results', 'figures');
    end
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    results_dir = fileparts(save_dir);
    scenarios = {'cosine_gaussian', 'cosine_sine', 'exp_gaussian'};
    titles = {'Cosine-Gaussian', 'Cosine-Sine', 'Exponential-Gaussian'};
    lw = 2.0;

    fig = figure('Name', 'Fig. 6 Tumor evolution', 'Position', [100, 100, 900, 900]);
    for k = 1:numel(scenarios)
        [sol, ~, prob] = solve_scenario_cached(scenarios{k}, results_dir);
        subplot(3, 1, k);
        N_avg = prob.param.e' * sol.states;
        plot(sol.time_vec, N_avg, 'k-', 'LineWidth', lw);
        title(titles{k});
        xlabel('Time (years)'); ylabel('Average tumor population');
        grid on; xlim([0, sol.T_final]);
    end
    sgtitle('Tumor population evolution under different treatment regimens', ...
        'FontSize', 14, 'FontWeight', 'bold');
    saveas(fig, fullfile(save_dir, 'figure_6_population_evolution.png'));
    fprintf('Saved Fig. 6 to %s\n', fullfile(save_dir, 'figure_6_population_evolution.png'));
end
