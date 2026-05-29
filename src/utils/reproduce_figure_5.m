function reproduce_figure_5(save_dir)
%REPRODUCE_FIGURE_5 Fig. 5 — final tumor distributions for all three strategies.

    if nargin < 1
        repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
        save_dir = fullfile(repo_root, 'results', 'figures');
    end
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    results_dir = fileparts(save_dir);
    scenarios = {'cosine_gaussian', 'cosine_sine', 'exp_gaussian'};
    labels = {'Cosine-Gaussian', 'Cosine-Sine', 'Exponential-Gaussian'};
    colors = [0.0, 0.45, 0.74; 0.85, 0.33, 0.10; 0.93, 0.69, 0.13];

    fig = figure('Name', 'Fig. 5 Final tumor distributions', 'Position', [120, 120, 800, 500]);
    hold on;
    for k = 1:numel(scenarios)
        [sol, ~, prob] = solve_scenario_cached(scenarios{k}, results_dir);
        plot(prob.param.x, sol.states(:, end), 'Color', colors(k, :), 'LineWidth', 2.5, ...
            'DisplayName', labels{k});
    end
    hold off;
    xlabel('Trait value'); ylabel('Final tumor population');
    title('Final tumor distributions for all therapy strategies');
    grid on; legend('Location', 'best');
    xlim([0, 1]);
    saveas(fig, fullfile(save_dir, 'figure_5_final_distributions.png'));
    fprintf('Saved Fig. 5 to %s\n', fullfile(save_dir, 'figure_5_final_distributions.png'));
end
