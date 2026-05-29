function reproduce_figure_4(save_dir)
%REPRODUCE_FIGURE_4 Fig. 4 — Exponential-Gaussian therapy dynamics + phi sensitivity.
%
%   reproduce_figure_4()
%
%   Top: optimal controls; bottom-left: average tumor; bottom-right: phi(x).
%   Also plots final trait distribution for phi_scale = 1, 2, 3 (paper text).

    if nargin < 1
        repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
        save_dir = fullfile(repo_root, 'results', 'figures');
    end
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    [solution, ~, problem] = solve_scenario_cached('exp_gaussian', ...
        fileparts(save_dir));
    plot_chemotherapy_results(solution, struct(), problem, save_dir);
    saveas(gcf, fullfile(save_dir, 'figure_4_exp_gaussian.png'));

    % Phi-scale sensitivity (Phi1 = Phi2 = 1, 2, 3)
    scales = [1, 2, 3];
    colors = lines(numel(scales));
    fig = figure('Name', 'Fig. 4 sensitivity (phi scale)', 'Position', [150, 150, 700, 500]);
    hold on;
    for k = 1:numel(scales)
        s = scales(k);
        prob = define_chemotherapy_problem('exp_gaussian', 'phi_scale', s);
        prob.grid.N = reproducibility_settings().N;
        prob.options.enable_pmp = false;
        prob.options.qpso.swarmSize = 80;
        prob.options.qpso.maxIter = 120;
        [sol, ~] = solveOCP(prob);
        plot(prob.param.x, sol.states(:, end), 'Color', colors(k, :), 'LineWidth', 2, ...
            'DisplayName', sprintf('\\Phi_1=\\Phi_2=%d', s));
    end
    hold off;
    xlabel('Trait value'); ylabel('Final tumor population');
    title('Final tumor distribution vs. efficacy scale (Exp-Gaussian)');
    grid on; legend('Location', 'best');
    saveas(fig, fullfile(save_dir, 'figure_4_phi_sensitivity.png'));
    fprintf('Saved Fig. 4 artifacts to %s\n', save_dir);
end
