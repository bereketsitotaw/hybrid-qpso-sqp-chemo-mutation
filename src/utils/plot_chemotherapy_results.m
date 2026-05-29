function plot_chemotherapy_results(solution, summary, problem, save_dir)
%PLOT_CHEMOTHERAPY_RESULTS Paper-style trait-model control and biology plots.

    if nargin < 4
        save_dir = '';
    end

    lw = 2.0;
    tv = solution.time_vec;
    X  = solution.states;
    U  = solution.controls;
    p  = problem.param;
    T  = solution.T_final;
    meta = problem.meta;

    fig = figure('Name', ['Chemotherapy: ' meta.scenario], ...
        'Position', [100, 100, 1200, 800]);
    sgtitle(sprintf('%s  (T=%.1f yr, J=%.1f)', problem.name, T, solution.cost), ...
        'FontSize', 14, 'FontWeight', 'bold');

    subplot(2, 2, 1);
    hold on;
    plot(tv, U(1, :), 'b-', 'LineWidth', lw, 'DisplayName', 'Drug 1 (u_1)');
    if size(U, 1) >= 2
        plot(tv, U(2, :), 'r--', 'LineWidth', lw, 'DisplayName', 'Drug 2 (u_2)');
    end
    hold off;
    title('Optimal Control Inputs');
    xlabel('Time (years)'); ylabel('Dosage');
    grid on; legend('Location', 'best');
    xlim([0, T]);

    subplot(2, 2, 2);
    N_avg = p.e' * X;
    plot(tv, N_avg, 'k-', 'LineWidth', lw);
    title('Average Tumor Population');
    xlabel('Time (years)'); ylabel('Average population');
    grid on; xlim([0, T]);

    subplot(2, 2, 3);
    hold on;
    plot(p.x, X(:, 1), 'k--', 'LineWidth', lw, 'DisplayName', 'Initial');
    plot(p.x, X(:, end), 'm-', 'LineWidth', lw, 'DisplayName', 'Final');
    hold off;
    title('Tumor Distribution Across Traits');
    xlabel('Resistance trait x'); ylabel('Population');
    grid on; legend('Location', 'best');

    subplot(2, 2, 4);
    phi1 = diag(p.Phi1);
    phi2 = diag(p.Phi2);
    hold on;
    plot(p.x, phi1, 'b-', 'LineWidth', lw, 'DisplayName', '\phi_1(x)');
    if any(phi2 ~= 0)
        plot(p.x, phi2, 'r--', 'LineWidth', lw, 'DisplayName', '\phi_2(x)');
    end
    hold off;
    title('Drug Efficacy Parameters');
    xlabel('Trait x'); ylabel('Parameter value');
    grid on; legend('Location', 'best');

    if ~isempty(save_dir)
        fname = fullfile(save_dir, ['case_' meta.scenario]);
        saveas(fig, [fname '.png']);
    end

    if isfield(solution, 'pmp') && ~isempty(fieldnames(solution.pmp))
        plot_pmp_diagnostics(solution, save_dir, meta.scenario, lw);
    end
end

%--------------------------------------------------------------------------
function plot_pmp_diagnostics(solution, save_dir, scenario, lw)
    pmp = solution.pmp;
    tv  = solution.time_vec;
    nu  = size(solution.controls, 1);

    fig = figure('Name', ['PMP: ' scenario], 'Position', [200, 200, 1100, 700]);
    sgtitle('PMP Optimality Verification', 'FontSize', 14, 'FontWeight', 'bold');

    subplot(2, 2, 1);
    if isfield(pmp, 'hamiltonian') && ~isempty(pmp.hamiltonian)
        plot(tv, pmp.hamiltonian, 'k-', 'LineWidth', lw);
        title('Hamiltonian H(t)'); xlabel('Time'); grid on;
    end

    if isfield(pmp, 'dH_du') && ~isempty(pmp.dH_du)
        subplot(2, 2, 2);
        for i = 1:nu
            plot(tv, pmp.dH_du(i, :), 'LineWidth', lw, ...
                'DisplayName', sprintf('dH/du_%d', i));
            hold on;
        end
        title('First-order optimality (dH/du)'); xlabel('Time');
        grid on; legend show;

        subplot(2, 2, 3);
        for i = 1:nu
            plot(tv, pmp.d2H_du2(i, :), 'LineWidth', lw, ...
                'DisplayName', sprintf('d^2H/du_%d^2', i));
            hold on;
        end
        title('Second-order (Legendre-Clebsch)'); xlabel('Time');
        grid on; legend show;

        subplot(2, 2, 4);
        for i = 1:nu
            plot(tv, pmp.dH_du_dt(i, :), 'LineWidth', lw, ...
                'DisplayName', sprintf('d/dt(dH/du_%d)', i));
            hold on;
        end
        title('Singular-arc indicator d/dt(dH/du)'); xlabel('Time');
        grid on; legend show;
    end

    if ~isempty(save_dir)
        saveas(fig, fullfile(save_dir, ['pmp_' scenario '.png']));
    end
end
