function problem = define_chemotherapy_problem(scenario, varargin)
%DEFINE_CHEMOTHERAPY_PROBLEM  Tumor-trait chemotherapy OCP (QPSO-SQP paper).
%
%   problem = define_chemotherapy_problem('single')
%   problem = define_chemotherapy_problem('cosine_gaussian')
%   problem = define_chemotherapy_problem('cosine_sine')
%   problem = define_chemotherapy_problem('exp_gaussian')
%
%   Optional name-value pairs:
%     'medication'  - 'single' | 'double' (default: inferred from scenario)
%     'mutation'    - true | false (default: false for 'single', true otherwise)
%     'phi_scale'   - scalar multiplier on phi1/phi2 (default: 1; see Fig. 4 sensitivity)
%
%   Paper: IEEE CDC 2025 — doi:10.1109/CDC57313.2025.11312610
%   Dynamics (Michaelis-Menten): Section II; numerical cases: Section V.
%
%   Legacy numeric IDs (define_chemo_trait_problem compatibility):
%     1 -> single,  2 -> cosine_sine,  3 -> exp_gaussian

    p = inputParser;
    p.addRequired('scenario', @(s) ischar(s) || isstring(s) || isnumeric(s));
    p.addParameter('medication', '', @(s) isempty(s) || ismember(lower(char(s)), {'single', 'double'}));
    p.addParameter('mutation', [], @(x) isempty(x) || islogical(x));
    p.addParameter('phi_scale', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
    p.parse(scenario, varargin{:});
    opts = p.Results;

    [scenario_name, legacy_id] = normalize_scenario(opts.scenario);
    spec = scenario_spec(scenario_name);

    if isempty(opts.mutation)
        use_mutation = spec.default_mutation;
    else
        use_mutation = opts.mutation;
    end

    if isempty(opts.medication)
        if spec.default_nu == 1
            med = 'single';
        else
            med = 'double';
        end
    else
        med = lower(char(opts.medication));
    end

  switch med
        case 'single'
            nu = 1;
        case 'double'
            nu = 2;
        otherwise
            error('medication must be ''single'' or ''double''.');
    end

    m = 21;
    fprintf('Loading: Chemotherapy OCP [%s, %s drug(s), mutation=%d]\n', ...
        spec.label, med, use_mutation);

    problem.name = sprintf('Chemotherapy (%s, %s)', spec.label, med);
    problem.nx = m;
    problem.nu = nu;

    problem.time.mode    = 'fixed';
    problem.time.T_fixed = 10.0;
    problem.grid.N       = 100;

    % --- Model parameters ---
    mp.m = m;
    mp.x = linspace(0, 1, m)';
    mp.e = ones(m, 1) / m;
    mp.r = 2 ./ (1 + 3 * mp.x.^4);
    mp.M = diag(0.5 * ones(m, 1));
    mp.sigma = 5;

    if use_mutation
        mp.theta = 0.2;
        mp.interaction_value = 0.2;
    else
        mp.theta = 0.0;
        mp.interaction_value = 0.0;
    end

    mp.P = build_P(mp.x, mp.sigma, m);
    mp.R = build_R(mp.r, mp.theta, mp.P, m);

    [phi1_vec, phi2_vec] = spec.phi_functions(mp.x);
    phi1_vec = opts.phi_scale * phi1_vec;
    phi2_vec = opts.phi_scale * phi2_vec;
    mp.Phi1  = diag(phi1_vec);
    mp.Phi2  = diag(phi2_vec);
    mp.Phi12 = mp.interaction_value * mp.Phi1;
    problem.param = mp;

    % --- Cost weights (Section V.A) ---
    w.alpha_vec     = 5   * ones(m, 1);
    w.beta_vec      = 400 * ones(m, 1);
    w.gamma1        = 4000;
    w.gamma2        = (nu > 1) * 4000;
    w.lambda_smooth = 1e2;
    problem.weights = w;

    problem.model.f            = @(t, x, u, p) ch_f(t, x, u, p, nu);
    problem.model.f_vectorized = @(t_vec, X, U, p) ch_f_vec(t_vec, X, U, p, nu);

    problem.cost.L            = @ch_L;
    problem.cost.L_vectorized = @ch_L_vec;
    problem.cost.Phi          = @ch_Phi;

    problem.bounds.x0   = 10 * ones(m, 1);
    problem.bounds.x_lb = zeros(m, 1);
    problem.bounds.x_ub = 500 * ones(m, 1);
    problem.bounds.u_lb = zeros(nu, 1);
    umax = 3 * ones(nu, 1);
    problem.bounds.u_ub = umax;

    problem.constraints.boundary_eq   = [];
    problem.constraints.boundary_ineq = [];
    problem.constraints.path_ineq     = [];

    problem.options.qpso.swarmSize    = 150;
    problem.options.qpso.maxIter      = 250;
    problem.options.smooth_qpso_guess = (nu == 2);
    problem.options.enable_pmp        = true;
    problem.options.state_init_mode   = 'rollout';

    problem.pmp_derivatives  = @(st, ct, co, p, w, tv) ch_pmp_derivatives(st, ct, co, p, w, tv, nu);
    problem.generate_summary = @ch_summary;

    problem.meta = struct();
    problem.meta.scenario      = scenario_name;
    problem.meta.legacy_id     = legacy_id;
    problem.meta.medication    = med;
    problem.meta.mutation      = use_mutation;
    problem.meta.phi_scale     = opts.phi_scale;
    problem.meta.paper_section = spec.paper_section;
end

%==========================================================================
function [name, legacy_id] = normalize_scenario(scenario)
    if isnumeric(scenario)
        legacy_id = scenario;
        switch scenario
            case 1, name = 'single';
            case 2, name = 'cosine_sine';
            case 3, name = 'exp_gaussian';
            otherwise
                error('Legacy scenario id must be 1, 2, or 3.');
        end
        return;
    end
    legacy_id = [];
    name = lower(strtrim(char(scenario)));
    switch name
        case {'1', 'fig1', 'pmp'},           name = 'single';
        case '2',                            name = 'cosine_sine';
        case '3',                            name = 'exp_gaussian';
        case {'cosine_gauss', 'cg'},         name = 'cosine_gaussian';
        case 'cs',                           name = 'cosine_sine';
        case 'eg',                           name = 'exp_gaussian';
    end
    valid = {'single', 'cosine_gaussian', 'cosine_sine', 'exp_gaussian'};
    if ~ismember(name, valid)
        error('Unknown scenario "%s". Use: %s', name, strjoin(valid, ', '));
    end
end

%==========================================================================
function spec = scenario_spec(name)
    spec = struct('name', name, 'default_nu', 2, 'default_mutation', true, ...
        'paper_section', '', 'label', '', 'phi_functions', []);
    switch name
        case 'single'
            spec.default_nu = 1;
            spec.default_mutation = false;
            spec.paper_section = 'Fig. 1 (single drug, no mutation)';
            spec.label = 'Single-drug PMP benchmark';
            spec.phi_functions = @phi_single_pmp;
        case 'cosine_gaussian'
            spec.paper_section = 'Sec. V.B.1 Cosine-Gaussian';
            spec.label = 'Cosine-Gaussian combination';
            spec.phi_functions = @phi_cosine_gaussian;
        case 'cosine_sine'
            spec.paper_section = 'Sec. V.B.2 Cosine-Sine';
            spec.label = 'Cosine-Sine combination';
            spec.phi_functions = @phi_cosine_sine;
        case 'exp_gaussian'
            spec.paper_section = 'Sec. V.B.3 Exponential-Gaussian';
            spec.label = 'Exponential-Gaussian combination';
            spec.phi_functions = @phi_exp_gaussian;
    end
end

%==========================================================================
function [phi1, phi2] = phi_single_pmp(x)
    phi1 = -sin(x - 1) + 1.5;
    phi2 = zeros(size(x));
end

function [phi1, phi2] = phi_cosine_gaussian(x)
    phi1 = 1 + cos(0.25 * pi * x).^2;
    phi2 = 1 + exp(-(x - 0.2).^2 / (2 * 0.25^2));
end

function [phi1, phi2] = phi_cosine_sine(x)
    phi1 = 1 + cos(0.25 * pi * x).^2;
    phi2 = 1 + sin(0.25 * pi * x).^2;
end

function [phi1, phi2] = phi_exp_gaussian(x)
    phi1 = 1 + exp(-3 * x);
    phi2 = 1 + exp(-(x - 0.5).^2 / (2 * 0.25^2));
end

%==========================================================================
function dNdt = ch_f(~, N, u, p, nu)
    N = max(N(:), 0);
    u1_MM  = u(1) / (1 + u(1));
    N_bar  = p.e' * N;
    G_N    = log(max(1 + N_bar, 1e-6));
    A = p.R - p.Phi1 * u1_MM - p.M * G_N;
    if nu >= 2
        u2_MM  = u(2) / (1 + u(2));
        u12_MM = u1_MM * u2_MM;
        A = A - p.Phi2 * u2_MM + p.Phi12 * u12_MM;
    end
    dNdt = A * N;
end

function F = ch_f_vec(~, X, U, p, nu)
    N_bar = p.e' * X;
    G_N   = log(max(1 + N_bar, 1e-6));
    u1_MM = U(1, :) ./ (1 + U(1, :));
    F = p.R * X - p.Phi1 * (X .* u1_MM) - p.M * (X .* G_N);
    if nu >= 2
        u2_MM  = U(2, :) ./ (1 + U(2, :));
        u12_MM = u1_MM .* u2_MM;
        F = F - p.Phi2 * (X .* u2_MM) + p.Phi12 * (X .* u12_MM);
    end
end

function val = ch_L(~, x, u, ~, ~, w)
    val = w.beta_vec' * x + w.gamma1 * u(1);
    if numel(u) >= 2
        val = val + w.gamma2 * u(2);
    end
end

function L_vec = ch_L_vec(~, X, U, ~, ~, w)
    L_vec = w.beta_vec' * X + w.gamma1 * U(1, :);
    if size(U, 1) >= 2
        L_vec = L_vec + w.gamma2 * U(2, :);
    end
end

function val = ch_Phi(~, xF, ~, ~, w)
    val = w.alpha_vec' * xF;
end

function [dH_du_all, d2H_du2_all, dH_du_dt_all] = ...
        ch_pmp_derivatives(states, controls, costates, p, w, time_vec, nu)

    N_opt = states;
    lam   = costates;
    u1    = controls(1, :);
    u1_MM = u1 ./ (1 + u1);
    phi1_vec  = diag(p.Phi1);
    phi12_vec = diag(p.Phi12);

    if nu >= 2
        u2    = controls(2, :);
        u2_MM = u2 ./ (1 + u2);
        phi2_vec = diag(p.Phi2);
        term_u1 = (-phi1_vec + phi12_vec .* u2_MM) ./ (1 + u1).^2;
        term_u2 = (-phi2_vec + phi12_vec .* u1_MM) ./ (1 + u2).^2;
        dH_du1 = w.gamma1 + sum(lam .* (term_u1 .* N_opt), 1);
        dH_du2 = w.gamma2 + sum(lam .* (term_u2 .* N_opt), 1);
        t2_u1  = 2 * (phi1_vec - phi12_vec .* u2_MM) ./ (1 + u1).^3;
        t2_u2  = 2 * (phi2_vec - phi12_vec .* u1_MM) ./ (1 + u2).^3;
        d2H_du1 = sum(lam .* (t2_u1 .* N_opt), 1);
        d2H_du2 = sum(lam .* (t2_u2 .* N_opt), 1);
        dt = time_vec(2) - time_vec(1);
        dH_du_all   = [dH_du1;  dH_du2];
        d2H_du2_all = [d2H_du1; d2H_du2];
        dH_du_dt_all = [gradient(dH_du1, dt); gradient(dH_du2, dt)];
    else
        term_u1 = -phi1_vec ./ (1 + u1).^2;
        dH_du1  = w.gamma1 + sum(lam .* (term_u1 .* N_opt), 1);
        t2_u1   = 2 * phi1_vec ./ (1 + u1).^3;
        d2H_du1 = sum(lam .* (t2_u1 .* N_opt), 1);
        dt = time_vec(2) - time_vec(1);
        dH_du_all   = dH_du1;
        d2H_du2_all = d2H_du1;
        dH_du_dt_all = gradient(dH_du1, dt);
    end
end

function ch_summary(solution, summary, problem)
    tv = solution.time_vec;
    X  = solution.states;
    U  = solution.controls;
    w  = problem.weights;
    T  = solution.T_final;
    meta = problem.meta;

    J_mayer = problem.cost.Phi([], X(:, end), T, problem.param, w);
    J_total = solution.cost;
    J_lag   = J_total - J_mayer;
    N_avg0  = problem.param.e' * X(:, 1);
    N_avgF  = problem.param.e' * X(:, end);

    fprintf(repmat('=', 1, 70));
    fprintf('\n   CHEMOTHERAPY OCP SUMMARY — %s\n', meta.scenario);
    fprintf('   %s | %s medication | mutation=%d\n', ...
        meta.paper_section, meta.medication, meta.mutation);
    fprintf(repmat('=', 1, 70));
    fprintf('\n');
    fprintf('  Duration:         %.3f (fixed)\n', T);
    fprintf('  Computation:      %.2f sec\n', summary.total_time);
    fprintf('  Lagrange cost:    %.2f (%.1f%%)\n', J_lag, 100 * J_lag / J_total);
    fprintf('  Mayer cost:       %.2f (%.1f%%)\n', J_mayer, 100 * J_mayer / J_total);
    fprintf('  Total cost:       %.2f\n', J_total);
    fprintf('  Avg tumor:        %.3f -> %.3f\n', N_avg0, N_avgF);
    fprintf('  Drug 1 dose:      %.3f\n', trapz(tv, U(1, :)));
    if size(U, 1) >= 2
        fprintf('  Drug 2 dose:      %.3f\n', trapz(tv, U(2, :)));
    end
    fprintf(repmat('=', 1, 70));
    fprintf('\n\n');
end

function P = build_P(x, sigma, m)
    P = zeros(m);
    for i = 1:m
        for j = 1:m
            P(i, j) = exp(-0.5 * ((x(i) - x(j)) / sigma)^2);
        end
        P(i, :) = P(i, :) / sum(P(i, :));
    end
end

function R = build_R(r, theta, P, m)
    R = zeros(m);
    for i = 1:m
        for j = 1:m
            if i == j
                R(i, j) = r(i) * (1 - theta);
            else
                R(i, j) = theta * r(j) * P(i, j);
            end
        end
    end
end
