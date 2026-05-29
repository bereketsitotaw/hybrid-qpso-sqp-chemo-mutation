function [solution, summary] = solveOCP(problem)
%SOLVEOCP  Solve a single-phase nonlinear optimal control problem.
%
%   [solution, summary] = solveOCP(problem)
%
%   Method: Hybrid QPSO global search  +  fmincon (SQP) local refinement
%   Transcription: Direct trapezoidal collocation
%
%   The problem struct API is documented in the define_*_problem files.

    fprintf('========================================\n');
    fprintf('  solveOCP: %s\n', problem.name);
    fprintf('========================================\n');

    % --- 1. Validate and normalize ---
    problem = validate_problem(problem);
    N  = problem.grid.N;
    nx = problem.nx;
    nu = problem.nu;
    debug = problem.options.debug_feasibility;

    % --- 2. QPSO global search ---
    fprintf('Phase 1: QPSO global search...\n');
    [u_qpso, x_qpso, qpso_history, qpso_time, T_qpso] = ocp_qpso(problem);

    % --- 3. Optionally smooth QPSO controls ---
    if problem.options.smooth_qpso_guess
        fprintf('   -> Smoothing QPSO controls...\n');
        for i = 1:nu
            u_qpso(i, :) = smoothdata(u_qpso(i, :), 'gaussian', 10);
        end
        u_qpso = max(bsxfun(@max, u_qpso, problem.bounds.u_lb), problem.bounds.u_lb);
        u_qpso = min(bsxfun(@min, u_qpso, problem.bounds.u_ub), problem.bounds.u_ub);
        
        fprintf('   -> Re-simulating trajectory for dynamic consistency...\n');
        t_vec_qpso = linspace(0, T_qpso, N + 1);
        [x_qpso_new, resim_ok] = ocp_forward_sim(u_qpso, T_qpso, t_vec_qpso, problem, false);
        if resim_ok
            x_qpso = x_qpso_new;
        else
            warning('solveOCP: re-simulation after smoothing failed. Falling back to unsmoothed states.');
        end
    else
        fprintf('   -> Smoothing disabled for this problem.\n');
    end

    % --- 4. Build NLP initial guess ---
    T0 = T_qpso;
    
    % State initialization mode
    switch problem.options.state_init_mode
        case 'rollout'
            % Use QPSO forward-simulated states (default)
            x_init = x_qpso;
            fprintf('   -> State init mode: rollout (QPSO forward sim)\n');
        case 'linear'
            % Linear interpolation from x0 to QPSO final state
            x0_col = problem.bounds.x0_lb(:);
            xf_col = x_qpso(:, end);
            x_init = x0_col + (xf_col - x0_col) * linspace(0, 1, N + 1);
            fprintf('   -> State init mode: linear (x0 -> xF_qpso)\n');
        otherwise
            x_init = x_qpso;
    end
    
    % --- Dynamics-aware defect scaling ---
    % Override the static 1/|x0| scaling with per-state magnitude estimates
    % from the actual QPSO trajectory. This balances Jacobian rows for SQP.
    t_vec_init = linspace(0, T0, N + 1);
    F_init = ocp_dynamics(t_vec_init, x_init, u_qpso, problem);
    defect_mag = max(abs(F_init), [], 2) * (T0 / N);  % [nx,1] typical defect size
    defect_mag = max(defect_mag, 1e-6);                % floor to avoid div-by-zero
    problem.scaling.defect_scale = 1.0 ./ defect_mag;  % override static guess
    fprintf('   -> Dynamics-aware defect scaling: min=%.2e  max=%.2e\n', ...
            min(problem.scaling.defect_scale), max(problem.scaling.defect_scale));

    z0 = pack_z(x_init, u_qpso, T0, problem);

    % --- 5. Build NLP bounds & Validate initial point ---
    [lb, ub] = get_nlp_bounds(problem);

    % NaN/Inf check
    if any(isnan(z0)) || any(isinf(z0))
        error('solveOCP: Initial guess z0 from QPSO contains non-finite values (NaN/Inf).');
    end

    % Clip z0 to NLP bounds
    finite_lb = isfinite(lb);
    finite_ub = isfinite(ub);
    z0(finite_lb) = max(z0(finite_lb), lb(finite_lb));
    z0(finite_ub) = min(z0(finite_ub), ub(finite_ub));

    assert(all(isfinite(z0)), 'solveOCP: z0 contains non-finite values after clipping.');
    
    lb_viol = 0; ub_viol = 0;
    if any(finite_lb), lb_viol = max(lb(finite_lb) - z0(finite_lb)); end
    if any(finite_ub), ub_viol = max(z0(finite_ub) - ub(finite_ub)); end
    
    fprintf('   Pre-fmincon z0 bound check | max LB viol = %.3e | max UB viol = %.3e\n', ...
        max(lb_viol,0), max(ub_viol,0));

    % Objective evaluation & dynamic scaling
    try
        J0_raw = nlp_objective_raw(z0, problem);
        if ~isfinite(J0_raw)
            error('solveOCP: Objective function is non-finite at initial point (J0=%.4e).', J0_raw);
        end
        
        % Set scaling factor such that initial objective is ~100
        if J0_raw > 1e2
            problem.scaling.obj_scale = 100 / J0_raw;
        else
            problem.scaling.obj_scale = 1.0;
        end
        
        J0 = J0_raw * problem.scaling.obj_scale;
        fprintf('   -> Initial objective (raw): J0 = %.4f\n', J0_raw);
        fprintf('   -> Objective scaling: %.2e (J0_scaled = %.2f)\n', ...
                problem.scaling.obj_scale, J0);
    catch ME
        error('solveOCP: Error evaluating objective at initial point: %s', ME.message);
    end

    % Feasibility diagnostics at z0
    print_feasibility_breakdown(z0, problem, 'z0 (pre-fmincon)');

    % --- 6. fmincon local refinement ---
    fprintf('Phase 2: fmincon (SQP) refinement...\n');
    tic;

    nlp_opts = optimoptions('fmincon', 'Algorithm', 'sqp', ...
        'Display', 'iter', ...
        'MaxIterations', 5000, ...
        'MaxFunctionEvaluations', 1e6, ...
        'ConstraintTolerance', 1e-6, ...
        'OptimalityTolerance', 1e-6, ...
        'FiniteDifferenceType', 'central', ...
        'UseParallel', true, ...
        'ScaleProblem', true);

    % Override with user NLP options
    if isfield(problem.options, 'nlp')
        fns = fieldnames(problem.options.nlp);
        for i = 1:length(fns)
            nlp_opts.(fns{i}) = problem.options.nlp.(fns{i});
        end
    end

    [z_opt, J_opt_scaled, exitflag, output, lambda_struct] = fmincon( ...
        @(z) nlp_objective(z, problem), ...
        z0, [], [], [], [], lb, ub, ...
        @(z) nlp_constraints(z, problem), ...
        nlp_opts);

    % Unscale the final objective result
    J_opt = J_opt_scaled / problem.scaling.obj_scale;

    nlp_time = toc;
    fprintf('   -> fmincon done in %.2f s.  J=%.4f  exitflag=%d\n', ...
            nlp_time, J_opt, exitflag);

    % Feasibility diagnostics at z_opt
    print_feasibility_breakdown(z_opt, problem, 'z_opt (post-fmincon)');

    % --- 7. Package solution ---
    [X_opt, U_opt, T_opt] = unpack_z(z_opt, problem);

    solution.problem_name = problem.name;
    solution.time_vec     = linspace(0, T_opt, N + 1);
    solution.states       = X_opt;
    solution.controls     = U_opt;
    solution.T_final      = T_opt;
    solution.cost         = J_opt;
    solution.exitflag     = exitflag;

    summary.problem_name  = problem.name;
    summary.qpso_time     = qpso_time;
    summary.nlp_time      = nlp_time;
    summary.total_time    = qpso_time + nlp_time;
    summary.exitflag      = exitflag;
    summary.nlp_output    = output;
    summary.lambda        = lambda_struct;
    summary.qpso_history  = qpso_history;

    % --- 8. Optional PMP analysis (guard on feasibility) ---
    if problem.options.enable_pmp
        % Check feasibility before running expensive PMP
        [~, ceq_final] = nlp_constraints(z_opt, problem);
        max_ceq_viol = max(abs(ceq_final));
        
        if exitflag > 0 || max_ceq_viol < 1e-3
            if ~debug
                fprintf('Post-processing: PMP analysis...\n');
                solution.pmp = ocp_pmp(solution, summary, problem);
            else
                fprintf('   PMP skipped (debug_feasibility mode).\n');
            end
        else
            fprintf('   Skipping PMP post-processing because NLP did not converge\n');
            fprintf('   to a feasible local minimum (exitflag=%d, max|ceq|=%.2e).\n', ...
                    exitflag, max_ceq_viol);
        end
    end

    fprintf('========================================\n');
    fprintf('  Solver complete.\n');
    fprintf('========================================\n');
end

%==========================================================================
%                         LOCAL FUNCTIONS
%==========================================================================

%--- Validate and normalize the problem struct ---
function prob = validate_problem(prob)
    assert(isfield(prob, 'name'),       'Missing problem.name');
    assert(isfield(prob, 'nx'),         'Missing problem.nx');
    assert(isfield(prob, 'nu'),         'Missing problem.nu');
    assert(isfield(prob, 'time'),       'Missing problem.time');
    assert(isfield(prob, 'grid'),       'Missing problem.grid');
    assert(isfield(prob, 'model'),      'Missing problem.model');
    assert(isfield(prob.model, 'f'),    'Missing problem.model.f');
    assert(isfield(prob, 'cost'),       'Missing problem.cost');
    assert(isfield(prob.cost, 'L'),     'Missing problem.cost.L');
    assert(isfield(prob, 'bounds'),     'Missing problem.bounds');

    % Normalize x0 shorthand
    if isfield(prob.bounds, 'x0') && ~isfield(prob.bounds, 'x0_lb')
        prob.bounds.x0_lb = prob.bounds.x0(:);
        prob.bounds.x0_ub = prob.bounds.x0(:);
    end

    % Default state bounds
    if ~isfield(prob.bounds, 'x_lb') || isempty(prob.bounds.x_lb)
        prob.bounds.x_lb = -inf(prob.nx, 1);
    end
    if ~isfield(prob.bounds, 'x_ub') || isempty(prob.bounds.x_ub)
        prob.bounds.x_ub = inf(prob.nx, 1);
    end

    % Default terminal cost
    if ~isfield(prob.cost, 'Phi'), prob.cost.Phi = []; end

    % Default constraints
    if ~isfield(prob, 'constraints'), prob.constraints = struct(); end
    if ~isfield(prob.constraints, 'boundary_eq'),   prob.constraints.boundary_eq   = []; end
    if ~isfield(prob.constraints, 'boundary_ineq'), prob.constraints.boundary_ineq = []; end
    if ~isfield(prob.constraints, 'path_ineq'),     prob.constraints.path_ineq     = []; end

    % Sanitize search bounds for QPSO (must be finite)
    u_lb = prob.bounds.u_lb;
    u_ub = prob.bounds.u_ub;
    u_lb_qpso = u_lb;
    u_ub_qpso = u_ub;

    for i = 1:prob.nu
        if ~isfinite(u_lb_qpso(i))
            u_lb_qpso(i) = -5.0;
        end
        if ~isfinite(u_ub_qpso(i))
            u_ub_qpso(i) = 5.0;
        end
        if u_ub_qpso(i) <= u_lb_qpso(i)
            u_ub_qpso(i) = u_lb_qpso(i) + 1.0;
        end
    end
    prob.bounds.u_lb_qpso = u_lb_qpso;
    prob.bounds.u_ub_qpso = u_ub_qpso;

    % Default options
    if ~isfield(prob, 'options'), prob.options = struct(); end
    if ~isfield(prob.options, 'smooth_qpso_guess'), prob.options.smooth_qpso_guess = true; end
    if ~isfield(prob.options, 'enable_pmp'),        prob.options.enable_pmp = false; end
    if ~isfield(prob.options, 'debug_feasibility'), prob.options.debug_feasibility = false; end
    if ~isfield(prob.options, 'state_init_mode'),   prob.options.state_init_mode = 'rollout'; end
    if ~isfield(prob.options, 'qpso')
        prob.options.qpso = struct('swarmSize', 50, 'maxIter', 100);
    end

    % Default param and weights
    if ~isfield(prob, 'param'),   prob.param   = struct(); end
    if ~isfield(prob, 'weights'), prob.weights = struct(); end

    % --- Auto-scaling ---
    % Compute scaling factors from initial conditions if not user-provided
    if ~isfield(prob, 'scaling'), prob.scaling = struct(); end
    
    x0_mag = abs(prob.bounds.x0_lb(:));
    
    if ~isfield(prob.scaling, 'state_scale')
        prob.scaling.state_scale = max(x0_mag, 1.0);  % per-state, at least 1
    end
    if ~isfield(prob.scaling, 'defect_scale')
        % Scale defects so each row is O(1) regardless of state magnitude
        prob.scaling.defect_scale = 1.0 ./ prob.scaling.state_scale;
    end
    if ~isfield(prob.scaling, 'boundary_eq_scale')
        % Scale boundary eq constraints by magnitude of related states
        prob.scaling.boundary_eq_scale = 1.0 ./ max(x0_mag, 1.0);
    end
    if ~isfield(prob.scaling, 'obj_scale')
        % Generic objective scaling: target magnitude ~ 10-100
        % We'll estimate this once later or use a safe default
        prob.scaling.obj_scale = 1.0;
    end
end

%--- Pack states, controls, (time) into decision vector ---
function z = pack_z(X, U, T, problem)
    z = [reshape(X, [], 1); reshape(U, [], 1)];
    if strcmp(problem.time.mode, 'free')
        z = [z; T];
    end
end

%--- Unpack decision vector ---
function [X, U, T] = unpack_z(z, problem)
    nx = problem.nx;
    nu = problem.nu;
    N  = problem.grid.N;

    n_x = nx * (N + 1);
    n_u = nu * (N + 1);

    X = reshape(z(1:n_x),            nx, N + 1);
    U = reshape(z(n_x+1 : n_x+n_u), nu, N + 1);

    if strcmp(problem.time.mode, 'free')
        T = z(end);
    else
        T = problem.time.T_fixed;
    end
end

%--- Build NLP bound vectors ---
function [lb, ub] = get_nlp_bounds(problem)
    nx = problem.nx;
    nu = problem.nu;
    N  = problem.grid.N;

    % State bounds [nx, N+1]
    lb_x = repmat(problem.bounds.x_lb(:), 1, N + 1);
    ub_x = repmat(problem.bounds.x_ub(:), 1, N + 1);

    % Pin initial state
    lb_x(:, 1) = problem.bounds.x0_lb;
    ub_x(:, 1) = problem.bounds.x0_ub;

    % Control bounds [nu, N+1]
    lb_u = repmat(problem.bounds.u_lb(:), 1, N + 1);
    ub_u = repmat(problem.bounds.u_ub(:), 1, N + 1);

    lb = [reshape(lb_x, [], 1); reshape(lb_u, [], 1)];
    ub = [reshape(ub_x, [], 1); reshape(ub_u, [], 1)];

    % Free final time
    if strcmp(problem.time.mode, 'free')
        lb = [lb; problem.time.T_bounds(1)];
        ub = [ub; problem.time.T_bounds(2)];
    end
end

%--- RAW NLP objective (for scaling estimation) ---
function J = nlp_objective_raw(z, problem)
    [X, U, T] = unpack_z(z, problem);
    t_vec = linspace(0, T, problem.grid.N + 1);
    J = ocp_cost(t_vec, X, U, T, problem);
end

%--- NLP objective (scaled for solver) ---
function J_scaled = nlp_objective(z, problem)
    J = nlp_objective_raw(z, problem);
    J_scaled = J * problem.scaling.obj_scale;
end

%--- NLP constraints (with scaling) ---
function [c, ceq] = nlp_constraints(z, problem)
    [X, U, T] = unpack_z(z, problem);
    N     = problem.grid.N;
    nx    = problem.nx;
    h     = T / N;
    t_vec = linspace(0, T, N + 1);

    % Equality: trapezoidal collocation defects
    F = ocp_dynamics(t_vec, X, U, problem);
    defects = X(:, 2:end) - X(:, 1:end-1) - (h/2) * (F(:, 1:end-1) + F(:, 2:end));
    
    % Scale defects: multiply each state-row by its defect_scale factor
    % This makes all defect constraint rows O(1) regardless of state magnitude
    ds = problem.scaling.defect_scale(:);  % [nx, 1]
    defects_scaled = bsxfun(@times, defects, ds);
    
    ceq = reshape(defects_scaled, [], 1);

    % Equality: boundary constraints (scaled)
    if ~isempty(problem.constraints.boundary_eq)
        ceq_bc = problem.constraints.boundary_eq(X(:,1), X(:,end), T, problem.param);
        ceq_bc = ceq_bc(:);
        
        % Scale boundary eq: use boundary_eq_scale truncated/padded to match length
        bes = problem.scaling.boundary_eq_scale(:);
        n_bc = length(ceq_bc);
        if length(bes) >= n_bc
            bc_scale = bes(1:n_bc);
        else
            bc_scale = [bes; ones(n_bc - length(bes), 1)];
        end
        ceq = [ceq; ceq_bc .* bc_scale];
    end

    % Inequality: boundary constraints
    c = [];
    if ~isempty(problem.constraints.boundary_ineq)
        c_bc = problem.constraints.boundary_ineq(X(:,1), X(:,end), T, problem.param);
        c = [c; c_bc(:)];
    end

    % Inequality: path constraints at each node
    if ~isempty(problem.constraints.path_ineq)
        c_path = [];
        for k = 1:(N + 1)
            c_k = problem.constraints.path_ineq(t_vec(k), X(:,k), U(:,k), problem.param);
            c_path = [c_path; c_k(:)]; %#ok<AGROW>
        end
        c = [c; c_path];
    end
end

%--- Per-family feasibility diagnostics ---
function print_feasibility_breakdown(z, problem, label)
    [X, U, T] = unpack_z(z, problem);
    N     = problem.grid.N;
    nx    = problem.nx;
    h     = T / N;
    t_vec = linspace(0, T, N + 1);

    fprintf('\n   --- Feasibility Breakdown: %s ---\n', label);

    % 1. Collocation defect residuals (unscaled, true physics)
    F = ocp_dynamics(t_vec, X, U, problem);
    defects = X(:, 2:end) - X(:, 1:end-1) - (h/2) * (F(:, 1:end-1) + F(:, 2:end));
    fprintf('   max defect residual      : %.3e\n', max(abs(defects(:))));
    
    % Per-state defect breakdown
    for i = 1:nx
        fprintf('     state %2d defect max    : %.3e\n', i, max(abs(defects(i,:))));
    end

    % 2. Boundary equality
    if ~isempty(problem.constraints.boundary_eq)
        ceq_bc = problem.constraints.boundary_eq(X(:,1), X(:,end), T, problem.param);
        fprintf('   max boundary eq residual : %.3e\n', max(abs(ceq_bc(:))));
    else
        fprintf('   max boundary eq residual : (none)\n');
    end

    % 3. Boundary inequality
    if ~isempty(problem.constraints.boundary_ineq)
        c_bc = problem.constraints.boundary_ineq(X(:,1), X(:,end), T, problem.param);
        fprintf('   max boundary ineq viol   : %.3e\n', max(max(c_bc(:)), 0));
    else
        fprintf('   max boundary ineq viol   : (none)\n');
    end

    % 4. Path inequality
    if ~isempty(problem.constraints.path_ineq)
        max_path = 0;
        for k = 1:(N+1)
            c_k = problem.constraints.path_ineq(t_vec(k), X(:,k), U(:,k), problem.param);
            max_path = max(max_path, max(max(c_k(:)), 0));
        end
        fprintf('   max path ineq violation  : %.3e\n', max_path);
    else
        fprintf('   max path ineq violation  : (none)\n');
    end

    % 5. State bound violations
    x_lb = problem.bounds.x_lb(:);
    x_ub = problem.bounds.x_ub(:);
    lb_viol = max(max(bsxfun(@minus, x_lb, X), 0), [], 'all');
    ub_viol = max(max(bsxfun(@minus, X, x_ub), 0), [], 'all');
    fprintf('   max state LB violation   : %.3e\n', lb_viol);
    fprintf('   max state UB violation   : %.3e\n', ub_viol);

    % 6. Control bound violations
    u_lb = problem.bounds.u_lb(:);
    u_ub = problem.bounds.u_ub(:);
    u_lb_viol = max(max(bsxfun(@minus, u_lb, U), 0), [], 'all');
    u_ub_viol = max(max(bsxfun(@minus, U, u_ub), 0), [], 'all');
    fprintf('   max control LB violation : %.3e\n', u_lb_viol);
    fprintf('   max control UB violation : %.3e\n', u_ub_viol);

    fprintf('   -------------------------------------------\n\n');
end