function pmp = ocp_pmp(solution, summary, problem)
%OCP_PMP  Optional PMP optimality diagnostics.
%
%   pmp = ocp_pmp(solution, summary, problem)
%
%   Extracts discrete costates from fmincon multipliers (eqnonlin),
%   computes the Hamiltonian, and calls problem.pmp_derivatives for
%   switching-function analysis.
%
%   This is post-processing only — never called during optimization.

    pmp = struct();

    if ~isfield(problem, 'pmp_derivatives') || isempty(problem.pmp_derivatives)
        fprintf('      PMP: no pmp_derivatives callback. Providing numerical placeholders.\n');
        % Still continue to provide zeroed placeholders so plots aren't empty
    end

    nx = problem.nx;
    nu = problem.nu;
    N  = problem.grid.N;
    T  = solution.T_final;
    h  = T / N;

    % --- Extract defect multipliers from fmincon ---
    lambda_all = summary.lambda.eqnonlin(:);
    n_defects  = nx * N;

    if length(lambda_all) < n_defects
        warning('PMP: not enough eqnonlin multipliers (%d < %d).', ...
                length(lambda_all), n_defects);
        return;
    end

    MU = reshape(lambda_all(1:n_defects), nx, N);

    % --- Undo defect scaling to recover physical multipliers ---
    % The NLP scales defects by problem.scaling.defect_scale (per-state),
    % so the KKT multipliers are on the scaled constraints.
    % Physical multiplier = scaled_multiplier * defect_scale_factor.
    ds = ones(nx, 1);
    if isfield(problem, 'scaling') && isfield(problem.scaling, 'defect_scale')
        ds = problem.scaling.defect_scale(:);
    end
    MU_physical = bsxfun(@times, MU, ds);  % [nx, N]

    % For trapezoidal collocation, the continuous-time costate is:
    %   lambda(t_k) approx -mu_k / h
    % We extend to N+1 nodes by repeating the last value.
    Lambda_raw = -(1/h) * [MU_physical, MU_physical(:, end)];   % [nx, N+1]

    if any(isnan(Lambda_raw(:))) || any(isinf(Lambda_raw(:)))
        warning('PMP: NaN/Inf in raw costates.');
    end

    % --- Hamiltonian: H = L + lambda' * f ---
    X     = solution.states;
    U     = solution.controls;
    t_vec = solution.time_vec;

    F = ocp_dynamics(t_vec, X, U, problem);

    Udot = zeros(size(U));
    for i = 1:nu
        Udot(i, :) = gradient(U(i, :), h);
    end

    % Running cost at each node
    if isfield(problem.cost, 'L_vectorized') && ~isempty(problem.cost.L_vectorized)
        L_vec = problem.cost.L_vectorized(t_vec, X, U, Udot, ...
                                           problem.param, problem.weights);
    else
        L_vec = zeros(1, N + 1);
        for k = 1:(N + 1)
            L_vec(k) = problem.cost.L(t_vec(k), X(:,k), U(:,k), Udot(:,k), ...
                                       problem.param, problem.weights);
        end
    end

    H = L_vec + sum(Lambda_raw .* F, 1);

    % Add smoothness penalty contribution to Hamiltonian
    % (the NLP cost includes lambda_smooth * integral(udot^2),
    %  so the pointwise Hamiltonian should include this term)
    if isfield(problem.weights, 'lambda_smooth') && problem.weights.lambda_smooth > 0
        for i = 1:nu
            H = H + problem.weights.lambda_smooth * Udot(i, :).^2;
        end
    end

    % Smooth costates for visualization
    Lambda_smooth = zeros(size(Lambda_raw));
    for i = 1:nx
        Lambda_smooth(i, :) = smoothdata(Lambda_raw(i, :), 'gaussian', ...
                                          max(3, round(N/20)));
    end

    % --- Problem-specific PMP derivatives ---
    % Initialize with placeholders of correct size [nu x N+1]
    dH_du    = zeros(nu, N + 1);
    d2H_du2  = zeros(nu, N + 1);
    dH_du_dt = zeros(nu, N + 1);

    if isfield(problem, 'pmp_derivatives') && ~isempty(problem.pmp_derivatives)
        try
            [dH_du_call, d2H_du2_call, dH_du_dt_call] = problem.pmp_derivatives( ...
                X, U, Lambda_raw, problem.param, problem.weights, t_vec);
            
            % Assign what we got, keeping placeholders if outputs are empty
            if ~isempty(dH_du_call),    dH_du    = dH_du_call;    end
            if ~isempty(d2H_du2_call),  d2H_du2  = d2H_du2_call;  end
            if ~isempty(dH_du_dt_call), dH_du_dt = dH_du_dt_call; end
        catch ME
            warning('PMP: pmp_derivatives failed: %s. Trying smoothed.', ME.message);
            try
                [dH_du_call, d2H_du2_call, dH_du_dt_call] = problem.pmp_derivatives( ...
                    X, U, Lambda_smooth, problem.param, problem.weights, t_vec);
                if ~isempty(dH_du_call),    dH_du    = dH_du_call;    end
                if ~isempty(d2H_du2_call),  d2H_du2  = d2H_du2_call;  end
                if ~isempty(dH_du_dt_call), dH_du_dt = dH_du_dt_call; end
            catch ME2
                warning('PMP: also failed with smoothed costates: %s', ME2.message);
            end
        end
    end

    % --- Numerical Fallbacks (Robustness for Plotting) ---
    % If dH_du_dt was not provided but dH_du exists, compute it numerically
    is_zero_dt = all(dH_du_dt(:) == 0);
    is_any_u   = any(dH_du(:) ~= 0);
    if is_zero_dt && is_any_u
        for i = 1:nu
            dH_du_dt(i, :) = gradient(dH_du(i, :), h);
        end
    end

    % If d2H_du2 was not provided (all zeros), we keep it as zeros 
    % but ensured it is the correct size [nu, N+1].

    % --- Package ---
    pmp.costates     = Lambda_smooth;
    pmp.costates_raw = Lambda_raw;
    pmp.hamiltonian  = H;
    pmp.dH_du        = dH_du;
    pmp.d2H_du2      = d2H_du2;
    pmp.dH_du_dt     = dH_du_dt;

    if ~isempty(dH_du)
        fprintf('      PMP: max|dH/du| = %.4e\n', max(abs(dH_du(:))));
    end
    fprintf('      Hamiltonian: mean=%.4f, std=%.4f\n', mean(H), std(H));
end
