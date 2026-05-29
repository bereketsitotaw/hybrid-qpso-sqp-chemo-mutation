function [u_best, x_best, history, elapsed, T_best] = ocp_qpso(problem)
%OCP_QPSO  QPSO global search over control space.
%
%   [u_best, x_best, history, elapsed, T_best] = ocp_qpso(problem)
%
%   Searches for good controls using Quantum Particle Swarm Optimization.
%   Forward-simulates each particle with ODE45 to evaluate cost.
%   Returns best controls, corresponding state trajectory, and real
%   convergence history.
%
%   Outputs:
%     u_best   [nu, N+1]      best control trajectory
%     x_best   [nx, N+1]      forward-simulated states from u_best
%     history  [maxIter, 1]   gBestCost at each iteration (real, not fake)
%     elapsed  scalar         wall-clock seconds
%     T_best   scalar         best final time length

    tic;

    N    = problem.grid.N;
    nx   = problem.nx;
    nu   = problem.nu;
    opts = problem.options.qpso;

    is_free_time = strcmp(problem.time.mode, 'free');

    % Setup Time, Variables, and Boundaries
    u_lb_qpso = problem.bounds.u_lb_qpso(:)';
    u_ub_qpso = problem.bounds.u_ub_qpso(:)';
    if is_free_time
        nvar = nu * (N + 1) + 1;
        T_fixed = [];
        t_vec_fixed = [];
        lb_flat = [repmat(u_lb_qpso, 1, N + 1), problem.time.T_bounds(1)];
        ub_flat = [repmat(u_ub_qpso, 1, N + 1), problem.time.T_bounds(2)];
    else
        T_fixed = problem.time.T_fixed;
        t_vec_fixed = linspace(0, T_fixed, N + 1);
        nvar = nu * (N + 1);
        lb_flat = repmat(u_lb_qpso, 1, N + 1);
        ub_flat = repmat(u_ub_qpso, 1, N + 1);
    end

    % Initialize swarm
    nPart     = opts.swarmSize;
    particles = lb_flat + (ub_flat - lb_flat) .* rand(nPart, nvar);
    
    % Fast-fail check: Evaluate particle 1 to surface structural errors
    fprintf('   -> Initializing QPSO & validating model dynamics...\n');
    [~, base_cost_est] = particle_cost(particles(1, :), T_fixed, t_vec_fixed, problem, is_free_time, true, 1e4);
    
    % Auto-scale penalties to avoid drowning out native cost
    penalty_scale = max(1e4, abs(base_cost_est) * 1000);

    pBest     = particles;
    pBestCost = inf(nPart, 1);
    gBest     = particles(1, :);
    gBestCost = inf;
    history   = zeros(opts.maxIter, 1);

    % Main QPSO loop
    for iter = 1:opts.maxIter
        parfor i = 1:nPart
            cost = particle_cost(particles(i, :), T_fixed, t_vec_fixed, problem, is_free_time, false, penalty_scale);
            if cost < pBestCost(i)
                pBestCost(i) = cost;
                pBest(i, :)  = particles(i, :);
            end
        end

        [minCost, idx] = min(pBestCost);
        if minCost < gBestCost
            gBestCost = minCost;
            gBest     = pBest(idx, :);
        end
        history(iter) = gBestCost;

        if mod(iter, 20) == 0
            fprintf('      QPSO iter %d/%d | best = %.4f\n', ...
                    iter, opts.maxIter, gBestCost);
        end

        % QPSO particle update
        mBest = mean(pBest, 1);
        for i = 1:nPart
            phi   = rand(1, nvar);
            p_i   = phi .* pBest(i, :) + (1 - phi) .* gBest;
            u_val = rand(1, nvar);
            L_val = 2 * abs(particles(i, :) - mBest);

            noise_amp   = 0.1;
            eta         = randn(1, nvar);
            unif        = rand(1, nvar) - 0.5;
            exploration = noise_amp * eta .* unif;

            particles(i, :) = p_i + (-1).^ceil(0.5 + rand(1, nvar)) .* L_val ...
                                .* log(1 ./ u_val) + exploration;
            
            % CLIP immediately to search bounds
            particles(i, :) = max(min(particles(i, :), ub_flat), lb_flat);
        end
    end

    % Recover best trajectory
    if is_free_time
        T_best = gBest(end);
        u_best = reshape(gBest(1:end-1), nu, N + 1);
        t_vec_best = linspace(0, T_best, N + 1);
    else
        T_best = T_fixed;
        u_best = reshape(gBest, nu, N + 1);
        t_vec_best = t_vec_fixed;
    end

    [x_best, ok] = ocp_forward_sim(u_best, T_best, t_vec_best, problem, false);
    if ~ok
        warning('ocp_qpso: best particle failed forward simulation.');
        x_best = repmat(problem.bounds.x0_lb, 1, N + 1);
    end

    elapsed = toc;
    fprintf('   -> QPSO complete in %.2f s.\n', elapsed);
end

%==========================================================================
%                         LOCAL FUNCTIONS
%==========================================================================

%--- Cost for one particle (with soft penalties) ---
function [cost, base_cost] = particle_cost(z_flat, fixed_T, fixed_t_vec, problem, is_free_time, fail_fast, penalty_scale)
    N  = problem.grid.N;
    nu = problem.nu;

    if is_free_time
        T = z_flat(end);
        t_vec = linspace(0, T, N + 1);
        U = reshape(z_flat(1:end-1), nu, N + 1);
    else
        T = fixed_T;
        t_vec = fixed_t_vec;
        U = reshape(z_flat, nu, N + 1);
    end

    if fail_fast
        % Evaluate without try-catch to surface syntax/dimension errors immediately
        [X, ok] = ocp_forward_sim(U, T, t_vec, problem, true);
        if ~ok
            cost = 1e12; base_cost = 1e12; return;
        end
        base_cost = ocp_cost(t_vec, X, U, T, problem);
        cost = base_cost;
        cost = add_penalties(cost, X, U, T, t_vec, problem, penalty_scale);
        return;
    end

    try
        % Forward simulate
        [X, ok] = ocp_forward_sim(U, T, t_vec, problem, false);
        if ~ok
            cost = 1e12; base_cost = 1e12;
            return;
        end

        % Base Bolza cost
        base_cost = ocp_cost(t_vec, X, U, T, problem);
        cost = base_cost;
        cost = add_penalties(cost, X, U, T, t_vec, problem, penalty_scale);

        if ~isfinite(cost)
            cost = 1e12;
        end
    catch
        cost = 1e12; base_cost = 1e12;
    end
end

%--- Utility to compute constraints penalties ---
function cost = add_penalties(cost, X, U, T, t_vec, problem, penalty_scale)
    N = problem.grid.N;
    
    % Soft penalty: boundary equality constraints
    if ~isempty(problem.constraints.boundary_eq)
        ceq = problem.constraints.boundary_eq(X(:,1), X(:,end), T, problem.param);
        cost = cost + penalty_scale * sum(ceq(:).^2);
    end

    % Soft penalty: boundary inequality constraints
    if ~isempty(problem.constraints.boundary_ineq)
        c = problem.constraints.boundary_ineq(X(:,1), X(:,end), T, problem.param);
        cost = cost + penalty_scale * sum(max(0, c(:)).^2);
    end

    % Soft penalty: path inequality constraints
    if ~isempty(problem.constraints.path_ineq)
        for k = 1:(N + 1)
            c_k = problem.constraints.path_ineq(t_vec(k), X(:,k), U(:,k), problem.param);
            cost = cost + (penalty_scale * 10) * sum(max(0, c_k(:)).^2);
        end
    end

    % Soft penalty: state bounds
    x_lb = problem.bounds.x_lb;
    x_ub = problem.bounds.x_ub;
    if any(isfinite(x_lb(:))) || any(isfinite(x_ub(:)))
        viol_lb = max(0, bsxfun(@minus, x_lb(:), X));
        viol_ub = max(0, bsxfun(@minus, X, x_ub(:)));
        cost = cost + penalty_scale * (sum(viol_lb(:).^2) + sum(viol_ub(:).^2));
    end
end
