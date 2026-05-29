function [X, ok] = ocp_forward_sim(U, T, t_vec, problem, throw_errors)
%OCP_FORWARD_SIM  Forward-simulate dynamics using ODE45.
%
%   [X, ok] = ocp_forward_sim(U, T, t_vec, problem, throw_errors)
%
%   Uses griddedInterpolant internally to achieve high performance compared
%   to pointwise interp1 calls inside the ODE step.
%
%   Inputs:
%     U            [nu, N+1] Matrix of control trajectories.
%     T            scalar    Final time.
%     t_vec        [1, N+1]  Time grid array.
%     problem      struct    Full problem configuration struct.
%     throw_errors bool      If true, bypasses try-catch so model errors throw immediately.
%
%   Outputs:
%     X            [nx, N+1] State trajectory outputs.
%     ok           bool      True if integration completed without throwing/failing.

    if nargin < 5
        throw_errors = false;
    end

    nx = problem.nx;
    nu = problem.nu;
    x0 = problem.bounds.x0_lb(:);  % fixed initial state

    % Pre-build interpolants for all control channels to eliminate overhead
    U_interp = cell(1, nu);
    for i = 1:nu
        U_interp{i} = griddedInterpolant(t_vec, U(i, :), 'linear', 'nearest');
    end

    ode_opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-4);

    if throw_errors
        [~, x_traj] = ode45(@(t, x) ode_rhs_fast(t, x, U_interp, problem), ...
                             t_vec, x0, ode_opts);
        X  = x_traj';   % [nx, N+1]
        ok = true;
    else
        try
            [~, x_traj] = ode45(@(t, x) ode_rhs_fast(t, x, U_interp, problem), ...
                                 t_vec, x0, ode_opts);
            X  = x_traj';   % [nx, N+1]
            ok = true;
        catch
            X  = zeros(nx, length(t_vec));
            ok = false;
        end
    end
end

function dxdt = ode_rhs_fast(t, x, U_interp, problem)
    nu = problem.nu;
    u  = zeros(nu, 1);
    for i = 1:nu
        u(i) = U_interp{i}(t);
    end

    % Clip to bounds (matching legacy rk4 behavior for stability)
    u = max(bsxfun(@max, u, problem.bounds.u_lb), problem.bounds.u_lb);
    u = min(bsxfun(@min, u, problem.bounds.u_ub), problem.bounds.u_ub);

    dxdt = problem.model.f(t, x, u, problem.param);
end
