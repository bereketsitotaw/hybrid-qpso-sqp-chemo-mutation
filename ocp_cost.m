function J = ocp_cost(t_vec, X, U, T, problem)
%OCP_COST  Evaluate the Bolza cost functional.
%
%   J = ocp_cost(t_vec, X, U, T, problem)
%
%   Computes J = integral_0^T L(t,x,u,udot,...) dt + Phi(x0,xF,T,...)
%   Uses problem.cost.L_vectorized if available for speed.

    N  = problem.grid.N;
    h  = T / N;
    nu = problem.nu;

    % Compute control rates (numerical gradient)
    Udot = zeros(size(U));
    for i = 1:nu
        Udot(i, :) = gradient(U(i, :), h);
    end

    % Running cost
    if isfield(problem.cost, 'L_vectorized') && ~isempty(problem.cost.L_vectorized)
        L_vec = problem.cost.L_vectorized(t_vec, X, U, Udot, problem.param, problem.weights);
    else
        Npts  = length(t_vec);
        L_vec = zeros(1, Npts);
        for k = 1:Npts
            L_vec(k) = problem.cost.L(t_vec(k), X(:,k), U(:,k), Udot(:,k), ...
                                       problem.param, problem.weights);
        end
    end
    J = trapz(t_vec, L_vec);

    % Smoothness penalty: integral of udot^2
    if isfield(problem.weights, 'lambda_smooth') && problem.weights.lambda_smooth > 0
        J_smooth = 0;
        for i = 1:nu
            J_smooth = J_smooth + trapz(t_vec, Udot(i, :).^2);
        end
        J = J + problem.weights.lambda_smooth * J_smooth;
    end

    % Terminal cost
    if ~isempty(problem.cost.Phi)
        J = J + problem.cost.Phi(X(:,1), X(:,end), T, problem.param, problem.weights);
    end
end
