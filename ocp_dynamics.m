function F = ocp_dynamics(t_vec, X, U, problem)
%OCP_DYNAMICS  Evaluate dynamics at all grid nodes.
%
%   F = ocp_dynamics(t_vec, X, U, problem)
%
%   Uses problem.model.f_vectorized if available for speed,
%   otherwise loops calling problem.model.f pointwise.
%
%   Inputs:
%     t_vec   [1, N+1]   time grid
%     X       [nx, N+1]  state matrix
%     U       [nu, N+1]  control matrix
%     problem struct     problem definition
%
%   Output:
%     F       [nx, N+1]  dynamics dxdt at each node

    if isfield(problem.model, 'f_vectorized') && ~isempty(problem.model.f_vectorized)
        F = problem.model.f_vectorized(t_vec, X, U, problem.param);
    else
        Npts = size(X, 2);
        nx   = problem.nx;
        F    = zeros(nx, Npts);
        for k = 1:Npts
            F(:, k) = problem.model.f(t_vec(k), X(:,k), U(:,k), problem.param);
        end
    end
end
