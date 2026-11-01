function w = paper_cost_weights(m, nu)
%PAPER_COST_WEIGHTS  CDC 2025 Section V.A cost weights (canonical defaults).
%
%   Paper: Kidane, Motayed, Wang — IEEE CDC 2025,
%   doi:10.1109/CDC57313.2025.11312610
%
%   Quoted numerical values (Section V.A / Hamiltonian discussion):
%     alpha = 5, beta = 400, gamma1 = 4000, gamma2 = 4000 (two-drug).
%
%   w = paper_cost_weights(m, nu)
%     m  - number of trait states (typically 21)
%     nu - number of controls (1 or 2)

    if nargin < 1 || isempty(m), m = 21; end
    if nargin < 2 || isempty(nu), nu = 1; end

    w = struct();
    w.alpha_vec     = 5   * ones(m, 1);
    w.beta_vec      = 400 * ones(m, 1);
    w.gamma1        = 4000;
    w.gamma2        = (nu > 1) * 4000;
    w.lambda_smooth = 1e2;
    w.profile       = 'paper';
    w.source        = 'CDC 2025 Section V.A';
end
