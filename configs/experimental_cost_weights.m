function w = experimental_cost_weights(m, nu)
%EXPERIMENTAL_COST_WEIGHTS  Non-paper exploratory trial weights (NOT validated).
%
%   These values were found as an uncommitted local edit on one checkout
%   (alpha 5→1, beta 400→0.1, gamma1 4000→20). They are objective-changing
%   and do NOT match CDC 2025 Section V.A. Kept only as a non-paper
%   exploratory / forensic archive — not validated for publication.
%   Default pipeline must use paper_cost_weights.
%
%   w = experimental_cost_weights(m, nu)

    if nargin < 1 || isempty(m), m = 21; end
    if nargin < 2 || isempty(nu), nu = 1; end

    w = struct();
    w.alpha_vec     = 1   * ones(m, 1);
    w.beta_vec      = 0.1 * ones(m, 1);
    w.gamma1        = 20;
    w.gamma2        = (nu > 1) * 4000;  % note: local drift left gamma2 unchanged
    w.lambda_smooth = 1e2;
    w.profile       = 'experimental_local_drift';
    w.source        = 'uncommitted local edit (NOT paper)';
end
