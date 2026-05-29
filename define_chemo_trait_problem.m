function problem = define_chemo_trait_problem(varargin)
%DEFINE_CHEMO_TRAIT_PROBLEM  Legacy wrapper — use define_chemotherapy_problem.
%
%   problem = define_chemo_trait_problem()     % Scenario 1 (single drug)
%   problem = define_chemo_trait_problem(2)    % Scenario 2 (cosine-sine)
%   problem = define_chemo_trait_problem(3)    % Scenario 3 (exp-gaussian)
%
%   Prefer: define_chemotherapy_problem('cosine_gaussian'), etc.

    if isempty(varargin)
        legacy_id = 1;
    else
        legacy_id = varargin{1};
    end

    map = {'single', 'cosine_sine', 'exp_gaussian'};
    if legacy_id < 1 || legacy_id > numel(map)
        error('Invalid scenario_id. Choose 1, 2, or 3.');
    end

    problem = define_chemotherapy_problem(map{legacy_id});
end
