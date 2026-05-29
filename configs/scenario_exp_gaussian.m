function cfg = scenario_exp_gaussian()
%SCENARIO_EXP_GAUSSIAN Sec. V.B.3 — exponential + Gaussian drugs.

    cfg = struct();
    cfg.id = 'exp_gaussian';
    cfg.label = 'Exponential-Gaussian combination (Sec. V.B.3)';
    cfg.scenario = 'exp_gaussian';
    cfg.medication = 'double';
    cfg.mutation = true;
end
