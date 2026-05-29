function cfg = scenario_cosine_gaussian()
%SCENARIO_COSINE_GAUSSIAN Sec. V.B.1 — cosine-squared + Gaussian drugs.

    cfg = struct();
    cfg.id = 'cosine_gaussian';
    cfg.label = 'Cosine-Gaussian combination (Sec. V.B.1)';
    cfg.scenario = 'cosine_gaussian';
    cfg.medication = 'double';
    cfg.mutation = true;
end
