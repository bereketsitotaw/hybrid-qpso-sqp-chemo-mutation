function cfg = scenario_cosine_sine()
%SCENARIO_COSINE_SINE Sec. V.B.2 — cosine-squared + sine-squared drugs.

    cfg = struct();
    cfg.id = 'cosine_sine';
    cfg.label = 'Cosine-Sine combination (Sec. V.B.2)';
    cfg.scenario = 'cosine_sine';
    cfg.medication = 'double';
    cfg.mutation = true;
end
