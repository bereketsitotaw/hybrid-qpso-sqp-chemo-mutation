function cfg = scenario_single()
%SCENARIO_SINGLE Fig. 1 — single medication, no mutation (PMP benchmark).

    cfg = struct();
    cfg.id = 'single';
    cfg.label = 'Single-drug PMP benchmark (Fig. 1)';
    cfg.scenario = 'single';
    cfg.medication = 'single';
    cfg.mutation = false;
end
