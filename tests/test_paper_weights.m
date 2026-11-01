function test_paper_weights()
%TEST_PAPER_WEIGHTS  Default profile must match CDC 2025 Section V.A.

    p1 = define_chemotherapy_problem('single');
    assert(strcmp(p1.weights.profile, 'paper'));
    assert(abs(p1.weights.alpha_vec(1) - 5) < 1e-12);
    assert(abs(p1.weights.beta_vec(1) - 400) < 1e-12);
    assert(abs(p1.weights.gamma1 - 4000) < 1e-12);
    assert(abs(p1.weights.gamma2 - 0) < 1e-12);

    p2 = define_chemotherapy_problem('cosine_gaussian');
    assert(abs(p2.weights.gamma1 - 4000) < 1e-12);
    assert(abs(p2.weights.gamma2 - 4000) < 1e-12);

    wcfg = paper_cost_weights(21, 1);
    assert(abs(wcfg.alpha_vec(1) - 5) < 1e-12);
end
