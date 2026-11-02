function test_experimental_profile_warns()
%TEST_EXPERIMENTAL_PROFILE_WARNS  Non-paper profile loads and is marked.

    warn_state = warning('query', 'define_chemotherapy_problem:NonPaperWeights');
    cleaner = onCleanup(@() warning(warn_state.state, 'define_chemotherapy_problem:NonPaperWeights'));
    warning('off', 'define_chemotherapy_problem:NonPaperWeights');

    p = define_chemotherapy_problem('single', 'weight_profile', 'experimental');
    assert(strcmp(p.weights.profile, 'experimental_local_drift'));
    assert(abs(p.weights.alpha_vec(1) - 1) < 1e-12);
    assert(abs(p.weights.beta_vec(1) - 0.1) < 1e-12);
    assert(abs(p.weights.gamma1 - 20) < 1e-12);

    % Paper remains default
    p0 = define_chemotherapy_problem('single');
    assert(strcmp(p0.weights.profile, 'paper'));
end
