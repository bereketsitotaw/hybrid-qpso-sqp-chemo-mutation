function test_scenario_aliases()
%TEST_SCENARIO_ALIASES  Legacy IDs and string aliases resolve.

    a = define_chemotherapy_problem(1);
    b = define_chemotherapy_problem('fig1');
    c = define_chemotherapy_problem('single');
    assert(strcmp(a.meta.scenario, 'single'));
    assert(strcmp(b.meta.scenario, 'single'));
    assert(strcmp(c.meta.scenario, 'single'));
    assert(a.nu == 1 && ~a.meta.mutation);

    d = define_chemotherapy_problem(2);
    assert(strcmp(d.meta.scenario, 'cosine_sine'));
    assert(d.nu == 2 && d.meta.mutation);

    e = define_chemotherapy_problem('cg');
    assert(strcmp(e.meta.scenario, 'cosine_gaussian'));
end
