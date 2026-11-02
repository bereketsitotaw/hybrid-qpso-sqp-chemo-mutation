function test_dynamics_finite()
%TEST_DYNAMICS_FINITE  One-step dynamics and costs are finite.

    p = define_chemotherapy_problem('cosine_gaussian');
    x0 = p.bounds.x0;
    u0 = [1; 1];
    f = p.model.f(0, x0, u0, p.param);
    assert(all(isfinite(f)), 'dynamics produced non-finite values');
    L = p.cost.L(0, x0, u0, p.param, p.param, p.weights);
    Phi = p.cost.Phi([], x0, 10, p.param, p.weights);
    assert(isfinite(L) && isfinite(Phi));
    assert(L > 0 && Phi > 0);
end
