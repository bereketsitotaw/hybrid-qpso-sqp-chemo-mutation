function test_problem_struct()
%TEST_PROBLEM_STRUCT  Required OCP fields exist for solveOCP.

    p = define_chemotherapy_problem('single');
    req_top = {'name','nx','nu','time','grid','param','weights','model','cost','bounds','options'};
    for i = 1:numel(req_top)
        assert(isfield(p, req_top{i}), 'missing %s', req_top{i});
    end
    assert(isa(p.model.f, 'function_handle'));
    assert(isa(p.cost.L, 'function_handle'));
    assert(isa(p.cost.Phi, 'function_handle'));
    assert(p.nx == 21);
    assert(p.grid.N == 100);
    assert(abs(p.time.T_fixed - 10) < 1e-12);
    assert(all(p.bounds.u_ub == 3));
end
