function test_reference_metrics_pin()
%TEST_REFERENCE_METRICS_PIN  Pinned metrics match reference .mat when present.

    root = fileparts(fileparts(mfilename('fullpath')));
    metrics = reference_metrics();
    mat_path = fullfile(root, metrics.mat_file);
    assert(exist(mat_path, 'file') == 2, 'missing reference mat: %s', mat_path);

    d = load(mat_path);
    assert(isfield(d, 'solution'));
    s = d.solution;
    assert(abs(s.T_final - metrics.T_final) < 1e-9);
    assert(abs(s.cost - metrics.cost) < 1e-2);
    assert(s.exitflag == metrics.exitflag);
    final_avg = mean(s.states(:, end));
    assert(abs(final_avg - metrics.final_avg_tumor) < 1e-5);
    dose1 = trapz(s.time_vec, s.controls(1, :));
    assert(abs(dose1 - metrics.dose1) < 1e-4);
end
