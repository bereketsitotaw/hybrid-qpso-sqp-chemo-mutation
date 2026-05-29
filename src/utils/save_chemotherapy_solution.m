function save_chemotherapy_solution(solution, summary, problem, filepath)
%SAVE_CHEMOTHERAPY_SOLUTION Persist solution struct to disk.

    payload = struct();
    payload.solution = solution;
    payload.summary = summary;
    payload.problem_meta = problem.meta;
    payload.saved_at = datetime('now');
    save(filepath, '-struct', 'payload');
end
