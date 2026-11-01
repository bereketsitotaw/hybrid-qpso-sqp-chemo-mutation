function main(arg1, arg2, varargin)
%MAIN Hybrid QPSO-SQP chemotherapy (trait model, single/double drug).
%
%   main                  % all scenarios
%   main('single')
%   main('double')        % three combination therapies
%   main('figures')       % Figs. 4–6 (uses cached results when present)
%   main('figures', 5)    % one figure only
%   main('smoke')         % fast non-paper solve (diagnostics only)
%   main('test')          % run MATLAB smoke tests in tests/

    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root, 'scripts'));
    addpath(fullfile(root, 'configs'));
    clc; close all;

    if nargin >= 1 && (ischar(arg1) || isstring(arg1))
        key = lower(char(arg1));
        if any(strcmp(key, {'figures', 'fig'}))
            if nargin < 2, which = [4, 5, 6]; else, which = arg2; end
            chemo_pipeline(root, 'figures', which);
            return;
        end
        if strcmp(key, 'test')
            run(fullfile(root, 'tests', 'run_smoke_tests.m'));
            return;
        end
        if strcmp(key, 'smoke')
            chemo_pipeline(root, 'smoke', 'single');
            return;
        end
    end

    if nargin < 1, arg1 = 'all'; end
    chemo_pipeline(root, 'run', arg1, varargin{:});
end
