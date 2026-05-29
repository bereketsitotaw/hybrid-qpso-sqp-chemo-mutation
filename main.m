function main(arg1, arg2, varargin)
%MAIN Hybrid QPSO-SQP chemotherapy (trait model, single/double drug).
%
%   main                  % all scenarios
%   main('single')
%   main('double')        % three combination therapies
%   main('figures')       % Figs. 4–6 (uses cached results when present)
%   main('figures', 5)    % one figure only

    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root, 'scripts'));
    clc; close all;

    if nargin >= 1 && (ischar(arg1) || isstring(arg1)) && any(strcmpi(char(arg1), {'figures', 'fig'}))
        if nargin < 2, which = [4, 5, 6]; else, which = arg2; end
        chemo_pipeline(root, 'figures', which);
    else
        if nargin < 1, arg1 = 'all'; end
        chemo_pipeline(root, 'run', arg1, varargin{:});
    end
end
