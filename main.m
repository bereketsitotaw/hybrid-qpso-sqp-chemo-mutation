function main(scenario_id, varargin)
%MAIN Hybrid QPSO-SQP chemotherapy optimal control (trait model).
%
% Usage (cd to repo root in MATLAB):
%   main                  % all paper scenarios
%   main('single')        % Fig. 1 — single drug, PMP benchmark
%   main('cosine_gaussian')
%   main('cosine_sine')
%   main('exp_gaussian')
%   main('double')        % all combination-therapy scenarios
%
%   reproduce_paper_figures('all')   % Figs. 1–6
%   reproduce_paper_figures([5 6])  % comparison figures
%
%   main('cosine_gaussian', 'medication', 'single')  % force 1 control
%
% Paper: Optimal Control for Cancer Chemotherapy Using Hybrid QPSO

    if nargin < 1 || isempty(scenario_id)
        scenario_id = 'all';
    end

    clc;
    close all;
    init_project();
    run_scenarios(scenario_id, varargin{:});
end
