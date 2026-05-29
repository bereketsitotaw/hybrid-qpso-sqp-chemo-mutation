function reproduce_paper_figures(which)
%REPRODUCE_PAPER_FIGURES Regenerate comparison figures from Section V–VI.
%
%   reproduce_paper_figures()      % Figs. 4–6 (uses cached solutions when available)
%   reproduce_paper_figures(4)     % Fig. 4 only
%   reproduce_paper_figures('all') % Figs. 1–6 (1–3 via main, 4–6 comparison)

    if nargin < 1 || isempty(which)
        which = [4, 5, 6];
    end

    clc;
    close all;
    init_project();

    repo_root = fileparts(mfilename('fullpath'));
    figures_dir = fullfile(repo_root, 'results', 'figures');
    if ~exist(figures_dir, 'dir')
        mkdir(figures_dir);
    end

    if ischar(which) || isstring(which)
        key = lower(char(which));
        if strcmp(key, 'all')
            fprintf('Running individual scenarios (Figs. 1–3)...\n');
            run_scenarios('all');
            which = [4, 5, 6];
        else
            which = str2double(key);
        end
    end

    fig_list = which(:)';
    for f = fig_list
        switch f
            case 1
                run_scenarios('single');
            case 2
                run_scenarios('cosine_gaussian');
            case 3
                run_scenarios('cosine_sine');
            case 4
                reproduce_figure_4(figures_dir);
            case 5
                reproduce_figure_5(figures_dir);
            case 6
                reproduce_figure_6(figures_dir);
            otherwise
                warning('Unknown figure index %d — skipped.', f);
        end
    end

    fprintf('\nFigure reproduction complete. See: %s\n', figures_dir);
end
