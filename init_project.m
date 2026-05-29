function repo_root = init_project()
%INIT_PROJECT Add chemotherapy OCP repo folders to the MATLAB path.

    repo_root = fileparts(mfilename('fullpath'));
    addpath(repo_root);
    addpath(fullfile(repo_root, 'configs'));
    addpath(fullfile(repo_root, 'src', 'utils'));
    addpath(fullfile(repo_root, 'tests'));
end
