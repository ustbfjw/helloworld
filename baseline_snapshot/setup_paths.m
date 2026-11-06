function setup_paths(root)
%SETUP_PATHS 显式加入模块，避免递归加入同名备份和结果目录。
if nargin < 1, root = fileparts(mfilename('fullpath')); end
addpath(root);
folders = {'00_config','01_trajectory','02_iscat_imaging','03_inta_analysis','04_results','utils','tests'};
for k=1:numel(folders)
    folder = fullfile(root,folders{k});
    if exist(folder,'dir'), addpath(folder); end
end
end
