function results = main_run_all(preset, output_dir, light_condition)
%MAIN_RUN_ALL 一键运行三组单分散样品：轨迹→相机→iNTA→中文结果。
% main_run_all             正式运行：30nm组20段，100/200nm各60段；每段5000帧。
% main_run_all('quick')    小样本演示，不用于精度验收。
% main_run_all('full',output_dir) 指定独立输出目录。
if nargin < 1, preset = 'full'; end
if nargin < 3, light_condition='static'; end
root = fileparts(mfilename('fullpath')); setup_paths(root);
p0 = sim_params(preset);
if strcmpi(light_condition,'periodic')
    p0.cam.flicker_amplitude=0.08;
    p0.out.dir=fullfile(root,'04_results','output_periodic');
elseif ~strcmpi(light_condition,'static')
    error('light_condition只支持static或periodic。');
end
if nargin >= 2 && ~isempty(output_dir), p0.out.dir = output_dir; end
if ~exist(p0.out.dir,'dir'), mkdir(p0.out.dir); end
print_setup(p0);
tab = ipsf_lookup_table(p0);
if p0.out.make_plots, plot_ipsf_examples(tab,p0); end
results = struct();
names = {'mono30','mono100','mono200'};
for k = 1:numel(names)
    p = apply_scenario(p0,names{k});
    fprintf('\n===== %s：%d段×%d帧 =====\n',names{k},p.sim.n_videos,p.sim.n_frames);
    [tracks,res,extra] = run_inta_pipeline(p,tab,true);
    print_summary(tracks,res,p,names{k});
    if p.out.make_plots
        ie = find(cellfun(@(e) ~isempty(e.frame),extra.ex),1);
        if ~isempty(ie)
            % 三维轨迹与相机示例来自同一视频，便于逐帧核对。
            plot_trajectories_3d(extra.traj{ie},p,names{k});
            plot_example_frames(extra.ex{ie},p,names{k});
            sample=struct('traj',extra.traj{ie},'ex',extra.ex{ie},'meta',extra.meta{ie},'video',ie);
            plot_validation_diagnostics(sample,tracks,p,names{k},res.sel);
            plot_illumination_diagnostics(sample,p,names{k});
        else
            plot_trajectories_3d(extra.traj{1},p,names{k});
        end
        plot_results(tracks,res,p,names{k});
    end
    results.(names{k}) = struct('tracks',tracks,'res',res,'p',p);
end
plot_group_comparison(results,p0);
save(fullfile(p0.out.dir,'inta_results.mat'),'results','p0','-v7.3');
fprintf('\n三组计算完成，结果目录：%s\n',p0.out.dir);
end
