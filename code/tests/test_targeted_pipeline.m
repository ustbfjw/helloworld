function metrics=test_targeted_pipeline(outdir)
%TEST_TARGETED_PIPELINE 定向runner的外置归一化与生产内部归一化数值等价。
if nargin<1, outdir=''; end
root=fileparts(fileparts(mfilename('fullpath'))); setup_paths(root);
p=apply_scenario(sim_params('quick'),'mono100'); p.sim.n_frames=500; p.sim.n_videos=1;
p.sim.L_xy=.5e-6; p.sim.z_init=[.9,1.1]*1e-6; p.ana.bootstrap_n=0;
p.cam.flicker_amplitude=.08; p.cam.flicker_frequency=423.7;
traj=simulate_trajectories(p,987654); [video,~]=render_iscat_video(traj,p,ipsf_lookup_table(p));
[gain,~]=targeted_method_gain(video,p,'full');
[internal,exi]=inta_track_video(video,p);
pm=p; pm.ana.background_mode='median';
corrected=single(video)./reshape(single(gain),1,1,[]);
[external,exe]=inta_track_video(corrected,pm);
assert(~isempty(internal),'测试必须实际检测到轨迹。');
assert(isequaln(internal,external),'外置归一化改变轨迹。');
assert(isequaln(exi.C,exe.C) && isequaln(exi.R,exe.R),'外置归一化改变分析图像。');
pn=targeted_method_gain(video,p,'pn');
expected=double(squeeze(sum(sum(single(video),1),2))); expected=expected/median(expected);
assert(isequal(pn,expected),'PN定义发生变化。');
assert(isequal(targeted_method_gain(video,p,'median'),ones(p.sim.n_frames,1)));
metrics=struct('passed',true,'n_tracks',numel(internal),'frames',p.sim.n_frames, ...
    'external_normalization_equal',true,'pn_definition_equal',true, ...
    'scope','fixture with centered initial particles, used only for numerical integration test, not formal sizing evidence');
if ~isempty(outdir)
    if ~exist(outdir,'dir'), mkdir(outdir); end
    fid=fopen(fullfile(outdir,'targeted_pipeline_tests.json'),'w','n','UTF-8');
    cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%s',jsonencode(metrics,PrettyPrint=true));
end
disp(metrics);
end
