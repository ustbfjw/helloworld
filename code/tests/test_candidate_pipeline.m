function metrics=test_candidate_pipeline(outdir)
%TEST_CANDIDATE_PIPELINE 实际追踪接口的内外校正等价；不作为物理精度证据。
root=fileparts(fileparts(mfilename('fullpath'))); setup_paths(root);
p=apply_scenario(sim_params('quick'),'mono100'); p.sim.n_frames=1000;
p.sim.n_videos=1; p.sim.L_xy=.5e-6; p.sim.z_init=[.9,1.1]*1e-6;
p.ana.bootstrap_n=0; p.ana.background_mode='spatial_validated';
traj=simulate_trajectories(p,981234); [raw,~]=render_iscat_video(traj,p,ipsf_lookup_table(p));
[X,Y]=meshgrid(linspace(-1,1,p.cam.Npx));
a=.045+.025*X+.015*Y.^2+.012*X.*Y;
t=(0:p.sim.n_frames-1)'/p.sim.fps; shape=sin(2*pi*423.7*t+.37);
% 此处ADU之后乘增益仅用于软件等价fixture；正式生成器在曝光内加照明。
video=single(raw).*single(1+a.*reshape(shape,1,1,[]));
[model,info]=estimate_validated_spatial_illumination(video,p);
assert(model.is_spatial && info.is_spatial,'等价测试必须覆盖实际空间分支。');
[inside,exi]=inta_track_video(video,p);
pm=p; pm.ana.background_mode='median';
[outside,exo]=inta_track_video(apply_candidate_illumination(video,model),pm);
assert(~isempty(inside),'必须检测到实际轨迹以检查接口。');
assert(isequaln(inside,outside) && isequaln(exi.C,exo.C) && isequaln(exi.R,exo.R), ...
    '内外空间校正未保持逐位一致。');
[fallback,~]=estimate_validated_spatial_illumination(raw,p);
assert(~fallback.is_spatial,'静态fixture应回退。');
[ti,ei]=inta_track_video(raw,p);
pm.ana.background_mode='spectral_common'; [to,eo]=inta_track_video(raw,pm);
assert(isequaln(ti,to) && isequaln(ei.C,eo.C) && isequaln(ei.R,eo.R));
metrics=struct('passed',true,'frames',p.sim.n_frames,'n_tracks',numel(inside), ...
    'spatial_internal_external_equal',true,'fallback_production_equal',true, ...
    'scope','software equivalence fixture only; post-ADC gain is not formal physical illumination');
if ~exist(outdir,'dir'), mkdir(outdir); end
fid=fopen(fullfile(outdir,'candidate_pipeline_tests.json'),'w','n','UTF-8');
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(metrics,PrettyPrint=true)); disp(metrics);
end
