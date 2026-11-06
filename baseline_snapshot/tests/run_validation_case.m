function record = run_validation_case(name,outdir,seed,nvideos,nframes,G0,exposure,nsub)
%RUN_VALIDATION_CASE 固定方案的单组端到端验证，可由独立MATLAB进程并行运行。
% 先估计再比较名义直径，不使用真值筛除轨迹或调整估计值。
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root); setup_paths(root);
p=apply_scenario(sim_params('full'),name);
if nargin<3, seed=10001; end
if nargin<4 || isempty(nvideos), nvideos=p.sim.n_videos; end
if nargin<5 || isempty(nframes), nframes=p.sim.n_frames; end
% 省略参数时保留该样品的采集设置，避免覆盖三组各自的固定方案。
if nargin<6 || isempty(G0), G0=p.mag.G0; end
if nargin<7 || isempty(exposure), exposure=p.cam.exposure; end
if nargin<8 || isempty(nsub), nsub=p.sim.n_sub; end
p.sim.seed=seed; p.sim.n_videos=nvideos; p.sim.n_frames=nframes;
p.mag.G0=G0; p.cam.exposure=exposure; p.sim.n_sub=nsub;
p.out.dir=outdir; p.out.make_plots=false; p.out.save_video=false;
p.out.save_diagnostics=false;
if ~exist(outdir,'dir'), mkdir(outdir); end
diary(fullfile(outdir,[name '_run.log'])); cleanup=onCleanup(@() diary('off')); %#ok<NASGU>
tic;
[tracks,res,extra]=run_inta_pipeline(p,[],true);
elapsed=toc;
g=res.group;
sat=sum(cellfun(@(m) m.n_sat,extra.meta));
valid=find(res.sel & isfinite([tracks.loc_rms_nm]));
loc_rms=nan;
if ~isempty(valid)
    w=[tracks(valid).n_loc];
    loc_rms=sqrt(sum(w.*[tracks(valid).loc_rms_nm].^2)/sum(w));
end
record=struct('scenario',name,'seed',seed,'n_videos',nvideos, ...
    'frames',nframes,'G0',G0,'exposure',exposure,'nsub',nsub, ...
    'nominal_nm',p.sample.nominal_nm,'diameter_nm',g.d_avgmsd, ...
    'bias_percent',100*(g.d_avgmsd/p.sample.nominal_nm-1), ...
    'ci_nm',g.d_ci95_interval,'n_tracks',g.n_tracks,'n_loc',g.n_loc, ...
    'videos_with_tracks',g.n_videos_with_tracks,'saturated_pixels',sat, ...
    'localization_rms_nm',loc_rms,'unmatched_selected',sum(res.sel & [tracks.truth_id]==0), ...
    'illumination_corrected_videos',sum(cellfun(@(e) e.illumination.applied,extra.ex)), ...
    'seconds',elapsed);
fprintf('\nVALIDATION %s d=%.6f bias=%+.4f%% CI=[%.3f %.3f] tracks=%d loc=%d videos=%d saturation=%d time=%.1fs\n', ...
    name,record.diameter_nm,record.bias_percent,record.ci_nm,g.n_tracks,g.n_loc,g.n_videos_with_tracks,sat,elapsed);
fprintf('MODELS linear=%.4f perp=%.4f waveform=%.4f\n', ...
    g.model_estimates.linear.d,g.model_estimates.perp.d,g.model_estimates.waveform.d);
% 保留首个有检测视频的真轨迹/示例帧，便于把3D、相机和检测位置对照。
ie=find(cellfun(@(e) ~isempty(e.frame),extra.ex),1);
sample=struct();
if ~isempty(ie), sample.traj=extra.traj{ie}; sample.ex=extra.ex{ie}; sample.meta=extra.meta{ie}; sample.video=ie; end
save(fullfile(outdir,[name '_validation.mat']),'tracks','res','p','record','sample','-v7.3');
fid=fopen(fullfile(outdir,[name '_validation.json']),'w','n','UTF-8');
fprintf(fid,'%s',jsonencode(record,PrettyPrint=true)); fclose(fid);
end
