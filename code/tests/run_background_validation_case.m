function record = run_background_validation_case(name,outdir,seed)
%RUN_BACKGROUND_VALIDATION_CASE 同一原始视频上的背景融合配对压力验收。
% 冻结条件：曝光内±8%共同照明、423.7 Hz、初相位0.37 rad。
% 每段仅生成一次轨迹和含噪相机视频，median与spectral_common共享输入。
% 真照明仅在两分支完成推断后评价增益误差，不传入背景估计器。
% 默认种子沿用上一版基准，属于配对压力复测，不是新的独立保留数据。
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root); setup_paths(root);
p=apply_scenario(sim_params('full'),name);
if nargin<3
    if strcmpi(name,'mono30'), seed=10001; else, seed=20001; end
end
p.sim.seed=seed;
p.cam.flicker_amplitude=0.08; p.cam.flicker_frequency=423.7; p.cam.flicker_phase=0.37;
p.ana.background_mode='spectral_common';
assert(strcmpi(p.ana.fit_model,'perp'),'压力验收的主估计器固定为perp。');
p.out.dir=outdir; p.out.make_plots=false; p.out.save_video=false; p.out.save_diagnostics=false;
p_baseline=p; p_baseline.ana.background_mode='median';
if ~exist(outdir,'dir'), mkdir(outdir); end
result_file=fullfile(outdir,[name '_background_validation.mat']);
assert(~exist(result_file,'file'),'该组压力验收结果已存在，请使用新的输出目录。');
diary(fullfile(outdir,[name '_background_run.log'])); cleanup=onCleanup(@() diary('off')); %#ok<NASGU>
fprintf('%s：配对压力复测，seed=%d，%d段×%d帧；不称为全新保留样本。\n', ...
    name,seed,p.sim.n_videos,p.sim.n_frames);
total_started=tic; tab=ipsf_lookup_table(p); tracks=empty_tracks(); tracks_baseline=empty_tracks(); sample=struct();
template=struct('video',0,'seed',0,'frames',0,'applied',false,'accepted_frequencies',[], ...
    'amplitudes',[],'baseline_gain_rms',nan,'residual_gain_rms',nan,'suppression_fraction',nan, ...
    'saturated_pixels',0,'tracks_baseline',0,'tracks_fused',0, ...
    'selected_tracks_baseline',0,'selected_tracks_fused',0, ...
    'selected_localizations_baseline',0,'selected_localizations_fused',0);
video_records=repmat(template,1,p.sim.n_videos);
for v=1:p.sim.n_videos
    started=tic; current_seed=seed+100*(v-1);
    traj=simulate_trajectories(p,current_seed);
    [video,meta]=render_iscat_video(traj,p,tab);
    [tb,~]=inta_track_video(video,p_baseline);
    tb=analyze_tracks(tb,traj,p_baseline,v);
    [tn,ex]=inta_track_video(video,p);
    tn=analyze_tracks(tn,traj,p,v);
    tracks_baseline=[tracks_baseline,tb]; tracks=[tracks,tn]; %#ok<AGROW>
    if ~isfield(sample,'traj') && ~isempty(ex.frame)
        sample=struct('traj',traj,'ex',ex,'meta',meta,'video',v);
    end
    % 只在评价阶段读取渲染器真增益。除去常数比例，评价残留的共同起伏。
    truth_gain=meta.illumination_gain(:); truth_gain=truth_gain/mean(truth_gain);
    estimated_gain=ex.illumination.gain(:);
    assert(numel(estimated_gain)==numel(truth_gain) && all(isfinite(estimated_gain) & estimated_gain>0));
    baseline_ratio=truth_gain; corrected_ratio=truth_gain./estimated_gain;
    baseline_rms=std(baseline_ratio/mean(baseline_ratio),1);
    residual_rms=std(corrected_ratio/mean(corrected_ratio),1);
    % 本段样本量使用与最终组统计相同的入选口径；关闭重复bootstrap以节省计算。
    pn=p; pn.ana.bootstrap_n=0; pb=p_baseline; pb.ana.bootstrap_n=0;
    rn=inta_population_analysis(tn,pn); rb=inta_population_analysis(tb,pb);
    row=template; row.video=v; row.seed=current_seed; row.frames=size(video,3);
    row.applied=ex.illumination.applied;
    row.accepted_frequencies=ex.illumination.accepted_frequencies;
    row.amplitudes=ex.illumination.amplitudes;
    row.baseline_gain_rms=baseline_rms; row.residual_gain_rms=residual_rms;
    row.suppression_fraction=1-residual_rms/baseline_rms;
    row.saturated_pixels=meta.n_sat; row.tracks_baseline=numel(tb); row.tracks_fused=numel(tn);
    row.selected_tracks_baseline=rb.group.n_tracks; row.selected_tracks_fused=rn.group.n_tracks;
    row.selected_localizations_baseline=rb.group.n_loc; row.selected_localizations_fused=rn.group.n_loc;
    video_records(v)=row;
    fprintf('视频 %d/%d：应用=%d，频率=%s Hz，残余RMS=%.6g，合格轨迹旧/新=%d/%d，%.1fs\n', ...
        v,p.sim.n_videos,row.applied,mat2str(row.accepted_frequencies,6),residual_rms, ...
        row.selected_tracks_baseline,row.selected_tracks_fused,toc(started));
    clear video
end
res_baseline=inta_population_analysis(tracks_baseline,p_baseline);
res=inta_population_analysis(tracks,p); elapsed=toc(total_started);
record=struct('scenario',name,'seed',seed,'n_videos',p.sim.n_videos,'frames',p.sim.n_frames, ...
    'nominal_nm',p.sample.nominal_nm,'baseline',branch_record(tracks_baseline,res_baseline,p), ...
    'fused',branch_record(tracks,res,p),'applied_fraction',mean([video_records.applied]), ...
    'baseline_gain_rms',sqrt(mean([video_records.baseline_gain_rms].^2)), ...
    'residual_gain_rms',sqrt(mean([video_records.residual_gain_rms].^2)), ...
    'maximum_video_residual_gain_rms',max([video_records.residual_gain_rms]), ...
    'saturated_pixels',sum([video_records.saturated_pixels]),'seconds',elapsed);
record.suppression_fraction=1-record.residual_gain_rms/record.baseline_gain_rms;
record.evaluation_scope='上一版基准种子的配对动态照明压力复测；不是全新独立保留样本';
record.gain_metric='每视频std((真增益/估计增益)/mean(真增益/估计增益),1)，全部视频按帧数相同的均方合并；基线估计增益为1';
fprintf('\nBACKGROUND %s old=%.6f new=%.6f nm bias_new=%+.4f%% residual=%.6g suppression=%.3f%% applied=%.1f%%\n', ...
    name,record.baseline.diameter_nm,record.fused.diameter_nm,record.fused.bias_percent, ...
    record.residual_gain_rms,100*record.suppression_fraction,100*record.applied_fraction);
save(result_file,'tracks','res','tracks_baseline','res_baseline','p','p_baseline', ...
    'record','video_records','sample','-v7.3');
write_json(fullfile(outdir,[name '_background_validation.json']),struct('record',record,'videos',video_records));
end

function r=branch_record(tracks,res,p)
g=res.group; selected=find(res.sel & isfinite([tracks.loc_rms_nm])); loc_rms=nan;
if ~isempty(selected)
    w=[tracks(selected).n_loc];
    loc_rms=sqrt(sum(w.*[tracks(selected).loc_rms_nm].^2)/sum(w));
end
r=struct('diameter_nm',g.d_avgmsd,'D_um2_s',g.D_avgmsd, ...
    'bias_percent',100*(g.d_avgmsd/p.sample.nominal_nm-1),'ci_nm',g.d_ci95_interval, ...
    'n_tracks',g.n_tracks,'n_loc',g.n_loc,'videos_with_tracks',g.n_videos_with_tracks, ...
    'unmatched_selected',sum(res.sel & [tracks.truth_id]==0),'localization_rms_nm',loc_rms, ...
    'fit_model',g.fit_model);
end

function write_json(path,value)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0,'无法写入JSON。');
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(value,PrettyPrint=true));
end
