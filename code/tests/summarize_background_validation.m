function acceptance=summarize_background_validation(outdir)
%SUMMARIZE_BACKGROUND_VALIDATION 汇总同视频配对背景压力测试，按冻结阈值验收。
% 这些种子用于上一版基准回归，不把重复使用的样本称为新独立保留数据。
% 旧法仅作比较；即使旧法粒径已准确，也不据此夸大新法的粒径改善。
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root); setup_paths(root);
names={'mono30','mono100','mono200'}; seeds=[10001 20001 20001];
records=struct([]); protocols=struct([]); results=struct(); baseline_results=struct();
matches=false(1,3); differences=cell(1,3);
oldvis=get(groot,'defaultFigureVisible'); set(groot,'defaultFigureVisible','off');
cleanup=onCleanup(@() set(groot,'defaultFigureVisible',oldvis)); %#ok<NASGU>
for k=1:3
    S=load(fullfile(outdir,[names{k} '_background_validation.mat']));
    if k==1, records=S.record; else, records(k)=S.record; end %#ok<AGROW>
    expected=apply_scenario(sim_params('full'),names{k}); expected.sim.seed=seeds(k);
    expected.cam.flicker_amplitude=0.08; expected.cam.flicker_frequency=423.7; expected.cam.flicker_phase=0.37;
    expected.ana.background_mode='spectral_common';
    e=protocol_fields(expected); a=protocol_fields(S.p); fields=fieldnames(e); different={};
    for j=1:numel(fields)
        if ~isequaln(a.(fields{j}),e.(fields{j})), different{end+1}=fields{j}; end %#ok<AGROW>
    end
    eb=expected; eb.ana.background_mode='median';
    if ~isequaln(protocol_fields(S.p_baseline),protocol_fields(eb))
        different{end+1}='baseline_protocol'; %#ok<AGROW>
    end
    matches(k)=isempty(different); differences{k}=different;
    pr=struct('scenario',names{k},'expected',e,'actual',a);
    if k==1, protocols=pr; else, protocols(k)=pr; end %#ok<AGROW>
    results.(names{k})=struct('tracks',S.tracks,'res',S.res,'p',S.p);
    baseline_results.(names{k})=struct('tracks',S.tracks_baseline,'res',S.res_baseline,'p',S.p_baseline);
    pp=S.p; pp.out.dir=outdir; pp.out.make_plots=true;
    print_summary(S.tracks,S.res,pp,names{k});
    plot_results(S.tracks,S.res,pp,names{k});
    if isfield(S.sample,'traj')
        plot_trajectories_3d(S.sample.traj,pp,names{k});
        plot_example_frames(S.sample.ex,pp,names{k});
        plot_validation_diagnostics(S.sample,S.tracks,pp,names{k},S.res.sel);
        plot_illumination_diagnostics(S.sample,pp,names{k});
    end
    close all;
end
new=[records.fused]; baseline=[records.baseline]; %#ok<NASGU>
checks=struct();
checks.three_separate_monodisperse_groups=numel(records)==3 && isequal({records.scenario},names);
checks.fixed_regression_seeds=isequal([records.seed],seeds);
checks.fixed_video_and_frame_counts=isequal([records.n_videos],[20 60 60]) && all([records.frames]==5000);
checks.frozen_acquisition_and_analysis_protocol=all(matches);
checks.fixed_perpendicular_estimator=all(strcmp({new.fit_model},'perp'));
checks.finite_estimates=all(isfinite([new.diameter_nm]));
checks.bias_within_five_percent=all(abs([new.bias_percent])<=5);
checks.at_least_ten_observed_videos=all([new.videos_with_tracks]>=10);
checks.at_least_5000_localizations=all([new.n_loc]>=5000);
checks.finite_video_bootstrap_ci=all(isfinite([new.ci_nm]));
checks.low_unmatched_fraction=all([new.unmatched_selected]./max([new.n_tracks],1)<=0.02);
checks.no_adc_saturation=all([records.saturated_pixels]==0);
checks.applied_in_at_least_95_percent_videos=all([records.applied_fraction]>=0.95);
checks.residual_common_illumination_below_0_001=all(isfinite([records.residual_gain_rms]) & [records.residual_gain_rms]<0.001);
checks.at_least_95_percent_rms_suppression=all(isfinite([records.suppression_fraction]) & [records.suppression_fraction]>=0.95);
acceptance=struct('passed',all(structfun(@(x) logical(x),checks)),'checks',checks, ...
    'records',records,'protocols',protocols,'protocol_differences',{differences}, ...
    'scope','三组单分散、同一原始压力视频的配对比较；沿用上一版基准种子，不是全新独立保留样本；不证明真实设备或全球首创');
acceptance.thresholds=struct('maximum_absolute_diameter_bias_percent',5,'minimum_active_videos',10, ...
    'minimum_localizations',5000,'maximum_unmatched_fraction',0.02,'minimum_applied_fraction',0.95, ...
    'maximum_residual_gain_rms_exclusive',0.001,'minimum_rms_suppression_fraction',0.95);
fid=fopen(fullfile(outdir,'background_acceptance.json'),'w','n','UTF-8');
assert(fid>=0,'无法写入验收JSON。'); fprintf(fid,'%s',jsonencode(acceptance,PrettyPrint=true)); fclose(fid);
write_csv(records,fullfile(outdir,'background_comparison.csv'));
save(fullfile(outdir,'background_results.mat'),'results','baseline_results','acceptance','-v7.3');
% 图5仅展示融合分支；旧/新配对数值完整写入CSV与JSON，不选择更接近真值者。
p0=sim_params('full'); p0.out.dir=outdir; p0.out.make_plots=true;
plot_group_comparison(results,p0); close all;
disp(checks);
assert(acceptance.passed,'背景融合压力验收未通过，详见background_acceptance.json；不得标为完成。');
fprintf('PASS background acceptance: fixed paired regression stress protocol; all frozen thresholds passed.\n');
end

function s=protocol_fields(p)
% 比较全部物理、采集、光学和估计设置；仅输出路径/保存开关不影响配对判定。
s=struct('phys',p.phys,'part',p.part,'mag',p.mag,'sim',p.sim,'opt',p.opt,'cam',p.cam,'ana',p.ana);
% JSON不能直接编码复数；保留材料折射率的实部/虚部以完整核对光学协议。
s.part.n_p=struct('real',real(p.part.n_p),'imag',imag(p.part.n_p));
s.opt.n_ref=struct('real',real(p.opt.n_ref),'imag',imag(p.opt.n_ref));
end

function write_csv(records,path)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0,'无法写入背景对照CSV。');
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,['scenario,nominal_nm,baseline_diameter_nm,fused_diameter_nm,baseline_bias_percent,fused_bias_percent,' ...
    'fused_ci_low_nm,fused_ci_high_nm,baseline_tracks,fused_tracks,baseline_localizations,fused_localizations,' ...
    'baseline_active_videos,fused_active_videos,total_videos,applied_fraction,baseline_gain_rms,' ...
    'residual_gain_rms,rms_suppression_fraction,maximum_video_residual_gain_rms,saturated_pixels\n']);
for k=1:numel(records)
    r=records(k); b=r.baseline; n=r.fused;
    fprintf(fid,'%s,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%d,%d,%d,%d,%d,%d,%d,%.9g,%.9g,%.9g,%.9g,%.9g,%d\n', ...
        r.scenario,r.nominal_nm,b.diameter_nm,n.diameter_nm,b.bias_percent,n.bias_percent,n.ci_nm, ...
        b.n_tracks,n.n_tracks,b.n_loc,n.n_loc,b.videos_with_tracks,n.videos_with_tracks,r.n_videos, ...
        r.applied_fraction,r.baseline_gain_rms,r.residual_gain_rms,r.suppression_fraction, ...
        r.maximum_video_residual_gain_rms,r.saturated_pixels);
end
end
