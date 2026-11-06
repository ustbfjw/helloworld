function acceptance = summarize_acceptance(outdir)
%SUMMARIZE_ACCEPTANCE 对预先规定的三组完整图像链路验收，不更换主估计器。
% 此函数只汇总已经运行的run_validation_case，不替代运行本身。
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root); setup_paths(root);
names={'mono30','mono100','mono200'}; results=struct(); records=struct([]);
protocol_matches=false(1,3); protocol_differences=cell(1,3); protocols=struct([]);
p0=sim_params('full'); p0.out.dir=outdir;
oldvis=get(groot,'defaultFigureVisible'); set(groot,'defaultFigureVisible','off');
cleanup=onCleanup(@() set(groot,'defaultFigureVisible',oldvis)); %#ok<NASGU>
for k=1:3
    S=load(fullfile(outdir,[names{k} '_validation.mat']));
    if k==1, records=S.record; else, records(k)=S.record; end %#ok<AGROW>
    expected=protocol_fields(apply_scenario(sim_params('full'),names{k}));
    actual=protocol_fields(S.p); fields=fieldnames(expected); different={};
    for j=1:numel(fields)
        if ~isequaln(actual.(fields{j}),expected.(fields{j}))
            different{end+1}=fields{j}; %#ok<AGROW>
        end
    end
    protocol_matches(k)=isempty(different); protocol_differences{k}=different;
    protocol=struct('scenario',names{k},'expected',expected,'actual',actual);
    if k==1, protocols=protocol; else, protocols(k)=protocol; end %#ok<AGROW>
    if ~isempty(different)
        fprintf('协议不一致 %s：%s\n',names{k},strjoin(different,', '));
    end
    results.(names{k})=struct('tracks',S.tracks,'res',S.res,'p',S.p);
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
summary=plot_group_comparison(results,p0); %#ok<NASGU>
close all;
% 精度及最小证据量预先固定；不通过时保留失败记录，禁止按真值调结果。
checks=struct();
checks.three_independent_groups=numel(records)==3;
expected_protocols=[protocols.expected];
checks.fixed_regression_protocol=all([records.seed]==[10001 20001 20001] & ...
    [records.n_videos]==[expected_protocols.n_videos] & [records.frames]==[expected_protocols.n_frames]);
checks.frozen_acquisition_protocol=all(protocol_matches);
checks.finite_estimates=all(isfinite([records.diameter_nm]));
checks.bias_within_five_percent=all(abs([records.bias_percent])<=5);
checks.at_least_ten_observed_videos=all([records.videos_with_tracks]>=10);
checks.at_least_5000_localizations=all([records.n_loc]>=5000);
checks.finite_video_bootstrap_ci=all(isfinite([records.ci_nm]));
checks.low_unmatched_fraction=all([records.unmatched_selected]./max([records.n_tracks],1)<=0.02);
checks.no_adc_saturation=all([records.saturated_pixels]==0);
checks.no_correction_on_static_illumination=all([records.illumination_corrected_videos]==0);
acceptance=struct('passed',all(structfun(@(x) logical(x),checks)), ...
    'checks',checks,'records',records,'scope', ...
    '三组独立单分散、默认静态照明、上一版固定基准种子回归；不等于全新保留样本或真实装置校准');
acceptance.protocols=protocols;
acceptance.protocol_differences=protocol_differences;
fid=fopen(fullfile(outdir,'acceptance.json'),'w','n','UTF-8');
fprintf(fid,'%s',jsonencode(acceptance,PrettyPrint=true)); fclose(fid);
save(fullfile(outdir,'inta_results.mat'),'results','p0','acceptance','-v7.3');
disp(checks);
assert(acceptance.passed,'端到端验收未通过，详见acceptance.json；不得把失败结果标为完成。');
fprintf('PASS acceptance: all three independent groups within 5%%, fixed regression seeds, full image pipeline.\n');
end

function s=protocol_fields(p)
% 仅核对会改变采集或估计的关键设置；输出目录与保存开关允许不同。
s=struct('fit_model',p.ana.fit_model,'loc_win',p.ana.loc_win, ...
    'loc_iter',p.ana.loc_iter,'det_k',p.ana.det_k,'min_len',p.ana.min_len, ...
    'n_fit',p.ana.n_fit,'n_lag',p.ana.n_lag,'edge_nsig',p.ana.edge_nsig, ...
    'link_r',p.ana.link_r,'link_gap',p.ana.link_gap, ...
    'exposure',p.cam.exposure,'photons',p.cam.photons,'gain',p.cam.gain, ...
    'pixel_size',p.cam.px,'image_size',p.cam.Npx, ...
    'n_per_size',p.part.n_per_size,'n_sub',p.sim.n_sub,'fps',p.sim.fps, ...
    'n_videos',p.sim.n_videos,'n_frames',p.sim.n_frames, ...
    'L_xy',p.sim.L_xy,'H_z',p.sim.H_z,'z_init',p.sim.z_init, ...
    'G0',p.mag.G0,'f',p.mag.f,'wave',p.mag.wave,'phase0',p.mag.phase0, ...
    'dir',p.mag.dir(:).','B0',p.mag.B0, ...
    'background_mode',p.ana.background_mode,'illumination_options',p.ana.illum, ...
    'flicker_amplitude',p.cam.flicker_amplitude, ...
    'flicker_frequency',p.cam.flicker_frequency,'flicker_phase',p.cam.flicker_phase);
end
