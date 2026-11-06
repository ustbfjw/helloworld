function record=test_candidate_summary(outdir)
%TEST_CANDIDATE_SUMMARY 140个合成统计MAT的端到端汇总测试，不生成仿真影片。
% 此fixture用于检查汇总代码，不是性能证据，不得并入正式实验结果。
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root); setup_paths(root);
if nargin<1, outdir=fullfile(fileparts(root),'candidate_statistics_fixture'); end
if ~exist(outdir,'dir'), mkdir(outdir); end
protocol=candidate_experiment_protocol(); protocol.synthetic_statistics_fixture=true;
write_text(fullfile(outdir,'protocol.json'),jsonencode(protocol,'PrettyPrint',true));
protocol_hash=targeted_sha256(fullfile(outdir,'protocol.json'));
write_text(fullfile(outdir,'protocol.sha256'),protocol_hash);
for ig=1:3
    group=protocol.groups(ig); folder=fullfile(outdir,group.name);
    if ~exist(folder,'dir'), mkdir(folder); end
    p=apply_scenario(sim_params('full'),group.name); p.ana.bootstrap_n=0;
    p.sim.n_videos=group.n_videos;
    empty=targeted_sufficient_stats(empty_tracks(),struct('sel',false(1,0)),p);
    tau=(1:p.ana.n_lag)'/p.sim.fps;
    D=p.phys.kB*p.phys.T/(3*pi*p.phys.eta*group.diameter_nm*1e-9)*1e12;
    for video_index=1:group.n_videos
        seed=group.seed_base+protocol.seed_stride*(video_index-1);
        s=empty; s.counts(:)=500;
        Dv=D*(1+0.03*sin(video_index));
        M=[4*Dv*tau+0.001,2*Dv*tau+0.0005,tau.^2,tau.^2,zeros(size(tau)),zeros(size(tau))];
        s.moments=M.*s.counts; s.driftA=diag([1 .1]); s.driftB=[2;.9];
        s.n_tracks=1; s.n_loc=600; s.matched_localizations=600;
        s.localization_squared_sum_nm2=600*2^2;
        metric=struct('gain_rms',0,'gain_rms_left',0,'gain_rms_right',0, ...
            'baseline_gain_rms',0,'baseline_gain_rms_left',0,'baseline_gain_rms_right',0, ...
            'applied',false,'saturated_pixels',0,'truth_velocity_um_s',9);
        items=repmat(struct('stats',s,'metric',metric,'illumination_info',struct(), ...
            'reused_from',''),6,4);
        for ic=2:6
            for im=1:4
                gain=.004; if im==4, gain=.002; end
                items(ic,im).metric.gain_rms=gain;
                items(ic,im).metric.gain_rms_left=gain;
                items(ic,im).metric.gain_rms_right=gain;
                items(ic,im).metric.baseline_gain_rms=.01;
                items(ic,im).metric.baseline_gain_rms_left=.01;
                items(ic,im).metric.baseline_gain_rms_right=.01;
                items(ic,im).metric.applied=true;
            end
        end
        for ic=1:6
            accepted=ismember(ic,[4 5]);
            items(ic,4).illumination_info=struct('is_spatial',accepted,'applied',ic>1);
        end
        % 保持有效位移对，但构造非正D，确认汇总如实保留而非判定代码失败。
        items(3,4).stats.moments(:,2)=(1-2*Dv*tau).*s.counts;
        % 仅一个视频有有效定位：不能伪造n=20/60的粒径精度证据。
        if video_index>1, items(4,4).stats=empty; end
        % 波形与恒流不可辨识，速度必须为NaN且不生成置信区间。
        items(5,4).stats.driftA=ones(2); items(5,4).stats.driftB=[2;2];
        % 完全无检测的视频仍保留在总体和配对索引中。
        if video_index==1, items(6,4).stats=empty; end
        save(fullfile(folder,sprintf('video_%03d.mat',video_index)), ...
            'items','p','video_index','seed','protocol_hash','-v7');
    end
end
state=rng; tic; report=summarize_candidate_experiment(outdir); elapsed=toc;
assert(isequal(state,rng),'汇总不能改变全局随机数状态。');
assert(numel(report.summary)==72 && numel(report.pairwise)==54);
assert(sum([report.pairwise.primary_comparison])==6);
s=report.summary;
pick=@(condition,method) strcmp({s.condition},condition) & strcmp({s.method},method);
r=s(pick('gaussian_patch','spatial_validated'));
assert(all([r.active_videos]==1) && all(isnan([r.d_ci95_low])));
assert(isequal([r.n_videos],[20 60 60]));
r=s(pick('same_frequency','spatial_validated'));
assert(all([r.D_um2_s]<0) && all(isnan([r.diameter_nm])) && all([r.bootstrap_valid]==0));
r=s(pick('tilted_smooth','spatial_validated'));
assert(~any([r.waveform_identifiable]) && all(isnan([r.magnetic_velocity_um_s])) && ...
    all(isnan([r.magnetic_velocity_ci95_low])));
r=s(pick('oblique_step','spatial_validated'));
assert(isequal([r.active_videos],[19 59 59]));
r=s(pick('static','spatial_validated')); assert(all(isnan([r.gain_rms_reduction_pct])));
r=s(pick('global_offbin','spatial_validated'));
assert(all(abs([r.gain_rms_pct]-.2)<1e-12) && all([r.waveform_identifiable]));
c=report.pairwise; r=c(strcmp({c.condition},'static'));
assert(all(strcmp({c.method},'spatial_validated')) && isfield(c,'active_candidate') && ~isfield(c,'active_full'));
assert(all([r.delta_D_um2_s]==0) && all([r.delta_D_ci95_low]==0) && ...
    all([r.delta_abs_bias_upper95]==0));
r=c([c.primary_comparison]); assert(all(abs([r.delta_gain_rms_pp]+.2)<1e-12));
assert(all(strcmp({r.reference},'full')) && ...
    all(ismember({r.condition},{'gaussian_patch','tilted_smooth'})));
r=s(pick('gaussian_patch','spatial_validated'));
assert(isequal([r.spatial_accepted_videos],[20 60 60]) && isequal([r.spatial_decision_videos],[20 60 60]));
r=s(pick('global_offbin','spatial_validated'));
assert(all([r.spatial_accepted_videos]==0),'不能将全局校正applied误读为空间模型接受。');
raw=readtable(fullfile(outdir,'video_metrics.csv')); assert(height(raw)==3360);
assert(sum(raw.spatial_accepted==1)==280);
assert(sum(raw.spatial_decision_available)==840);
assert(~any(raw.spatial_accepted(~raw.spatial_decision_available)));
% 错误seed须在计算开始前拒绝，测试后还原fixture，不能静默混用重复数据。
badfile=fullfile(outdir,'mono30','video_001.mat'); original=load(badfile,'seed');
seed=original.seed+1; save(badfile,'seed','-append'); rejected=false;
try, summarize_candidate_experiment(outdir); catch, rejected=true; end
seed=original.seed; save(badfile,'seed','-append');
assert(rejected,'不一致的seed必须导致完整性检查失败。');
% 保留旧hash，只改变仍合法的JSON字节，公开汇总入口也必须立即拒绝。
protocol_path=fullfile(outdir,'protocol.json'); original_text=fileread(protocol_path);
write_text(protocol_path,[original_text,sprintf('\n')]); protocol_rejected=false;
try
    summarize_candidate_experiment(outdir);
catch error_info
    protocol_rejected=strcmp(error_info.identifier,'candidate:ProtocolHashMismatch');
end
write_text(protocol_path,original_text);
assert(protocol_rejected,'协议JSON改变而hash未变必须由SHA256检查拒绝。');
assert(strcmp(targeted_sha256(protocol_path),protocol_hash),'测试后必须完整还原协议。');
record=struct('test','candidate_summary_synthetic_fixture','passed',true, ...
    'input_mat_files',140,'summary_rows',72,'pairwise_rows',54,'video_metric_rows',3360, ...
    'bootstrap_replicates',2000,'elapsed_seconds',elapsed, ...
    'empty_videos_retained',true,'one_active_video_has_no_size_ci',true, ...
    'negative_D_retained',true,'unidentifiable_velocity_is_nan',true, ...
    'seed_mismatch_rejected',rejected,'spatial_decision_separated_from_global_applied',true, ...
    'protocol_content_hash_mismatch_rejected',protocol_rejected, ...
    'synthetic_fixture_not_performance_evidence',true);
write_text(fullfile(outdir,'test_record.json'),jsonencode(record,'PrettyPrint',true));
fprintf('test_candidate_summary PASS: 140 fixture MAT, 72 summary rows, 54 comparisons; %.1f s.\n',elapsed);
end

function write_text(path,text)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0);
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',text);
end
