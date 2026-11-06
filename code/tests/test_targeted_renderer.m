function metrics = test_targeted_renderer(outdir)
%TEST_TARGETED_RENDERER 共享渲染与生产渲染的逐位回归；不代表粒径精度验收。
%   3个单分散组，各100帧、完整128x128相机、4点有限曝光及光子/读出噪声。
%   所有global条件逐位比较视频、生产元数据字段和调用后的随机流状态。
%   额外检查局部照明左右半图、无噪声/瞬时模型、零调制及非法输入。
if nargin<1, outdir=''; end
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root); setup_paths(root);
saved_rng=rng; restore_rng=onCleanup(@() rng(saved_rng)); %#ok<NASGU>
timer=tic; p0=sim_params('full'); p0.sim.n_frames=100;
conditions=struct('name',{'static','global','near','outside','same','local_left_half'}, ...
    'amplitude',{0,0.08,0.08,0.08,0.08,0.08}, ...
    'frequency',{423.7,423.7,11.4,13.4,10,423.7}, ...
    'phase',{0.37,0.37,0.37,0.37,0.37,0.37}, ...
    'spatial_mode',{'global','global','global','global','global','left_half'});
tab=ipsf_lookup_table(p0);
metrics=struct('scope','Renderer equivalence only, not end-to-end sizing acceptance', ...
    'frames_per_group',100,'camera_Npx',p0.cam.Npx,'exposure_samples',p0.sim.n_sub, ...
    'shared_standard_normal_noise',true,'global_exact_cases',0);
groups={'mono30','mono100','mono200'};
for ig=1:3
    p=apply_scenario(p0,groups{ig});
    traj=simulate_trajectories(p,810000+ig);
    rng(910000+ig,'twister'); start_rng=rng;
    [videos,metas]=render_targeted_conditions(traj,p,tab,conditions); shared_rng=rng;
    assert(isequal(size(videos),[1,6]) && isequal(size(metas),[1,6]));
    for k=1:5
        pp=p; pp.cam.flicker_amplitude=conditions(k).amplitude;
        pp.cam.flicker_frequency=conditions(k).frequency; pp.cam.flicker_phase=conditions(k).phase;
        rng(start_rng); [reference,refmeta]=render_iscat_video(traj,pp,tab); reference_rng=rng;
        assert(isequal(reference,videos{k}),'targeted:ImageMismatch', ...
            '组%s条件%s未与生产渲染逐位一致。',groups{ig},conditions(k).name);
        names=fieldnames(refmeta);
        for j=1:numel(names)
            assert(isequaln(refmeta.(names{j}),metas{k}.(names{j})), ...
                'targeted:MetadataMismatch','生产元数据字段%s不一致。',names{j});
        end
        assert(isequal(reference_rng,shared_rng),'targeted:RngMismatch','随机流消耗不同。');
        metrics.global_exact_cases=metrics.global_exact_cases+1;
    end
    verify_local(videos,metas,traj,p,conditions);
    fprintf('TARGETED_RENDERER %s: 5 global videos/metadata/RNG exact; noisy local halves exact.\n',groups{ig});
end
metrics.global_video_metadata_rng_exact=true;
metrics.noisy_local_halves_exact=true;

% 额外无噪声回归：局部照明只能改变相机左半边，背景仍为同一次抽样。
p.cam.noise=false; p.sim.n_frames=12;
traj=simulate_trajectories(p,810010); rng(910010,'twister'); start_rng=rng;
[videos,metas]=render_targeted_conditions(traj,p,tab,conditions);
shared_rng=rng; verify_local(videos,metas,traj,p,conditions);
pp=p; pp.cam.flicker_amplitude=0; rng(start_rng);
reference=render_iscat_video(traj,pp,tab); reference_rng=rng;
assert(isequal(reference,videos{1}) && isequal(shared_rng,reference_rng));
metrics.noiseless_local_halves_exact=true;
metrics.noiseless_rng_exact=true;

% 零幅值时，全部空间/频率条件必须退化成相同静态视频和单位真增益。
zero=conditions;
for k=1:numel(zero), zero(k).amplitude=0; end
rng(start_rng); [zeros_video,zeros_meta]=render_targeted_conditions(traj,p,tab,zero);
for k=1:numel(zero)
    assert(isequal(zeros_video{k},reference) && all(zeros_meta{k}.gain_offset==0) && ...
        all(zeros_meta{k}.illumination_gain==1));
end
assert(isequal(rng,reference_rng));
metrics.zero_amplitude_identity=true;

% 兼容只有瞬时pos_wrap的外部轨迹；验证无t_exposure分支及nq=1。
instant=struct('d',traj.d,'n_frames',traj.n_frames,'pos_wrap',traj.pos_wrap);
rng(start_rng); [instant_video,instant_meta]=render_targeted_conditions(instant,p,tab,conditions(2));
shared_rng=rng; pp=p; pp.cam.flicker_amplitude=conditions(2).amplitude;
pp.cam.flicker_frequency=conditions(2).frequency; pp.cam.flicker_phase=conditions(2).phase;
rng(start_rng); [reference,refmeta]=render_iscat_video(instant,pp,tab);
assert(isequal(instant_video{1},reference) && isequal(shared_rng,rng) && ...
    isequal(instant_meta{1}.illumination_gain,refmeta.illumination_gain));
metrics.instantaneous_compatibility=true;

% 非法照明参数必须在抽随机数之前报错，不得静默生成全黑影像。
invalid_values={nan,inf,1i,[0 1],[],'0'};
keys={'amplitude','frequency','phase'}; invalid_cases=0;
for ik=1:numel(keys)
    for iv=1:numel(invalid_values)
        bad=conditions(1); bad.(keys{ik})=invalid_values{iv};
        verify_invalid(traj,p,tab,bad); invalid_cases=invalid_cases+1;
    end
end
for a=[-1,1,1.1]
    bad=conditions(1); bad.amplitude=a; verify_invalid(traj,p,tab,bad); invalid_cases=invalid_cases+1;
end
bad=conditions(1); bad.frequency=-1; verify_invalid(traj,p,tab,bad); invalid_cases=invalid_cases+1;
bad=conditions(1); bad.spatial_mode='right_half'; verify_invalid(traj,p,tab,bad); invalid_cases=invalid_cases+1;
bad=conditions(1); bad.name=''; verify_invalid(traj,p,tab,bad); invalid_cases=invalid_cases+1;
verify_invalid(traj,p,tab,struct([])); invalid_cases=invalid_cases+1;
verify_invalid(traj,p,tab,rmfield(conditions(1),'frequency')); invalid_cases=invalid_cases+1;
metrics.invalid_conditions_rejected_before_rng=invalid_cases;
metrics.pass=true; metrics.elapsed_seconds=toc(timer);
if ~isempty(outdir)
    if ~exist(outdir,'dir'), mkdir(outdir); end
    fid=fopen(fullfile(outdir,'targeted_renderer_tests.json'),'w','n','UTF-8');
    assert(fid>=0,'测试指标JSON无法写入。');
    close_file=onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%s',jsonencode(metrics,PrettyPrint=true));
end
fprintf('PASS test_targeted_renderer: %d exact global cases; local/zero/instantaneous checks; %d invalid inputs (%.2f s).\n', ...
    metrics.global_exact_cases,invalid_cases,metrics.elapsed_seconds);
end

function verify_local(videos,metas,traj,p,conditions)
N=p.cam.Npx; left=1:floor(N/2); right=floor(N/2)+1:N;
assert(isequal(videos{6}(:,left,:),videos{2}(:,left,:)), ...
    'targeted:LocalLeftMismatch','局部照明左半图必须等于全局同幅照明。');
assert(isequal(videos{6}(:,right,:),videos{1}(:,right,:)), ...
    'targeted:LocalRightMismatch','局部照明右半图必须等于静态照明。');
mask=zeros(N); mask(:,left)=1;
assert(isequal(metas{6}.gain_mask,mask) && isequal(metas{2}.gain_mask,ones(N)));
truth=mean(conditions(6).amplitude*sin(2*pi*conditions(6).frequency*traj.t_exposure+conditions(6).phase),2);
assert(max(abs(metas{6}.gain_offset-truth))<5e-16);
assert(max(abs(metas{6}.illumination_gain-(1+mean(mask(:))*truth)))<5e-16);
assert(max(abs(1+metas{2}.gain_offset-metas{2}.illumination_gain))<5e-16);
assert(isequal(metas{6}.bg_static,metas{1}.bg_static));
end

function verify_invalid(traj,p,tab,conditions)
before=rng; rejected=false;
try
    render_targeted_conditions(traj,p,tab,conditions);
catch err
    rejected=strcmp(err.identifier,'targeted:InvalidConditions');
end
assert(rejected,'非法条件未被预期输入检查拒绝。');
assert(isequal(before,rng),'非法输入不得污染随机流。');
end
