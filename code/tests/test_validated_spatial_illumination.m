function metrics=test_validated_spatial_illumination(outdir)
%TEST_VALIDATED_SPATIAL_ILLUMINATION 局部增益候选的合成机制测试。
% 不替代真实iSCAT渲染、盲定位/链接/粒径和独立保留条件验收。
if nargin<1, outdir=''; end
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root); setup_paths(root);
p=sim_params('full'); N=64; F=5000; t=(0:F-1)'/p.sim.fps;
[X,Y]=meshgrid(linspace(-1,1,N));
B=800*(0.7+0.3*exp(-2*(X.^2+Y.^2)));
rs=RandStream('mt19937ar','Seed',202611051);
clean=repmat(single(B),1,1,F)+single(0.2*randn(rs,N,N,F));
before=rng;
[static_model,static_info]=estimate_validated_spatial_illumination(clean,p);
assert_fallback(static_model,clean,p,static_info);
assert(isequal(static_model.global_gain,ones(F,1)));
partition=static_info.partition;
for k=1:16
    assert(isempty(intersect(partition.training_linear_indices{k},partition.validation_linear_indices{k})));
    assert(all(partition.valid_pixel_counts(k,:)>=4));
end
assert(numel(partition.reference_frames)<=256);

wave=sin(2*pi*423.7*t+0.37);
global_video=clean.*reshape(single(1+0.08*wave),1,1,F);
[global_model,global_info]=estimate_validated_spatial_illumination(global_video,p);
assert_fallback(global_model,global_video,p,global_info);

% 有明确二阶空间幅度变化的场；形状不作为估计器输入。
amplitude=0.04+0.035*X+0.025*Y+0.012*X.^2+0.006*X.*Y+0.008*Y.^2;
video=clean.*(1+single(amplitude).*reshape(single(wave),1,1,F));
[smooth_model,smooth_info]=estimate_validated_spatial_illumination(video,p);
assert(smooth_model.is_spatial && smooth_info.applied,'二阶空间输入未通过局部模型验证。');
assert(isequal(smooth_info.local_frequencies,smooth_info.accepted_frequencies));
assert(isequal(smooth_model.global_gain,ones(F,1)),'空间模型不得叠加原全局增益。');
assert(max(abs(mean(smooth_model.temporal_basis,1)))<1e-13);
assert(isscalar(smooth_info.accepted_frequencies) && abs(smooth_info.accepted_frequencies-423.7)<0.02);
assert(smooth_info.score.relative_energy_reduction>=0.5);
assert(smooth_info.score.global_residual_amplitude_rms>=0.002);
assert(contains(smooth_info.score.amplitude_definition,'not time-domain RMS'));
assert(all(smooth_info.score.local_residual_amplitude<=smooth_info.score.global_residual_amplitude+0.002));
assert(all(smooth_info.score.local_residual_amplitude<=smooth_info.score.uncorrected_amplitude+0.002));
assert(smooth_info.gain_range(1)>0.5 && smooth_info.gain_range(2)<=1.5);
identity=struct('global_gain',ones(F,1),'spatial_coefficients',[],'temporal_basis',[],'is_spatial',false);
uncorrected_rms=gain_error(identity,amplitude,wave);
corrected_rms=gain_error(smooth_model,amplitude,wave);
assert(corrected_rms<0.001 && corrected_rms<0.1*uncorrected_rms, ...
    '二阶场的后验照明误差未充分降低。');
fprintf('SPATIAL smooth: before=%.8g after=%.8g validation_reduction=%.8g f=%.8f\n', ...
    uncorrected_rms,corrected_rms,smooth_info.score.relative_energy_reduction,smooth_info.accepted_frequencies);

% 多频接口必须联合生成2K列系数/基函数，不能重复叠除各频率增益。
amplitude2=0.02-0.025*X+0.012*Y; wave2=sin(2*pi*731.2*t-0.25);
multivideo=clean.*(1+single(amplitude).*reshape(single(wave),1,1,F)+ ...
    single(amplitude2).*reshape(single(wave2),1,1,F));
[multi_model,multi_info]=estimate_validated_spatial_illumination(multivideo,p);
assert(multi_model.is_spatial && numel(multi_info.accepted_frequencies)==2);
assert(size(multi_model.spatial_coefficients,3)==4 && size(multi_model.temporal_basis,2)==4);
multi_rms=gain_error(multi_model,cat(3,amplitude,amplitude2),[wave wave2]);
assert(multi_rms<0.001);

% 只改留出像素：训练频率、训练幅相和多项式系数必须完全不变，但必须拒绝。
bad=reshape(video,N*N,F); original=reshape(clean,N*N,F);
anti_amplitude=0.08-amplitude(:);
for k=1:16
    ids=partition.validation_linear_indices{k};
    bad(ids,:)=original(ids,:).*single(1+anti_amplitude(ids)*wave.');
end
bad=reshape(bad,N,N,F);
[bad_model,bad_info]=estimate_validated_spatial_illumination(bad,p);
assert_fallback(bad_model,bad,p,bad_info);
assert(isequal(smooth_info.local_frequencies,bad_info.local_frequencies));
assert(isequal(smooth_info.training_coefficients,bad_info.training_coefficients));
assert(isequal(smooth_info.polynomial_coefficients,bad_info.polynomial_coefficients));
assert(bad_info.score.local_energy>0.5*bad_info.score.global_energy || ...
    any(bad_info.score.local_residual_amplitude>bad_info.score.global_residual_amplitude+0.002));

% 磁同频与邻近保护区都必须保留原结果，不允许“空间验证”解除频带保护。
protected_frequencies=[10 11.4];
for k=1:numel(protected_frequencies)
    magnetic_video=clean.*(1+single(amplitude).* ...
        reshape(single(sin(2*pi*protected_frequencies(k)*t+0.37)),1,1,F));
    [magnetic_model,magnetic_info]=estimate_validated_spatial_illumination(magnetic_video,p);
    assert_fallback(magnetic_model,magnetic_video,p,magnetic_info);
    assert(isempty(magnetic_info.accepted_frequencies));
    assert(any(contains(magnetic_info.rejection_reasons,'磁驱动保护区')));
end
% 辅助覆盖细化后进入保护区的路径：13Hz粗峰通过，12.55Hz细峰被保护。
pp=p; pp.mag.f=10.6;
refined_video=clean.*(1+single(amplitude).*reshape(single(sin(2*pi*12.55*t+0.37)),1,1,F));
[refined_model,refined_info]=estimate_validated_spatial_illumination(refined_video,pp);
assert_fallback(refined_model,refined_video,pp,refined_info);
assert(any(contains(refined_info.rejection_reasons,'细化频率位于磁驱动保护区')));

% 生成器、粒径、真值字段不是估计输入。
pt=p; pt.cam.flicker_amplitude=0.99; pt.cam.flicker_frequency=777; pt.cam.flicker_phase=-3;
pt.cam.flicker_spatial_shape='unrelated'; pt.part.d_list=[1 2 3]*1e-9;
pt.truth=struct('illumination_gain',zeros(F,1),'positions',zeros(F,2));
[same_model,same_info]=estimate_validated_spatial_illumination(video,pt);
assert(isequaln(smooth_model,same_model) && isequaln(smooth_info,same_info));
assert(isequal(before,rng),'候选改变了全局仿真随机流。');

% 大幅空间场即使通过相对评分，仍必须满足全图全时正增益约束。
large_amplitude=0.35+0.25*X+0.15*Y;
large_video=clean.*(1+single(large_amplitude).*reshape(single(wave),1,1,F));
[large_model,large_info]=estimate_validated_spatial_illumination(large_video,p);
assert_fallback(large_model,large_video,p,large_info);
assert(contains(large_info.reason,'安全正值范围'));

zero_video=zeros(32,32,32,'single');
[zero_model,zero_info]=estimate_validated_spatial_illumination(zero_video,p);
assert_fallback(zero_model,zero_video,p,zero_info);
assert(contains(zero_info.reason,'支持不足'));

invalid_count=0;
invalid_inputs={ones(32),ones(16,16,32),complex(ones(32,32,32)), ...
    -ones(32,32,32),nan(32,32,32),inf(32,32,32)};
for k=1:numel(invalid_inputs)
    must_fail(@() estimate_validated_spatial_illumination(invalid_inputs{k},p));
    invalid_count=invalid_count+1;
end
invalid_parameters={};
pp=p; pp.sim.fps=NaN; invalid_parameters{end+1}=pp;
pp=p; pp.sim.fps=[5000 5000]; invalid_parameters{end+1}=pp;
pp=p; pp.mag.f=-1; invalid_parameters{end+1}=pp;
pp=p; pp=rmfield(pp,'mag'); invalid_parameters{end+1}=pp;
pp=p; pp.ana.illum.enable_magnetic_guard=false; invalid_parameters{end+1}=pp;
pp=p; pp.ana.illum.protect_harmonics=0; invalid_parameters{end+1}=pp;
pp=p; pp.ana.illum.protect_bins=0; invalid_parameters{end+1}=pp;
pp=p; pp.ana.illum.min_amplitude=NaN; invalid_parameters{end+1}=pp;
pp=p; pp.ana.illum.reference_frames=0; invalid_parameters{end+1}=pp;
for k=1:numel(invalid_parameters)
    must_fail(@() estimate_validated_spatial_illumination(ones(32,32,32,'single'),invalid_parameters{k}));
    invalid_count=invalid_count+1;
end
metrics=struct('scope','Synthetic candidate mechanism tests; no end-to-end sizing claim', ...
    'matlab_version',version,'seed',202611051,'n_frames',F,'fps',p.sim.fps, ...
    'static_exact_fallback',true,'global_exact_fallback',true,'actual_correction_metadata',true, ...
    'smooth_field_accepted',true, ...
    'uncorrected_gain_rms',uncorrected_rms,'corrected_gain_rms',corrected_rms, ...
    'multiple_frequency_gain_rms',multi_rms,'multiple_frequency_interface',true, ...
    'validation_energy_reduction',smooth_info.score.relative_energy_reduction, ...
    'global_validation_residual_amplitude_rms',smooth_info.score.global_residual_amplitude_rms, ...
    'score_amplitude_definition',smooth_info.score.amplitude_definition, ...
    'heldout_damage_rejected',true,'training_unchanged_after_heldout_change',true, ...
    'magnetic_guard_retained',true,'refined_magnetic_guard',true,'rng_unchanged',true,'no_truth_leakage',true, ...
    'gain_range_guard',true,'invalid_cases_rejected',invalid_count,'pass',true, ...
    'static_reason',static_info.reason,'global_reason',global_info.reason, ...
    'heldout_reason',bad_info.reason,'large_gain_reason',large_info.reason);
if ~isempty(outdir)
    if ~exist(outdir,'dir'), mkdir(outdir); end
    fid=fopen(fullfile(outdir,'validated_spatial_illumination_tests.json'),'w','n','UTF-8');
    assert(fid>=0); cleanup=onCleanup(@() fclose(fid));
    fprintf(fid,'%s',jsonencode(metrics,PrettyPrint=true));
    save(fullfile(outdir,'validated_spatial_illumination_tests.mat'),'metrics','smooth_info','bad_info','large_info');
end
fprintf('PASS test_validated_spatial_illumination: fallback, spatial model, heldout, magnetic guard and inputs.\n');
end

function assert_fallback(model,video,p,info)
[gain,fallback_info]=estimate_common_illumination(video,p);
assert(~model.is_spatial && isequal(model.global_gain,gain));
assert(isempty(model.spatial_coefficients) && isempty(model.temporal_basis));
assert(~info.is_spatial && info.applied==fallback_info.applied);
assert(isequal(info.accepted_frequencies,fallback_info.accepted_frequencies));
end

function rms=gain_error(model,amplitude,wave)
amplitude=reshape(amplitude,[],size(wave,2)); P=size(amplitude,1); F=size(wave,1);
total=zeros(P,1); squares=zeros(P,1);
if model.is_spatial, coefficients=reshape(model.spatial_coefficients,P,[]); end
for start=1:128:F
    idx=start:min(F,start+127);
    if model.is_spatial
        gain=1+coefficients*model.temporal_basis(idx,:).';
    else
        gain=model.global_gain(idx).';
    end
    ratio=(1+amplitude*wave(idx,:).')./gain;
    total=total+sum(ratio,2); squares=squares+sum(ratio.^2,2);
end
variance=max(squares/F./(total/F).^2-1,0); rms=sqrt(mean(variance));
end

function must_fail(fun)
failed=false;
try
    fun();
catch
    failed=true;
end
assert(failed,'非法输入未被拒绝。');
end
