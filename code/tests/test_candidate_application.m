function metrics = test_candidate_application(outdir)
%TEST_CANDIDATE_APPLICATION 实际增益、分块校正与事后评价的一致性回归。
%   明确采用single增益；空间指标与直接小数组逐像素std比较。
if nargin<1, outdir=''; end
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root); setup_paths(root);
before_rng=rng; ny=7; nx=9; F=235; t=(0:F-1)'/1000;
video=reshape(uint16(mod(0:ny*nx*F-1,1500)+500),ny,nx,F);
global_gain=1+0.08*sin(2*pi*47.3*t+0.37);
global_model=struct('global_gain',global_gain,'spatial_coefficients',[], ...
    'temporal_basis',[],'is_spatial',false);
actual=apply_candidate_illumination(video,global_model);
baseline=single(video)./reshape(single(global_gain),1,1,[]);
assert(isequal(actual,baseline),'全局回退与原single归一化算序不同。');
idx=[1,101,F,3]; assert(isequal(candidate_gain_block(global_model,idx),reshape(single(global_gain(idx)),1,1,[])));
identity=global_model; identity.global_gain=ones(F,1);
assert(isequal(apply_candidate_illumination(video,identity),single(video)));
metrics=struct('scope','Application and post-inference metrics only; no sizing acceptance', ...
    'frames',F,'block_frames',100,'global_single_exact',true,'identity_exact',true);

% 两频率、四时间基，不同空间系数；显式矩阵是分块入口的独立参照。
[X,Y]=meshgrid(linspace(-1,1,nx),linspace(-1,1,ny));
basis=[sin(2*pi*47.3*t),cos(2*pi*47.3*t),sin(2*pi*83.1*t),cos(2*pi*83.1*t)];
coefficients=cat(3,0.04+0.01*X,0.02*Y,0.015*(X.*Y),0.01*(X.^2));
model=struct('global_gain',ones(F,1),'spatial_coefficients',coefficients, ...
    'temporal_basis',basis,'is_spatial',true);
explicit_gain=reshape(single(1+reshape(coefficients,ny*nx,4)*basis.'),ny,nx,F);
assert(isequal(candidate_gain_block(model,idx),explicit_gain(:,:,idx)));
assert(isequal(apply_candidate_illumination(video,model),single(video)./explicit_gain), ...
    '分块空间校正与显式矩阵不一致。');
metrics.spatial_single_matrix_exact=true; metrics.multifrequency_basis=true;

meta=struct('gain_mask',0.4+0.2*X+0.1*Y,'gain_offset',0.07*basis(:,1)+0.03*basis(:,4)+0.01);
measured=candidate_gain_metrics(model,meta);
direct=direct_metrics(double(explicit_gain),meta);
error=max(abs(metric_vector(measured)-metric_vector(direct)));
assert(error<1e-10,'分块二阶矩与直接逐像素std不一致。');
assert(isequal(fieldnames(measured),fieldnames(direct)));
fallback_metrics=candidate_gain_metrics(global_model,meta);
assert(isequal(fallback_metrics,targeted_gain_metrics(double(single(global_gain)),meta)));
metrics.maximum_metric_error=error; metrics.global_metric_exact=true;

% 理想模型仍实际转single，故相对double真值有浮点残差；近零容差明确给出。
perfect=model; perfect.spatial_coefficients=cat(3,0.07*meta.gain_mask,zeros(ny,nx), ...
    zeros(ny,nx),0.03*meta.gain_mask);
perfect_meta=meta; perfect_meta.gain_offset=0.07*basis(:,1)+0.03*basis(:,4);
near_zero=candidate_gain_metrics(perfect,perfect_meta);
direct_zero=direct_metrics(double(candidate_gain_block(perfect)),perfect_meta);
assert(near_zero.gain_rms<1e-7 && ...
    max(abs(metric_vector(near_zero)-metric_vector(direct_zero)))<1e-7);
metrics.near_perfect_rms=near_zero.gain_rms;
metrics.near_zero_numerical_tolerance=1e-7;
static_meta=meta; static_meta.gain_offset=zeros(F,1);
zero_model=model; zero_model.spatial_coefficients=zeros(size(coefficients));
static_metric=candidate_gain_metrics(zero_model,static_meta);
assert(all(metric_vector(static_metric)==0),'恒等空间模型的静态基线应严格为0。');
metrics.static_zero=true;

invalid_cases=0;
bad_models={struct(),rmfield(model,'temporal_basis')};
bad=global_model; bad.is_spatial=0; bad_models{end+1}=bad;
bad_gains={[],zeros(F,1),nan(F,1),inf(F,1),1i*ones(F,1),ones(5,47),realmax*ones(F,1)};
for k=1:numel(bad_gains), bad=global_model; bad.global_gain=bad_gains{k}; bad_models{end+1}=bad; end
bad=model; bad.global_gain(1)=1.1; bad_models{end+1}=bad;
bad=model; bad.temporal_basis=basis(1:end-1,:); bad_models{end+1}=bad;
bad=model; bad.temporal_basis(2)=nan; bad_models{end+1}=bad;
bad=model; bad.spatial_coefficients=coefficients(:,:,1:3); bad_models{end+1}=bad;
bad=model; bad.spatial_coefficients(1)=inf; bad_models{end+1}=bad;
bad=model; bad.spatial_coefficients=-2*ones(ny,nx,2); bad.temporal_basis=ones(F,2); bad_models{end+1}=bad;
bad=global_model; bad.spatial_coefficients=zeros(ny,nx,2); bad_models{end+1}=bad;
for k=1:numel(bad_models)
    verify_invalid(@() candidate_gain_block(bad_models{k},1)); invalid_cases=invalid_cases+1;
end
bad_indices={[],0,F+1,1.5,nan,1i,[1 2;3 4]};
for k=1:numel(bad_indices)
    verify_invalid(@() candidate_gain_block(model,bad_indices{k})); invalid_cases=invalid_cases+1;
end
bad_videos={[],nan(ny,nx,F),-ones(ny,nx,F),1i*ones(ny,nx,F), ...
    ones(ny,nx,F,2),zeros(ny,nx,F+1),zeros(ny+1,nx,F)};
for k=1:numel(bad_videos)
    verify_invalid(@() apply_candidate_illumination(bad_videos{k},model)); invalid_cases=invalid_cases+1;
end
bad_truths={struct(),rmfield(meta,'gain_offset')};
bad=meta; bad.gain_offset=zeros(F-1,1); bad_truths{end+1}=bad;
bad=meta; bad.gain_offset(1)=nan; bad_truths{end+1}=bad;
bad=meta; bad.gain_mask=zeros(ny+1,nx); bad_truths{end+1}=bad;
bad=meta; bad.gain_mask(1)=inf; bad_truths{end+1}=bad;
bad=meta; bad.gain_offset=-100*ones(F,1); bad_truths{end+1}=bad;
for k=1:numel(bad_truths)
    verify_invalid(@() candidate_gain_metrics(model,bad_truths{k})); invalid_cases=invalid_cases+1;
end
assert(isequal(before_rng,rng),'校正及指标计算不能改变随机流。');
metrics.invalid_inputs_rejected=invalid_cases; metrics.rng_unchanged=true; metrics.pass=true;
if ~isempty(outdir)
    if ~exist(outdir,'dir'), mkdir(outdir); end
    fid=fopen(fullfile(outdir,'candidate_application_tests.json'),'w','n','UTF-8');
    assert(fid>=0,'测试JSON无法写入。'); cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%s',jsonencode(metrics,PrettyPrint=true));
end
fprintf('PASS test_candidate_application: global/spatial single exact, metric error %.3g, %d invalid cases.\n',error,invalid_cases);
end

function metric=direct_metrics(gain,meta)
F=numel(meta.gain_offset); truth=1+meta.gain_mask.*reshape(meta.gain_offset,1,1,F);
ratio=truth./gain; residual=std(ratio./mean(ratio,3),1,3);
base=std(truth./mean(truth,3),1,3); nleft=floor(size(truth,2)/2);
metric=struct('gain_rms',rms_pixels(residual),'gain_rms_left',rms_pixels(residual(:,1:nleft)), ...
    'gain_rms_right',rms_pixels(residual(:,nleft+1:end)), ...
    'baseline_gain_rms',rms_pixels(base),'baseline_gain_rms_left',rms_pixels(base(:,1:nleft)), ...
    'baseline_gain_rms_right',rms_pixels(base(:,nleft+1:end)));
end

function value=rms_pixels(values)
value=sqrt(mean(values(:).^2));
end

function values=metric_vector(metric)
values=[metric.gain_rms,metric.gain_rms_left,metric.gain_rms_right, ...
    metric.baseline_gain_rms,metric.baseline_gain_rms_left,metric.baseline_gain_rms_right];
end

function verify_invalid(callback)
rejected=false;
try, callback(); catch err, rejected=startsWith(err.identifier,'iNTA:candidate:'); end
assert(rejected,'非法输入未被预期校验拒绝。');
end
