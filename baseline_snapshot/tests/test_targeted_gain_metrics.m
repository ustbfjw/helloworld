function metrics = test_targeted_gain_metrics(outdir)
%TEST_TARGETED_GAIN_METRICS 照明指标的解析量、空间分区与输入校验回归。
%   只验证评价指标定义，不代表背景估计器或粒径反演已通过验收。
if nargin<1, outdir=''; end
F=1000; A=0.08; offset=A*sin(2*pi*(0:F-1)'/F);
meta=struct('gain_mask',ones(4,6),'gain_offset',offset);
names={'gain_rms','gain_rms_left','gain_rms_right','baseline_gain_rms', ...
    'baseline_gain_rms_left','baseline_gain_rms_right'};
uncorrected=targeted_gain_metrics(ones(F,1),meta);
assert(isequal(fieldnames(uncorrected),names'));
for k=1:6, assert(abs(uncorrected.(names{k})-A/sqrt(2))<2e-15); end
perfect=targeted_gain_metrics(1+offset,meta);
assert(perfect.gain_rms==0 && perfect.gain_rms_left==0 && perfect.gain_rms_right==0);
metrics=struct('scope','Illumination metric correctness only, not estimator or sizing acceptance', ...
    'global_uncorrected_rms',uncorrected.gain_rms,'global_perfect_correction_rms',perfect.gain_rms);

% 左半边扰动时，全图基线RMS为A/2，右半边真实照明恒定。
local=meta; local.gain_mask(:,4:6)=0;
baseline=targeted_gain_metrics(ones(F,1),local);
assert(abs(baseline.baseline_gain_rms-A/2)<2e-15 && ...
    abs(baseline.baseline_gain_rms_left-A/sqrt(2))<2e-15 && baseline.baseline_gain_rms_right==0);
assert(baseline.gain_rms==baseline.baseline_gain_rms);
% 若将空间均值的半幅扰动误当全局增益，原本恒定的右半边会被引入人工闪烁。
half_gain=1+offset/2;
local_result=targeted_gain_metrics(half_gain,local);
expected_right=sqrt(1/sqrt(1-(A/2)^2)-1);
assert(abs(local_result.gain_rms_right-expected_right)<2e-14);
assert(local_result.gain_rms_right>0.02 && local_result.baseline_gain_rms_right==0);
verify_direct(local_result,half_gain,local);
metrics.local_baseline_rms=baseline.baseline_gain_rms;
metrics.local_right_artificial_rms=local_result.gain_rms_right;

% 非二值、可重复掩膜值、负掩膜值及奇数列也按像素数正确加权。
arbitrary=struct('gain_mask',[0,0.2,0.2,0.7,1;1.3,-0.2,0.7,0.7,1], ...
    'gain_offset',offset+0.03);
arbitrary_gain=1+0.04*cos(2*pi*(0:F-1)'/F+0.4);
result=targeted_gain_metrics(arbitrary_gain,arbitrary);
verify_direct(result,arbitrary_gain,arbitrary);
metrics.nonbinary_weighted_matches_pixelwise=true;

% 估计增益的常量倍数不应改变相对波动；非零时间均值也须逐像素归一。
scaled=targeted_gain_metrics(7.3*arbitrary_gain,arbitrary);
for k=1:6, assert(abs(scaled.(names{k})-result.(names{k}))<2e-15); end
shifted=meta; shifted.gain_offset=offset+0.05;
shifted_result=targeted_gain_metrics(ones(F,1),shifted);
assert(abs(shifted_result.gain_rms-A/(sqrt(2)*1.05))<2e-15);
row_result=targeted_gain_metrics(arbitrary_gain',arbitrary);
for k=1:6, assert(row_result.(names{k})==result.(names{k})); end
metrics.constant_gain_scale_invariant=true;
metrics.temporal_mean_normalization=true;

% 静态及恒定空间照明的基线都必须严格为0；不输出抑制率字段。
static=arbitrary; static.gain_offset=zeros(F,1);
static_result=targeted_gain_metrics(ones(F,1),static);
for k=1:6, assert(static_result.(names{k})==0); end
static.gain_offset=0.1*ones(F,1);
constant_result=targeted_gain_metrics(3.7*ones(F,1),static);
for k=1:6, assert(constant_result.(names{k})==0); end
one_frame=struct('gain_mask',logical([1,0]),'gain_offset',0.1);
single_result=targeted_gain_metrics(2,one_frame);
for k=1:6, assert(single_result.(names{k})==0); end
metrics.static_baseline_strict_zero=true;

% 错误输入不得返回NaN指标或用不完整空间信息代替完整真增益。
bad_gain={[],zeros(F,1),-ones(F,1),nan(F,1),inf(F,1),1i*ones(F,1), ...
    ones(F-1,1),ones(2,F/2),'1'};
invalid_cases=0;
for k=1:numel(bad_gain)
    verify_invalid(bad_gain{k},meta); invalid_cases=invalid_cases+1;
end
bad_masks={[],nan(2),inf(2),1i*ones(2),ones(2,1),ones(2,2,2),'aa'};
for k=1:numel(bad_masks)
    bad=meta; bad.gain_mask=bad_masks{k}; verify_invalid(ones(F,1),bad); invalid_cases=invalid_cases+1;
end
bad_offsets={[],nan(F,1),inf(F,1),1i*ones(F,1),ones(F-1,1),ones(2,F/2),'a',-ones(F,1)};
for k=1:numel(bad_offsets)
    bad=meta; bad.gain_offset=bad_offsets{k}; verify_invalid(ones(F,1),bad); invalid_cases=invalid_cases+1;
end
verify_invalid(ones(F,1),rmfield(meta,'gain_mask')); invalid_cases=invalid_cases+1;
verify_invalid(ones(F,1),[meta,meta]); invalid_cases=invalid_cases+1;
verify_invalid(ones(F,1),0); invalid_cases=invalid_cases+1;
metrics.invalid_inputs_rejected=invalid_cases; metrics.pass=true;
if ~isempty(outdir)
    if ~exist(outdir,'dir'), mkdir(outdir); end
    fid=fopen(fullfile(outdir,'targeted_gain_metrics_tests.json'),'w','n','UTF-8');
    assert(fid>=0,'测试指标JSON无法写入。');
    close_file=onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%s',jsonencode(metrics,PrettyPrint=true));
end
fprintf('PASS test_targeted_gain_metrics: analytic/global/local/nonbinary/scale/static checks; %d invalid inputs.\n',invalid_cases);
end

function verify_direct(result,gain,meta)
% 小型逐像素独立参考计算，用于检验按unique值加权实现及左右分区。
residual=zeros(size(meta.gain_mask)); baseline=residual;
for i=1:numel(residual)
    truth=1+meta.gain_mask(i)*meta.gain_offset;
    ratio=truth./gain;
    residual(i)=std(ratio/mean(ratio),1);
    baseline(i)=std(truth/mean(truth),1);
end
nleft=floor(size(residual,2)/2);
left=residual(:,1:nleft); right=residual(:,nleft+1:end);
left0=baseline(:,1:nleft); right0=baseline(:,nleft+1:end);
expected=[sqrt(mean(residual(:).^2)),sqrt(mean(left(:).^2)),sqrt(mean(right(:).^2)), ...
    sqrt(mean(baseline(:).^2)),sqrt(mean(left0(:).^2)),sqrt(mean(right0(:).^2))];
actual=struct2array_local(result);
assert(max(abs(actual-expected))<2e-15,'unique掩膜加权结果与逐像素参考不一致。');
end

function values=struct2array_local(result)
% 不依赖工具箱中的struct2array。
values=[result.gain_rms,result.gain_rms_left,result.gain_rms_right, ...
    result.baseline_gain_rms,result.baseline_gain_rms_left,result.baseline_gain_rms_right];
end

function verify_invalid(gain,meta)
rejected=false;
try
    targeted_gain_metrics(gain,meta);
catch err
    rejected=strcmp(err.identifier,'targeted:InvalidGainInput');
end
assert(rejected,'非法输入未被预期校验拒绝。');
end
