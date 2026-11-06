function metric = targeted_gain_metrics(gain, meta)
%TARGETED_GAIN_METRICS 逐像素相对照明残差的时间标准差，再做空间均方根。
%   曝光平均真增益为true_gain(x,y,f)=1+gain_mask(x,y)*gain_offset(f)。
%   对每个像素计算ratio=true_gain/gain，将ratio除以自身时间均值，
%   再计算总体标准差（分母F）；最后在全图、左半图及右半图分别取RMS。
%   baseline使用gain=1，静态照明的baseline为0，不计算无定义的抑制百分比。
%   左半图为前floor(列数/2)列，剩余列归入右半图，与x<0掩膜约定一致。
%   按掩膜不同取值逐一计算，只保存F长向量，无需展开三维真值影像。
assert(isstruct(meta) && isscalar(meta) && all(isfield(meta,{'gain_mask','gain_offset'})), ...
    'targeted:InvalidGainInput','meta必须为含gain_mask和gain_offset的标量结构体。');
mask=meta.gain_mask; offset=meta.gain_offset;
assert(isnumeric(gain) && isreal(gain) && isvector(gain) && ~isempty(gain) && ...
    all(isfinite(gain(:))) && all(gain(:)>0), ...
    'targeted:InvalidGainInput','估计增益必须为有限正实数向量。');
assert((isnumeric(mask) || islogical(mask)) && isreal(mask) && ismatrix(mask) && ...
    ~isempty(mask) && size(mask,2)>=2 && all(isfinite(mask(:))), ...
    'targeted:InvalidGainInput','空间掩膜必须为有限实矩阵，且至少有两列。');
assert(isnumeric(offset) && isreal(offset) && isvector(offset) && ...
    numel(offset)==numel(gain) && all(isfinite(offset(:))), ...
    'targeted:InvalidGainInput','曝光平均偏移必须为与估计增益等长的有限实向量。');
gain=double(gain(:)); offset=double(offset(:));
[values,~,indices]=unique(double(mask(:)));
counts=accumarray(indices,1,[numel(values),1]);
indices=reshape(indices,size(mask)); nleft=floor(size(mask,2)/2);
left_indices=indices(:,1:nleft);
left_counts=accumarray(left_indices(:),1,[numel(values),1]);
right_counts=counts-left_counts;
residual_var=zeros(numel(values),1); baseline_var=zeros(numel(values),1);
for k=1:numel(values)
    true_gain=1+values(k)*offset;
    assert(all(isfinite(true_gain)) && all(true_gain>0), ...
        'targeted:InvalidGainInput','由掩膜与偏移构成的真照明必须为有限正强度。');
    ratio=true_gain./gain;
    assert(all(isfinite(ratio)) && all(ratio>0), ...
        'targeted:InvalidGainInput','真增益与估计增益比值超出有效数值范围。');
    residual_var(k)=relative_std(ratio)^2;
    baseline_var(k)=relative_std(true_gain)^2;
end
metric=struct('gain_rms',weighted_rms(residual_var,counts), ...
    'gain_rms_left',weighted_rms(residual_var,left_counts), ...
    'gain_rms_right',weighted_rms(residual_var,right_counts), ...
    'baseline_gain_rms',weighted_rms(baseline_var,counts), ...
    'baseline_gain_rms_left',weighted_rms(baseline_var,left_counts), ...
    'baseline_gain_rms_right',weighted_rms(baseline_var,right_counts));
end

function value=relative_std(series)
if all(series==series(1))
    value=0; % 恒定序列的理论标准差严格为零，避免常量求均值的舍入残差。
else
    % 先缩放至不超过1，避免有限正数求和溢出；不改变相对标准差。
    scaled=series/max(series);
    value=std(scaled/mean(scaled),1);
end
end

function value=weighted_rms(variances,weights)
value=sqrt(sum(variances.*weights)/sum(weights));
end
