function metric = candidate_gain_metrics(model, meta)
%CANDIDATE_GAIN_METRICS 推断后的照明误差评价，使用与校正完全相同的single增益。
%   每像素r=true_gain/actual_gain，除以自身时间均值后求总体标准差；
%   再在全图、左半图、右半图分别求像素RMS。只返回六个约定指标。
%   空间模式按100帧累计sum(r)、sum(r^2)，无需整段double真值视频。
%   二阶矩相减存在取消误差：负方差钳0，近完美结果有约sqrt(eps)
%   的数值地板，不能把该精度下的数值排序解释为技术效果。
candidate_gain_block(model,1);
if ~model.is_spatial
    metric=targeted_gain_metrics(double(single(model.global_gain)),meta); return;
end
assert(isstruct(meta) && isscalar(meta) && all(isfield(meta,{'gain_mask','gain_offset'})), ...
    'iNTA:candidate:InvalidTruth','meta必须含gain_mask和gain_offset。');
mask=meta.gain_mask; offset=meta.gain_offset; F=numel(model.global_gain);
ny=size(model.spatial_coefficients,1); nx=size(model.spatial_coefficients,2);
assert((isnumeric(mask) || islogical(mask)) && isreal(mask) && ismatrix(mask) && ...
    isequal(size(mask),[ny,nx]) && nx>=2 && all(isfinite(mask(:))), ...
    'iNTA:candidate:InvalidTruth','真值掩膜必须为匹配空间模型的有限实矩阵，且至少两列。');
assert(isnumeric(offset) && isreal(offset) && isvector(offset) && numel(offset)==F && ...
    all(isfinite(offset(:))),'iNTA:candidate:InvalidTruth','曝光平均偏移必须为F长有限实向量。');
mask=double(mask(:)); offset=double(offset(:)); P=ny*nx;
sum_ratio=zeros(P,1); sum_squared=zeros(P,1);
for first=1:100:F
    idx=first:min(first+99,F);
    truth=1+mask*offset(idx).';
    assert(all(isfinite(truth(:))) && all(truth(:)>0), ...
        'iNTA:candidate:InvalidTruth','真照明强度必须有限且为正。');
    actual=reshape(double(candidate_gain_block(model,idx)),P,[]);
    ratio=truth./actual;
    assert(all(isfinite(ratio(:))) && all(ratio(:)>0), ...
        'iNTA:candidate:InvalidTruth','真照明与实际增益的比值无效。');
    sum_ratio=sum_ratio+sum(ratio,2); sum_squared=sum_squared+sum(ratio.^2,2);
end
assert(all(isfinite(sum_ratio)) && all(isfinite(sum_squared)) && all(sum_ratio>0), ...
    'iNTA:candidate:InvalidTruth','照明误差累计超出数值范围。');
mean_ratio=sum_ratio/F;
variance=max(0,sum_squared/F-mean_ratio.^2);
residual_variance=reshape(variance./mean_ratio.^2,ny,nx);
% 真值具有1+mask*offset形式，基线方差可直接求解，静态offset严格为0。
offset_std=std(offset,1); if all(offset==offset(1)), offset_std=0; end
baseline_variance=reshape((mask*offset_std./(1+mask*mean(offset))).^2,ny,nx);
nleft=floor(nx/2);
metric=struct('gain_rms',spatial_rms(residual_variance), ...
    'gain_rms_left',spatial_rms(residual_variance(:,1:nleft)), ...
    'gain_rms_right',spatial_rms(residual_variance(:,nleft+1:end)), ...
    'baseline_gain_rms',spatial_rms(baseline_variance), ...
    'baseline_gain_rms_left',spatial_rms(baseline_variance(:,1:nleft)), ...
    'baseline_gain_rms_right',spatial_rms(baseline_variance(:,nleft+1:end)));
end

function value=spatial_rms(variances)
value=sqrt(mean(variances(:)));
end
