function corrected = apply_candidate_illumination(video, model)
%APPLY_CANDIDATE_ILLUMINATION 按100帧块对原始视频施加实际单精度增益。
%   输入只能为相机视频和推断模型，不接受照明真值或真实轨迹。
%   生产追踪入口可直接逐块调用candidate_gain_block，避免分配整段输出。
assert(isnumeric(video) && isreal(video) && ~isempty(video) && ndims(video)<=3, ...
    'iNTA:candidate:InvalidVideo','video必须为非空实数灰度视频。');
if ~isa(video,'uint16') && ~isa(video,'uint8')
    assert(all(isfinite(video(:))) && all(video(:)>=0), ...
        'iNTA:candidate:InvalidVideo','相机视频必须有限且非负。');
end
candidate_gain_block(model,1); % 统一模型校验；不生成整段增益。
[ny,nx,F]=size(video);
assert(numel(model.global_gain)==F,'iNTA:candidate:InvalidVideo','模型与视频帧数不一致。');
if model.is_spatial
    assert(size(model.spatial_coefficients,1)==ny && size(model.spatial_coefficients,2)==nx, ...
        'iNTA:candidate:InvalidVideo','空间模型与视频尺寸不一致。');
end
corrected=zeros(size(video),'single');
for first=1:100:F
    idx=first:min(first+99,F); gain=candidate_gain_block(model,idx);
    block=single(video(:,:,idx));
    if any(gain(:)~=1), block=block./gain; end
    assert(all(isfinite(block(:))),'iNTA:candidate:InvalidVideo','单精度校正结果超出数值范围。');
    corrected(:,:,idx)=block;
end
end
