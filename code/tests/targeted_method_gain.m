function [gain,info]=targeted_method_gain(video,p,method)
%TARGETED_METHOD_GAIN 六种方法均只读取同一原始相机视频及固定分析参数。
% PN使用每帧全视野总功率，覆盖全部帧，不作频率筛选或先用本方法去噪。
F=size(video,3);
switch method
    case 'median'
        gain=ones(F,1); info=struct('applied',false,'accepted_frequencies',[]);
    case 'pn'
        gain=double(squeeze(sum(sum(single(video),1),2))); gain=gain(:);
        gain=gain/median(gain);
        assert(all(isfinite(gain) & gain>0),'PN估计增益必须为有限正值。');
        info=struct('applied',any(gain~=1),'accepted_frequencies',[]);
    otherwise
        assert(any(strcmp(method,{'full','no_spatial_gate','no_magnetic_guard','no_refinement'})), ...
            '未知对照方法。');
        p.ana.illum.enable_spatial_gate=~strcmp(method,'no_spatial_gate');
        p.ana.illum.enable_magnetic_guard=~strcmp(method,'no_magnetic_guard');
        p.ana.illum.refine_frequency=~strcmp(method,'no_refinement');
        [gain,info]=estimate_common_illumination(video,p);
end
info.method=method;
end
