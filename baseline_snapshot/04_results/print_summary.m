function print_summary(tracks, res, p, tag)
%PRINT_SUMMARY 单分散组的主估计、模型对照及统计量；真值仅用于最后的评价。
% 主结果来自组内扩散系数；个体直径平均仅作诊断，不等于组代表直径。
if nargin < 4, tag = p.scenario; end
if ~exist(p.out.dir,'dir'), mkdir(p.out.dir); end
fname = fullfile(p.out.dir, ['results_summary_' tag '.txt']);
fid = fopen(fname, 'w', 'n', 'UTF-8');
if fid < 0, error('不能写入汇总文件：%s',fname); end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fids = [1 fid]; r = res.group;
pr(fids, '\n========== 单分散样品 %s ==========\n', tag);
pr(fids, '主模型：%s；合格轨迹 %d 条，定位 %d 点；视频 %d 段，其中有合格轨迹 %d 段。\n', ...
    r.fit_model, r.n_tracks, r.n_loc, r.n_videos, r.n_videos_with_tracks);
pr(fids, '组内扩散系数 D = %.6f µm²/s\n', r.D_avgmsd);
pr(fids, '组代表水动力直径 dH = %.3f nm\n', r.d_avgmsd);
ci = r.d_ci95_interval;
if numel(ci)==2 && all(isfinite(ci))
    pr(fids, '按独立视频 bootstrap 的 95%% CI：[%.3f, %.3f] nm\n',ci(1),ci(2));
else
    pr(fids, '95%% CI：不可用（独立视频数量不足或有效 bootstrap 估计不足）。\n');
end
pr(fids, '组内场同步速度估计 = %.4f µm/s，沿梯度恒定流速估计 = %.4f µm/s\n', ...
    r.v_mean, r.flow_mean);
pr(fids, '\n模型对照（同一批合格轨迹，禁止事后按真值选择模型）：\n');
pr(fids, '%12s %16s %16s\n','模型','D(µm²/s)','dH(nm)');
models = {'linear','perp','waveform'};
for i = 1:numel(models)
    e = r.model_estimates.(models{i});
    pr(fids, '%12s %16.6f %16.3f\n',models{i},e.D,e.d);
end
pr(fids, '\n轨迹级诊断：先逐条求直径再平均 = %.3f nm，标准差 = %.3f nm。\n', ...
    r.d_mean,r.d_std);
pr(fids, '此数值受倒数变换与轨迹长度影响，不作为本组主结果。\n');
if ~isempty(tracks)
    n_invalid = sum(res.sel & (~isfinite([tracks.d_mes]) | [tracks.d_mes]<=0));
    pr(fids, '合格轨迹中没有有效正个体直径的数量：%d（仍保留其有效位移参与组内统计）。\n',n_invalid);
end
% 参考尺寸仅用于评价，不改变 res.sel 或任何扩散拟合。
active = find(p.part.n_per_size>0);
if numel(active)==1
    nominal = (p.part.d_list(active)+2*p.part.t_shell)*1e9;
    pr(fids, '\n仿真验证：设置的 dH = %.3f nm，主估计相对误差 = %+.3f%%。\n', ...
        nominal,100*(r.d_avgmsd/nominal-1));
    pr(fids, '该参考尺寸只用于误差评价，不参与轨迹筛选或估计。\n');
end
pr(fids, '保存位置：%s\n',fname);
end

function pr(fids,varargin)
for f = fids, fprintf(f,varargin{:}); end
end
