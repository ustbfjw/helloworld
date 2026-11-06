function [tau, msd, msd_dim, npairs] = compute_msd(pos, frames, n_lag, dt, bounds, margin)
%COMPUTE_MSD  含丢帧轨迹的MSD；边界余量用于减轻有限视场截断偏差
%   pos：n×dim 位置 [µm]；frames：n×1 帧号（可有间隙）；n_lag：最大滞后帧数；dt：帧间隔 [s]
%   bounds：检测区域半宽 [µm]（标量或 1×dim；省略则不做边界筛选）
%   margin：n_lag×1 各滞后的边界余量 [µm]
%
%   有限视场偏差：粒子大步跨出检测区域时轨迹终止，留在轨迹内的位移偏向小值，
%   使 MSD 斜率（扩散系数）被低估（对快速扩散的小粒子可达 5%~10%）。
%   消除方法：对滞后 k 只统计"起点"距边界大于 margin(k) ≈ 2.5·σ_k 的位移对，
%   这样终点跨出区域的概率可忽略，位移统计不再受截断影响。
%   输出 msd(k) = <|Δr(k·dt)|²>（各维之和），msd_dim(k,:) 各维分量，npairs(k) 参与平均的位移对数
dim = size(pos, 2);
if nargin < 5, bounds = []; end
if nargin < 6, margin = zeros(n_lag, 1); end
fr  = frames - frames(1) + 1; Lf = fr(end);
P   = nan(Lf, dim); P(fr, :) = pos;             % 用 NaN 填充丢帧
msd = nan(n_lag, 1); msd_dim = nan(n_lag, dim); npairs = zeros(n_lag, 1);
for k = 1:n_lag
    if k >= Lf, break; end
    P0 = P(1:end-k, :);
    dP = P(1+k:end, :) - P0;
    ok = ~any(isnan(dP), 2);
    if ~isempty(bounds)
        ok = ok & all(bsxfun(@lt, abs(P0), bounds(:).' - margin(k)), 2);
    end
    npairs(k) = sum(ok);
    if npairs(k) > 0
        msd_dim(k, :) = mean(dP(ok, :).^2, 1);
        msd(k) = sum(msd_dim(k, :));
    end
end
tau = (1:n_lag).'*dt;
end
