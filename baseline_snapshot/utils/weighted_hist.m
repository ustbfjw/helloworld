function [centers, counts] = weighted_hist(x, w, edges)
%WEIGHTED_HIST  加权直方图（等宽 bin）。论文中直方图按轨迹长度加权（Supplementary 3.2）
%   x: 数据向量；w: 权重向量；edges: 等间距 bin 边界
x = x(:); w = w(:);
nb  = numel(edges) - 1;
dx  = edges(2) - edges(1);
idx = floor((x - edges(1))/dx) + 1;
ok  = idx >= 1 & idx <= nb & isfinite(x);
counts  = accumarray(idx(ok), w(ok), [nb 1]);
centers = (edges(1:end-1) + edges(2:end))/2;
end
