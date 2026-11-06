function H = ipsf_eval(tab, rho, u)
%IPSF_EVAL  在查找表上插值得到复数 iPSF 幅值 H(ρ, u)
%   rho：任意形状数组 [m]；u：标量离焦量 [m]
%   先沿 u 线性插值得到一条径向剖面 Hu(ρ)，再沿 ρ 做向量化线性插值（比 interp1 快得多）
fu = (u - tab.u(1))/tab.du + 1;
i0 = min(max(floor(fu), 1), numel(tab.u) - 1);
wu = min(max(fu - i0, 0), 1);
Hu = (1 - wu)*tab.H(:, i0) + wu*tab.H(:, i0 + 1);      % Nρ×1 复数剖面
nr = numel(tab.rho);
fr = rho(:)/tab.drho + 1;                               % 以 1 为起点的分数索引（列向量）
j0 = min(max(floor(fr), 1), nr - 1);
wr = min(max(fr - j0, 0), 1);
H  = reshape((1 - wr).*Hu(j0) + wr.*Hu(j0 + 1), size(rho));   % 与 rho 同形状
end
