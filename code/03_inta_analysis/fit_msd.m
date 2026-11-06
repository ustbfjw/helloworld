function [D, b, v] = fit_msd(tau, msd_xy, msd_perp, p)
%FIT_MSD  由 MSD 曲线拟合扩散系数 D [µm²/s]（论文：对单条轨迹的 MSD 做线性拟合）
%   'linear' : MSD_xy = 4Dτ + b             论文方法；b 包含定位误差等贡献（截距自由）
%   'drift'  : MSD_xy = 4Dτ + v²τ² + b      显式考虑恒定漂移（交变场半周期内近似恒速）
%   'perp'   : MSD_⊥  = 2Dτ + b             仅用垂直于磁场梯度方向的分量，对磁泳漂移免疫
%   注：对于方波交变梯度场，在拟合所用的短滞后（≈1 ms）内漂移项 v²τ² 远小于 4Dτ，
%       因此论文的线性拟合仍然适用；'drift' / 'perp' 可用于对比验证。
nf = p.ana.n_fit;
D = nan; b = nan; v = nan;
switch lower(p.ana.fit_model)
    case 'linear', y = msd_xy;   dim = 2; deg = 1;
    case 'drift',  y = msd_xy;   dim = 2; deg = 2;
    case 'perp',   y = msd_perp; dim = 1; deg = 1;
    otherwise, error('未知拟合模型: %s', p.ana.fit_model);
end
tau=tau(:); y=y(:);
ok = isfinite(y) & isfinite(tau); ok(nf+1:end) = false;
if sum(ok) < deg + 1, return; end
t_ms = tau(ok)*1e3;                        % 用 ms 作时间单位，避免病态拟合
c = polyfit(t_ms, y(ok), deg);
if deg == 1
    D = c(1)/(2*dim)*1e3; b = c(2);
else
    if c(1)<0
        % v²不能为负；退回v=0的边界解，不能只把显示速度截零。
        c=polyfit(t_ms,y(ok),1);
        D=c(1)/(2*dim)*1e3; b=c(2); v=0;
    else
        D = c(2)/(2*dim)*1e3; b = c(3); v = sqrt(c(1))*1e3;
    end
end
end
