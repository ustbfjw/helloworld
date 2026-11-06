function m = magnetic_moment(d, B, p)
%MAGNETIC_MOMENT  超顺磁团簇粒子在场强 B 下的平均磁矩 m(B) [A·m^2]
%   模型：粒子由直径 d_crys 的超顺磁晶粒组成（团簇型 Fe3O4 粒子的常见情形），
%   每个晶粒的磁矩 m_c = Ms·V_crys，热平衡下沿场方向的平均值由 Langevin 函数给出：
%       <m_c> = m_c · L(ξ),   ξ = m_c·B/(kB·T),   L(ξ) = coth(ξ) - 1/ξ
%   粒子总磁矩（磁性相体积分数 f_mag）：
%       m(B) = f_mag · V_core · Ms · L(ξ)
%   输入：d 可为向量 [m]；B 为标量 [T]
V_core = pi/6 * d(:).^3;                 % 粒子核体积
V_crys = pi/6 * p.part.d_crys^3;         % 单个晶粒体积
xi = p.part.Ms * V_crys * B / (p.phys.kB * p.phys.T);
if abs(xi) < 1e-3
    L = xi/3;                            % 小 ξ 级数展开，避免 coth 数值误差
else
    L = coth(xi) - 1/xi;
end
m = p.part.f_mag * V_core * p.part.Ms * L;
end
