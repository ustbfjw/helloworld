function [eps_s, phi_s, C_peak] = particle_contrast(d, n_p, p)
%PARTICLE_CONTRAST  由 Mie 理论 + 装置标定常数得到各粒子的 iSCAT 散射场参数
%   论文 Methods：C = (I_det - I_bg)/I_bg ≈ 2|E_sca/E_ref|cosθ，E_sca ∝ 后向散射振幅。
%   本函数以参考粒子（p.opt.d_ref, p.opt.n_ref, 最大正对比度 p.opt.C_ref）标定装置常数 κ，
%   在参考粒子H=1且同相条件下，C_ref=2ε_ref+ε_ref²，保留散射强度项。
%   C_peak返回2ε+ε²：|H|=1且同相时的理想上界，不等于各材料实际中心峰值。
%   实际对比度由材料相位、离焦及曝光积分共同决定，不能一律用C∝d³。
%   散射相位 φ = arg S1(π) - arg S1_ref(π) + φ0（论文 Methods 中"材料相关散射相位"）。
%   d 可为向量；n_p 为该材料的（复）折射率。
S_ref = mie_backscatter(p.opt.d_ref, p.opt.n_ref, p.opt.n_m, p.opt.lambda);
eps_ref = sqrt(1+p.opt.C_ref)-1;
kappa = eps_ref / abs(S_ref);
d = d(:);
S = zeros(size(d));
for k = 1:numel(d)
    S(k) = mie_backscatter(d(k), n_p, p.opt.n_m, p.opt.lambda);
end
eps_s = kappa*abs(S);
C_peak = 2*eps_s+eps_s.^2;
phi_s  = angle(S) - angle(S_ref) + p.opt.phi0;
end
