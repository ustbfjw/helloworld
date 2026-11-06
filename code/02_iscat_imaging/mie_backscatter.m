function [S1pi, Qsca, Qext] = mie_backscatter(d, n_p, n_m, lambda)
%MIE_BACKSCATTER  均匀球的 Mie 散射：后向散射振幅 S1(180°) 及效率因子
%   按 Bohren & Huffman《Absorption and Scattering of Light by Small Particles》
%   （论文参考文献 24）的 bhmie 递推实现；支持复折射率（吸收性粒子，如 Fe3O4、Au）。
%   iSCAT 探测的是后向散射场，其幅值 ∝ |S1(180°)|；在 Rayleigh 极限
%   S1(π) ≈ -i·x³·(m²-1)/(m²+2) ∝ 极化率 α（论文式 (2)），即 C ∝ d³。
%   输入：d 直径 [m]，n_p 粒子折射率（可复数），n_m 介质折射率，lambda 真空波长 [m]
%   输出：S1pi 复振幅（无量纲），Qsca / Qext 散射 / 消光效率
x  = pi*d*n_m/lambda;                 % 尺寸参数
m  = n_p/n_m;                         % 相对折射率
mx = m*x;
nmax = ceil(x + 4*x^(1/3) + 2);       % 级数截断项数（Wiscombe 判据）
nst  = ceil(max(nmax, abs(mx)) + 15);

% 对数导数 D_n(mx) = ψ'_n(mx)/ψ_n(mx)，向下递推（数值稳定）
Dn = zeros(nst, 1);
for n = nst:-1:2
    Dn(n-1) = n/mx - 1/(Dn(n) + n/mx);
end

% Riccati–Bessel 函数 ψ_n(x) 与 ξ_n(x) = ψ_n - iχ_n，向上递推
psi_prev = cos(x);  psi_cur = sin(x);     % ψ_{-1}, ψ_0
chi_prev = -sin(x); chi_cur = cos(x);     % χ_{-1}, χ_0
xi_cur   = psi_cur - 1i*chi_cur;          % ξ_0

S1pi = 0; Qsca = 0; Qext = 0;
for n = 1:nmax
    psi_new = (2*n-1)/x*psi_cur - psi_prev;   % ψ_n
    chi_new = (2*n-1)/x*chi_cur - chi_prev;   % χ_n
    xi_new  = psi_new - 1i*chi_new;           % ξ_n
    an = ((Dn(n)/m + n/x)*psi_new - psi_cur) / ((Dn(n)/m + n/x)*xi_new - xi_cur);
    bn = ((m*Dn(n) + n/x)*psi_new - psi_cur) / ((m*Dn(n) + n/x)*xi_new - xi_cur);
    % θ = 180°：π_n(-1) = (-1)^(n-1)·n(n+1)/2，τ_n(-1) = -π_n(-1)
    S1pi = S1pi + 0.5*(2*n+1)*(-1)^(n+1)*(an - bn);
    Qsca = Qsca + (2*n+1)*(abs(an)^2 + abs(bn)^2);
    Qext = Qext + (2*n+1)*real(an + bn);
    psi_prev = psi_cur; psi_cur = psi_new;
    chi_prev = chi_cur; chi_cur = chi_new;
    xi_cur   = xi_new;
end
Qsca = 2/x^2*Qsca;
Qext = 2/x^2*Qext;
end
