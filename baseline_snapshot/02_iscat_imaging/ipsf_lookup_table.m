function tab = ipsf_lookup_table(p, rho_max)
%IPSF_LOOKUP_TABLE  反射式 iSCAT 干涉点扩散函数 (iPSF) 查找表 H(ρ, u)
%
%   物理图像（论文 Fig.1a/1b，参考文献 18, 25 的标量简化模型）：
%   - 参考光：照明平面波在盖玻片/水界面反射，成像后为均匀背景 E_ref；
%   - 散射光：高度 z 处粒子发出的球面波，经高 NA 物镜成像为 iPSF。
%   以参考场为单位，散射场在像面（物面尺度）上的归一化分布为 Debye 积分：
%       H(ρ,u) = (1/N)·∫_0^θmax sqrt(cosθ)·sinθ·J0(k_m ρ sinθ)·exp(i·k_m·u·(1+cosθ)) dθ
%   其中 ρ 为到粒子像中心的径向距离，u = z - z_f 为离焦量，k_m = 2π n_m/λ。
%   相位 exp(i k_m u) 来自照明光到达粒子的额外光程（相对参考光），exp(i k_m u cosθ) 为离焦相位；
%   二者共同给出 iPSF 随轴向位置的周期性明暗交替（论文"对比度沿轴向调制"）与离焦环。
%   相机强度：I ∝ |1 + ε·e^{iφ}·H(ρ,u)|²，ε = |E_sca/E_ref|（见 render_iscat_video）。
%
%   R_max仅为初始查表尺度，不是PSF的物理截止半径。可传入所需rho_max。
%   输出 tab.rho [1×Nρ]、tab.u [1×Nu]、tab.H [Nρ×Nu]（复数，H(0,0)=1）
if nargin<2, rho_max=p.opt.R_max*sqrt(2); end
k     = 2*pi*p.opt.n_m/p.opt.lambda;
smax  = min(p.opt.NA/p.opt.n_m, 0.98);         % 水中可被收集的最大 sinθ
thmax = asin(smax);
Nth   = 800;
th    = linspace(0, thmax, Nth);
dth   = th(2) - th(1);
w     = sqrt(cos(th)).*sin(th)*dth;            % 积分权重（含 aplanatic 切趾因子）
rho   = 0 : p.opt.drho : (ceil(rho_max/p.opt.drho)+2)*p.opt.drho;
% 覆盖整个仿真盒，不能在仍有可见信号的固定离焦距离硬删除粒子。
umax = max([p.opt.u_max,abs(p.opt.z_f),abs(p.sim.H_z-p.opt.z_f)]);
umax = ceil(umax/p.opt.du)*p.opt.du;
u = -umax : p.opt.du : umax;
J     = besselj(0, k*(sin(th).')*rho);         % Nth×Nρ（只需计算一次）
H     = zeros(numel(rho), numel(u));
for iu = 1:numel(u)
    ph = exp(1i*k*u(iu)*(1 + cos(th)));        % 1×Nth
    H(:, iu) = ((w.*ph)*J).';                  % 注意用 .'（非共轭转置）
end
H = H / sum(w);                                % 归一化：H(ρ=0,u=0) = 1
tab.rho = rho; tab.u = u; tab.H = H;
tab.drho = p.opt.drho; tab.du = p.opt.du; tab.k = k;
tab.effective_NA = smax*p.opt.n_m;   % 传播光标量近似；非完整油镜/界面矢量模型
tab.Hmax_u = max(abs(H), [], 1);              % 各离焦量下的最大 |H|（用于渲染截止判断）
end
