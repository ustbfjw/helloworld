function p = sim_params(preset)
%SIM_PARAMS 仿真与分析的基础参数；分组粒子数、视频数与定位窗见apply_scenario。
%   p = sim_params()          正式测量：30nm组20段，100/200nm组各60段（apply_scenario设置）
%   p = sim_params('quick')   功能演示，不用于精度验收
%
%   单位约定：所有物理量使用 SI 单位（m, s, kg, T, A/m），
%   仅在绘图与结果汇总时换算为 nm / µm。
%
%   依据论文：Kashkanova et al., Nature Methods 19, 586–593 (2022)

if nargin < 1, preset = 'full'; end

%% ---------- 1. 物理常数与介质（纯水） ----------
p.phys.kB    = 1.380649e-23;  % 玻尔兹曼常数 [J/K]
p.phys.T     = 294.15;        % 温度 [K]（21 °C，与论文 Fig.2b 一致）
p.phys.eta   = 0.98e-3;       % 水的动力黏度 @21 °C [Pa·s]
p.phys.rho_f = 998;           % 水的密度 [kg/m^3]
p.phys.g     = 9.81;          % 重力加速度 [m/s^2]

%% ---------- 2. 磁纳米粒子（Fe3O4 团簇型，超顺磁） ----------
p.part.d_list     = [30 100 200]*1e-9; % 核直径 [m]（题目设定：30/100/200 nm）
p.part.n_per_size = [8 0 0];           % 默认30 nm稀样品，apply_scenario分别设置三组
p.part.t_shell    = 0;                 % 水化/包覆层厚度 [m]（流体力学直径 d_H = d + 2 t_shell）
p.part.rho_p      = 5200;              % Fe3O4 密度 [kg/m^3]（用于沉降）
p.part.Ms         = 4.0e5;             % 饱和磁化强度 [A/m]（纳米 Fe3O4 典型 3~4.8e5）
p.part.d_crys     = 12e-9;             % 超顺磁晶粒直径 [m]（Langevin 函数参数）
p.part.f_mag      = 0.6;               % 磁性相体积分数（团簇/包覆导致 <1）
p.part.n_p        = 2.40 + 0.50i;      % Fe3O4 复折射率 @525 nm（文献近似值，可修改）

%% ---------- 3. 交变梯度磁场 ----------
% 场模型：|B|(r,t) = B0 + G(t)·(r·e_g)，  G(t) = G0·s(t)，s(t) 为 ±1 方波或正弦
% 粒子受力：F_mag = m(|B|)·∇|B| ≈ m(B0)·G(t)·e_g   （视场内 G·r << B0）
p.mag.B0     = 20e-3;    % 均匀偏置磁场 [T]（用于磁化粒子）
p.mag.G0     = 50;       % 梯度幅值 [T/m]
p.mag.f      = 10;       % 交变频率 [Hz]
p.mag.wave   = 'square'; % 'square' | 'sine'
p.mag.dir    = [1 0 0];  % 梯度方向（单位向量；默认沿成像面 x 轴，便于在图像中观察往复漂移）
p.mag.phase0 = 0;        % 初相位 [rad]

%% ---------- 4. 轨迹仿真 ----------
p.sim.fps      = 5000;   % 相机帧率 [fps]（论文：5000 fps）
p.sim.n_sub    = 4;      % 每次曝光的强度积分采样数，增加至8可检查收敛性
p.sim.n_frames = 5000;   % 每段视频帧数（1 s，与论文一致）
p.sim.n_videos = 20;     % 30nm组20段；apply_scenario将强散射组设为60段，用于按视频bootstrap
p.sim.L_xy     = 16e-6;  % 仅为初始xy均匀采样区域边长；之后开域运动，不循环注入
p.sim.boundary_xy = 'open'; % xy无反射/周期边界，所有粒子完整影响相机视场
p.sim.H_z      = 4e-6;   % 盒高[m]，球心反射边界；扩散仍采用体相近似
p.sim.z_init   = [];     % 初始 z 范围 [zmin zmax]（空 = 盒内均匀）；慢速大粒子可设为探测区附近以提高短视频利用率
p.sim.gravity  = true;   % 是否包含重力-浮力沉降（对 ≤200 nm 粒子影响很小）
p.sim.seed     = 1;      % 随机种子

%% ---------- 5. 反射式 iSCAT 光学模型 ----------
p.opt.lambda = 525e-9;   % 波长 [m]（论文：525 nm）
p.opt.NA     = 1.46;     % 物镜数值孔径（论文：63× 油镜，NA 1.46）
p.opt.n_m    = 1.334;    % 水的折射率
p.opt.z_f    = 1.0e-6;   % 焦平面距盖玻片高度 [m]（论文：~1 µm）
p.opt.phi0   = 0;        % 装置附加相位 [rad]
p.opt.u_max  = 3.0e-6;   % 查找表最小离焦覆盖范围，自动扩展以覆盖整个盒厚
p.opt.R_max  = 24e-6;    % 初始查表径向尺度；不是PSF截断，渲染按实际轨迹自动扩展
p.opt.C_min  = 0;        % 默认不丢弃任何已有粒子的尾场；不按信噪比截断光场
p.opt.drho   = 20e-9;    % iPSF 查找表径向步长 [m]
p.opt.du     = 20e-9;    % iPSF 查找表轴向步长 [m]
% 对比度标定（论文 Supplementary 2.1 思路：一个装置标定常数）：
% 以论文中 30 nm 金纳米粒子的最大正对比度（∛C ≈ 0.65 → C ≈ 0.275）为参考，
% 由完整强度公式求ε_ref；其它粒子的场幅值按Mie后向散射振幅比例缩放
p.opt.C_ref  = 0.275;          % 参考粒子最大正对比度
p.opt.d_ref  = 30e-9;          % 参考粒子直径 [m]
p.opt.n_ref  = 0.63 + 2.07i;   % 参考粒子（金）折射率 @525 nm（论文采用值）

%% ---------- 6. 相机模型 ----------
p.cam.Npx     = 128;          % 像素数（128×128，论文设置）
p.cam.px      = 7e-6/128;     % 物面像素尺寸 [m]（7×7 µm² 视场）
p.cam.exposure= 80e-6;       % 曝光时间[s]；0为瞬时模型，须≤1/fps
p.cam.photons = 5e4;          % 背景每像素每帧光子数（决定散粒噪声 ~0.45%；论文对小粒子用更高功率）
p.cam.read_e  = 10;           % 读出噪声 [e- rms]
p.cam.bits    = 12;           % ADC 位数（论文相机为 12 bit；若出现饱和提示可改 16）
p.cam.gain    = 800/5e4;      % 增益 [ADU/e-]，使背景约 800 ADU（为强散射粒子留出 ~5 倍动态范围）
p.cam.illum_w = 12e-6;        % 照明高斯光斑 1/e² 半径 [m]（视场内轻微不均匀）
p.cam.bg_amp  = 0.03;         % 静态背景结构（盖玻片粗糙度等）相对幅度
p.cam.bg_corr = 0.4e-6;       % 静态背景相关长度 [m]
% 可关闭的周期照明压力测试；默认静态照明保持上一版物理条件。
% 此处是真值生成参数，背景校正算法不读取它们。
p.cam.flicker_amplitude = 0; % 相对强度振幅；例如0.08表示±8%
p.cam.flicker_frequency = 423.7; % 非整数Hz用于检验非DFT栅格周期
p.cam.flicker_phase = 0.37;  % 初相位[rad]

%% ---------- 7. iNTA 分析参数 ----------
p.ana.bg_block  = 1000;  % 中值背景的时间块长度 [帧]
p.ana.background_mode = 'spectral_common'; % 新融合；'median'恢复原预处理
% 只识别共同周期照明，不对逐像素粒子视频作单频投影。
p.ana.illum = struct('min_amplitude',0.002,'noise_factor',8, ...
    'min_coherence',0.95,'max_amplitude_cv',0.3,'max_frequencies',3, ...
    'protect_harmonics',5,'protect_bins',2);
p.ana.bg_masked = true;  % 检测到强散射粒子时，用带掩膜的中值重算背景（消除慢速强粒子的鬼影）
p.ana.strong_factor = 50;% "强检测"判据：RVT > strong_factor × 阈值
p.ana.mask_r    = 30;    % 掩膜半径 [px]（覆盖强粒子的离焦环）
p.ana.rvt_radii = [2 3 4 5 6 8 10 13 16]; % RVT 环半径 [px]（较大的 r_max 可覆盖离焦粒子的外环）
p.ana.det_win   = 5;     % 局部极大搜索半宽 [px]
p.ana.det_k     = 16;    % 检测阈值 = median + det_k × 1.4826·MAD（纯噪声帧的 RVT 最大值约为 10~14 MAD）
p.ana.loc_win   = 7;     % 细定位窗口半宽 [px]（在对比度图上做径向对称中心）
p.ana.loc_iter  = 3;     % 细定位迭代次数
p.ana.min_sep   = 10;    % 同帧两检测的最小间距 [px]（细定位后合并同一粒子的重复检测）
p.ana.ring_excl = 20;    % 环伪影排除半径 [px]：粗候选若落在更强候选的此范围内且 RVT < ring_ratio×后者则剔除
p.ana.ring_ratio= 1.0;   % 取 1.0 即"每个环系统只保留最强候选"（稀样品推荐；两粒子相距 <1 µm 时会合并）
p.ana.link_r    = 6;     % 轨迹链接最大位移 [px]
p.ana.link_gap  = 3;     % 允许的最大丢帧数
p.ana.min_len   = 40;    % 参与统计的最短轨迹长度（论文：25 或 100）
p.ana.n_lag     = 10;    % 计算 MSD 的最大滞后帧数
p.ana.n_fit     = 4;     % 拟合 MSD 使用的滞后数（短滞后可减小有限视场带来的截断偏差）
p.ana.edge_nsig = 2.5;   % MSD 边界余量系数（起点距检测区域边界 > edge_nsig·σ_k 才计入）
p.ana.fit_model = 'perp'; % 主结果用垂直磁梯度方向；linear/drift作诊断
p.ana.weight_by_len = true;  % 个体直方图权重；主估计用有效位移对数
p.ana.min_pairs = 5;        % 可参与组统计的最少位移对
p.ana.bootstrap_n = 500;    % 独立视频重抽样次数
p.ana.bootstrap_seed = 31415; % 不改变运动/相机随机数序列

%% ---------- 8. 输出 ----------
p.out.dir        = fullfile(fileparts(mfilename('fullpath')), '..', '04_results', 'output');
p.out.make_plots = true;
p.out.save_video = false;    % 是否把相机视频存为 .mat（每段约 160 MB）
p.out.save_diagnostics = true; % 保存真轨迹、示例帧、成像元数据供追溯

%% ---------- 快速测试预设 ----------
p.preset = preset;
if strcmpi(preset, 'quick')
    p.sim.n_frames = 1500;
    p.sim.n_videos = 3;
    p.ana.min_len  = 25;
elseif ~strcmpi(preset,'full')
    error('预设只能是full或quick。');
end
end
