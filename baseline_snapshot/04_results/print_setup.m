function print_setup(p)
%PRINT_SETUP 显示三组单分散实验的物理参数与明确的模型适用范围。
d = p.part.d_list(:); dH = d + 2*p.part.t_shell;
gamma = 3*pi*p.phys.eta*dH;
D = p.phys.kB*p.phys.T./gamma;
m = magnetic_moment(d, p.mag.B0, p);
v = m*p.mag.G0./gamma;
[eps_s, ~, Cupper] = particle_contrast(d, p.part.n_p, p);
fprintf('样品：三组独立单分散磁纳米颗粒，不进行混合样品分群。\n');
fprintf('介质：纯水，T = %.2f K，η = %.3f mPa·s，n = %.3f\n', ...
    p.phys.T, p.phys.eta*1e3, p.opt.n_m);
fprintf('磁场：固定偏置 B0 = %.1f mT，交变梯度 G0 = %.1f T/m，f = %.1f Hz，%s\n', ...
    p.mag.B0*1e3, p.mag.G0, p.mag.f, p.mag.wave);
fprintf('梯度方向：[%g %g %g]；主估计：%s\n', p.mag.dir, p.ana.fit_model);
fprintf('背景模式：%s；模拟周期照明振幅 %.1f%%，频率 %.3f Hz（仅生成器真值）。\n', ...
    p.ana.background_mode,100*p.cam.flicker_amplitude,p.cam.flicker_frequency);
fprintf('相机：%d×%d px，%.2f nm/px，%g fps；曝光 %.1f µs，曝光积分点数 %d\n', ...
    p.cam.Npx, p.cam.Npx, p.cam.px*1e9, p.sim.fps, p.cam.exposure*1e6, p.sim.n_sub);
fprintf('背景光子数：每像素每帧 %.0f；三组使用相同曝光与光子设置。\n',p.cam.photons);
fprintf('采集协议按已知单分散样品预先设置，不按反演结果或真值误差自动选参。\n');
if isfield(p,'scenario')
    names={p.scenario};
else
    names={'mono30','mono100','mono200'};
end
for k=1:numel(names)
    if isfield(p,'scenario'), ps=p; else, ps=apply_scenario(p,names{k}); end
    fprintf('  %s：每段%d个初始粒子，%d段×%d帧，定位半窗%d px（%d×%d）。\n', ...
        names{k},sum(ps.part.n_per_size),ps.sim.n_videos,ps.sim.n_frames, ...
        ps.ana.loc_win,2*ps.ana.loc_win+1,2*ps.ana.loc_win+1);
end
fprintf('强散射组采用低浓度，减少其他粒子远场干涉；不代表已解决高浓度重叠定位。\n');
fprintf('光学：λ = %.0f nm，标称 NA = %.2f，标量模型有效传播 NA = %.3f，焦高 %.2f µm\n', ...
    p.opt.lambda*1e9, p.opt.NA, min(p.opt.NA/p.opt.n_m,0.98)*p.opt.n_m, p.opt.z_f*1e6);
fprintf('模型范围：体相各向同性扩散、独立超顺磁团簇、简化标量反射 iSCAT。\n');
fprintf('下表仅为模拟真值参考；反演不使用已知直径、磁矩或理论磁速度。\n');
fprintf('%10s %10s %14s %15s %13s %16s\n', ...
    '核径(nm)', 'dH(nm)', 'D(µm²/s)', 'v_mag(µm/s)', '场幅比ε', '同相对比度上界');
for k = 1:numel(d)
    fprintf('%10.1f %10.1f %14.5f %15.5f %13.5f %16.5f\n', ...
        d(k)*1e9, dH(k)*1e9, D(k)*1e12, v(k)*1e6, eps_s(k), Cupper(k));
end
fprintf('同相对比度上界 = 2ε+ε²；实际对比度还取决于材料相位及离焦。\n');
end
