function p = apply_scenario(p, name)
%APPLY_SCENARIO 三组独立单分散样品，按散射强度设置稀释采集协议。
% 名义直径只用于生成和验收比较，不作为扩散反演的约束。
% 强散射粒子的远处干涉环也会影响定位，因此每段只模拟一个粒子，
% 用更多独立视频保持总粒子时长；不是把轨迹按真实误差筛掉。
% 初始xy仍在16um区域均匀分布，不向视场中心注入粒子。
switch lower(name)
    case 'mono30'
        p.part.n_per_size = [8 0 0]; k = 1;
        p.ana.loc_win = 7;
    case 'mono100'
        p.part.n_per_size = [0 1 0]; k = 2;
        p.ana.loc_win = 15;
    case 'mono200'
        p.part.n_per_size = [0 0 1]; k = 3;
        p.ana.loc_win = 15;
    otherwise, error('只支持mono30、mono100、mono200三组独立样品。');
end
% 较大窗口由离焦定位测试确定，覆盖更多环结构；不改变MSD拟合滞后。
% quick仅为3段功能演示；正式强散射组60段与原3粒子×20段总量一致。
if k>1 && ~strcmpi(p.preset,'quick'), p.sim.n_videos=60; end
p.scenario = name;
p.sample.size_index = k;
p.sample.nominal_nm = (p.part.d_list(k)+2*p.part.t_shell)*1e9;
end
