function [tracks, res, extra] = run_inta_pipeline(p, tab, verbose)
%RUN_INTA_PIPELINE  完整流程：轨迹仿真 → iSCAT 成像 → iNTA 追踪与反演，循环多段视频
%   tracks：所有视频段的轨迹（含反演结果）；res：群体统计；extra：各段的 traj / 示例帧 / 成像元数据
if nargin < 2 || isempty(tab), tab = ipsf_lookup_table(p); end
if nargin < 3, verbose = true; end
tracks = empty_tracks();
extra.traj = cell(1, p.sim.n_videos); extra.ex = cell(1, p.sim.n_videos); extra.meta = cell(1, p.sim.n_videos);
for v = 1:p.sim.n_videos
    t0 = tic;
    traj = simulate_trajectories(p, p.sim.seed + 100*(v-1));         % [1] 三维轨迹
    t1 = toc(t0);
    [video, meta] = render_iscat_video(traj, p, tab);                   % [2] 相机图像序列
    t2 = toc(t0);
    [tr, ex] = inta_track_video(video, p);                              % [3a] 背景扣除/RVT/定位/链接
    tr = analyze_tracks(tr, traj, p, v);                                % [3b] MSD → D → d，∛C，漂移，真值匹配
    t3 = toc(t0);
    tracks = [tracks, tr]; %#ok<AGROW>
    extra.traj{v} = traj; extra.ex{v} = ex; extra.meta{v} = meta;
    if p.out.save_video || p.out.save_diagnostics
        if ~exist(p.out.dir,'dir'), mkdir(p.out.dir); end
        seed = p.sim.seed + 100*(v-1);
        prefix = sprintf('%s_seed%d_video%02d',p.scenario,seed,v);
        if p.out.save_video
            save(fullfile(p.out.dir,[prefix '_camera.mat']),'video','p','seed','-v7.3');
        end
        if p.out.save_diagnostics
            save(fullfile(p.out.dir,[prefix '_diagnostics.mat']),'traj','ex','meta','tr','p','seed','-v7.3');
        end
    end
    if verbose
        fprintf('视频段 %d/%d：轨迹仿真 %.1f s，成像 %.1f s，iNTA 分析 %.1f s；得到 %d 条轨迹（≥%d 点：%d 条）\n', ...
            v, p.sim.n_videos, t1, t2-t1, t3-t2, numel(tr), p.ana.min_len, sum([tr.n_loc] >= p.ana.min_len));
    end
end
res = inta_population_analysis(tracks, p);
end
