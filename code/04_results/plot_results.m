function plot_results(tracks, res, p, tag)
%PLOT_RESULTS 单分散组的 MSD 与轨迹诊断；不使用真值标签选择轨迹。
% MSD 图将垂直分量乘2，与二维曲线放在同一尺度显示；主拟合用真实维数。
if nargin < 4, tag = p.scenario; end
if isempty(tracks) || ~any(res.sel)
    fprintf('[%s] 没有合格轨迹，跳过结果图。\n',tag); return;
end
r = res.group; tau = r.tau(:); idx = find(res.sel);
blue = [0.10 0.38 0.66]; orange = [0.85 0.39 0.10]; green = [0.15 0.58 0.36];
outdir = p.out.dir;
if ~exist(outdir,'dir'), mkdir(outdir); end

fig = figure('Name',['MSD ' tag],'Color','w','Position',[70 70 870 530]);
hold on;
plot(tau*1e3,r.msd_avg,'o-','Color',orange,'LineWidth',1.5,'MarkerSize',5);
plot(tau*1e3,2*r.msd_perp_avg,'s-','Color',blue,'LineWidth',1.8,'MarkerSize',5);
plot(tau*1e3,r.msd_corrected_avg,'^-','Color',green,'LineWidth',1.5,'MarkerSize',5);
labels = {'二维 MSD（含驱动）','2 × 垂直方向 MSD','二维 MSD（波形去漂移）'};
model = r.fit_model;
if strcmpi(model,'perp')
    yfit = 2*(2*r.D_avgmsd*tau+r.b_avgmsd);
elseif strcmpi(model,'drift')
    yfit = 4*r.D_avgmsd*tau+r.model_estimates.drift.v^2*tau.^2+r.b_avgmsd;
else
    yfit = 4*r.D_avgmsd*tau+r.b_avgmsd;
end
nf = min(p.ana.n_fit,numel(tau));
plot(tau(1:nf)*1e3,yfit(1:nf),'k--','LineWidth',2);
labels{end+1} = sprintf('主拟合 %s：D=%.3f µm²/s',model,r.D_avgmsd);
legend(labels,'Location','northwest','Interpreter','none');
xlabel('滞后时间 τ (ms)'); ylabel('MSD 或 2 × MSD⊥ (µm²)');
title(sprintf('%s：组代表 dH = %.2f nm；有效轨迹 %d 条',tag,r.d_avgmsd,r.n_tracks), ...
    'Interpreter','none');
grid on; box on;
print(fig,fullfile(outdir,['fig3_msd_' tag '.png']),'-dpng','-r150');

% 个体分布不截掉尾部，且明示其不是主要估计量。
individual = [tracks(idx).d_mes]; nloc = [tracks(idx).n_loc];
valid = isfinite(individual) & individual>0;
fig = figure('Name',['Track diagnostics ' tag],'Color','w','Position',[90 90 1050 440]);
subplot(1,2,1); hold on;
if any(valid)
    d = individual(valid); weights = nloc(valid);
    lo = max(0,min(d)*0.9); hi = max(d)*1.1;
    if hi<=lo, hi=lo+1; end
    [centers,counts] = weighted_hist(d,weights,linspace(lo,hi,32));
    counts = counts/sum(counts);
    bar(centers,counts,1,'FaceColor',[0.73 0.83 0.92],'EdgeColor',blue);
    yl=ylim;
    plot([r.d_avgmsd r.d_avgmsd],yl,'k--','LineWidth',1.5);
    legend('个体直径分布','组代表直径','Location','best');
else
    text(0.5,0.5,'无有效正个体直径','HorizontalAlignment','center');
end
xlabel('个体水动力直径 (nm)'); ylabel('相对频数（按定位点数加权）');
title('诊断分布：不作为本组主结果'); grid on; box on;
subplot(1,2,2); hold on;
if any(valid)
    semilogy(nloc(valid),individual(valid),'o','Color',blue,'MarkerFaceColor',[0.73 0.83 0.92]);
    set(gca,'YScale','log');
    xl=xlim;
    plot(xl,[r.d_avgmsd r.d_avgmsd],'k--','LineWidth',1.5);
end
xlabel('每条轨迹的定位点数'); ylabel('个体直径 (nm，对数轴)');
title(sprintf('共 %d 条合格轨迹；%d 条无有效个体直径',numel(idx),sum(~valid)));
grid on; box on;
print(fig,fullfile(outdir,['fig4_track_diagnostics_' tag '.png']),'-dpng','-r150');
end
