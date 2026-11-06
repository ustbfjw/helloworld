function summary = plot_group_comparison(results,p0)
%PLOT_GROUP_COMPARISON 三个独立样品的主估计、视频 bootstrap CI 和参考误差。
% 用法：summary=plot_group_comparison(results,p0)。仅汇总，绝不重估或合并样品。
% 即使关闭绘图，也写出 UTF-8 文本及 CSV，便于验收与重复运行对照。
names = {'mono30','mono100','mono200'};
summary = repmat(struct('scenario','','reference_dH_nm',nan,'estimated_dH_nm',nan, ...
    'ci_low_nm',nan,'ci_high_nm',nan,'relative_error_percent',nan,'D_um2_s',nan, ...
    'n_tracks',0,'n_localizations',0,'n_videos',0,'n_videos_with_tracks',0,'model',''),1,3);
for k=1:3
    item=results.(names{k}); r=item.res.group; p=item.p;
    active=find(p.part.n_per_size>0);
    assert(numel(active)==1,'三组汇总仅接受单分散样品。');
    nominal=(p.part.d_list(active)+2*p.part.t_shell)*1e9;
    summary(k).scenario=names{k};
    summary(k).reference_dH_nm=nominal;
    summary(k).estimated_dH_nm=r.d_avgmsd;
    summary(k).ci_low_nm=r.d_ci95_interval(1); summary(k).ci_high_nm=r.d_ci95_interval(2);
    summary(k).relative_error_percent=100*(r.d_avgmsd/nominal-1);
    summary(k).D_um2_s=r.D_avgmsd;
    summary(k).n_tracks=r.n_tracks; summary(k).n_localizations=r.n_loc;
    summary(k).n_videos=r.n_videos; summary(k).n_videos_with_tracks=r.n_videos_with_tracks;
    summary(k).model=r.fit_model;
end
if ~exist(p0.out.dir,'dir'), mkdir(p0.out.dir); end
write_summary(summary,p0.out.dir);
if ~p0.out.make_plots, return; end
cols=size_colors(3); ref=[summary.reference_dH_nm]; est=[summary.estimated_dH_nm];
lo=[summary.ci_low_nm]; hi=[summary.ci_high_nm];
fig=figure('Name','Three independent groups','Color','w','Position',[110 110 1080 470]);
subplot(1,2,1); hold on;
upper=max([ref est(isfinite(est)) hi(isfinite(hi))])*1.15;
plot([0 upper],[0 upper],'k--','LineWidth',1.2,'DisplayName','理想一致线');
for k=1:3
    if isfinite(est(k))
        plot(ref(k),est(k),'o','Color',cols(k,:),'MarkerFaceColor',cols(k,:),'MarkerSize',8, ...
            'DisplayName',sprintf('%g nm 独立组',ref(k)));
        % 直接绘制置信区间端点，避免人为截断不对称 bootstrap 区间。
        if all(isfinite([lo(k),hi(k)]))
            plot([ref(k) ref(k)],[lo(k) hi(k)],'-','Color',cols(k,:),'LineWidth',1.7,'HandleVisibility','off');
            plot(ref(k)+[-3 3],[lo(k) lo(k)],'-','Color',cols(k,:),'HandleVisibility','off');
            plot(ref(k)+[-3 3],[hi(k) hi(k)],'-','Color',cols(k,:),'HandleVisibility','off');
        end
    end
end
axis([0 upper 0 upper]); axis square; grid on; box on;
xlabel('设置的水动力直径，仅作参考 (nm)'); ylabel('反演水动力直径 (nm)');
title('主估计与按独立视频 bootstrap 的 95% CI'); legend('Location','northwest');
if any(~isfinite(lo) | ~isfinite(hi))
    text(0.04,0.05,'部分组 CI 不可用，详见文本汇总','Units','normalized','FontSize',9);
end
subplot(1,2,2); hold on;
bias=[summary.relative_error_percent];
for k=1:3, bar(k,bias(k),0.6,'FaceColor',cols(k,:),'EdgeColor','none'); end
plot([0.5 3.5],[5 5],'k--'); plot([0.5 3.5],[-5 -5],'k--');
plot([0.5 3.5],[0 0],'k-');
set(gca,'XTick',1:3,'XTickLabel',arrayfun(@(x)sprintf('%g nm',x),ref,'UniformOutput',false));
xlim([0.5 3.5]); bound=max([7 abs(bias(isfinite(bias)))*1.2]); ylim([-bound bound]);
ylabel('相对参考尺寸的误差 (%)'); title('三组分别评价；虚线为 ±5% 参考目标');
grid on; box on;
print(fig,fullfile(p0.out.dir,'fig5_group_comparison.png'),'-dpng','-r150');
end

function write_summary(s,outdir)
fid=fopen(fullfile(outdir,'summary_three_groups.csv'),'w','n','UTF-8');
if fid<0, error('无法写入三组 CSV 汇总。'); end
clean=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'scenario,reference_dH_nm,estimated_dH_nm,ci_low_nm,ci_high_nm,relative_error_percent,D_um2_s,n_tracks,n_localizations,n_videos,n_videos_with_tracks,model\n');
for k=1:numel(s)
    r=s(k);
    fprintf(fid,'%s,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%d,%d,%d,%d,%s\n',r.scenario, ...
        r.reference_dH_nm,r.estimated_dH_nm,r.ci_low_nm,r.ci_high_nm,r.relative_error_percent, ...
        r.D_um2_s,r.n_tracks,r.n_localizations,r.n_videos,r.n_videos_with_tracks,r.model);
end
fid2=fopen(fullfile(outdir,'summary_three_groups.txt'),'w','n','UTF-8');
if fid2<0, error('无法写入三组文本汇总。'); end
clean2=onCleanup(@() fclose(fid2)); %#ok<NASGU>
fprintf(fid2,'三个独立单分散样品：每组单独估计 D，再换算水动力直径。\n');
fprintf(fid2,'参考尺寸仅用于验证，不输入估计器。CI 按独立视频 bootstrap；NaN 表示统计量不足。\n\n');
fprintf(fid2,'%10s %12s %12s %23s %12s %10s %12s\n','样品','参考dH(nm)','反演dH(nm)','95%CI(nm)','偏差(%)','轨迹数','有效视频数');
for k=1:numel(s)
    r=s(k);
    fprintf(fid2,'%10s %12.3f %12.3f [%9.3f,%9.3f] %12.3f %10d %12d\n',r.scenario, ...
        r.reference_dH_nm,r.estimated_dH_nm,r.ci_low_nm,r.ci_high_nm, ...
        r.relative_error_percent,r.n_tracks,r.n_videos_with_tracks);
end
fprintf(fid2,'\n单个结果接近参考尺寸不构成完整验收；请结合多种子验收报告和图像/轨迹检查。\n');
end
