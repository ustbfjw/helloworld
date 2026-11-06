function plot_validation_diagnostics(sample, tracks, p, tag, selected)
%PLOT_VALIDATION_DIAGNOSTICS 真位置、检测位置和定位残差的独立验证图。
% sample含traj/ex/meta/video，必须来自同一视频；selected通常传res.sel。
% 真实粒子ID仅供本图诊断，不修改轨迹、入选规则或任何粒径估计。
% 曝光平均真位置只是几何诊断基准，不是复杂iPSF中心的严格物理定义。
if nargin<4, tag='sample'; end
if nargin<5, selected=true(1,numel(tracks)); end
selected=logical(selected(:).');
assert(numel(selected)==numel(tracks),'selected长度必须与tracks一致。');
if ~exist(p.out.dir,'dir'), mkdir(p.out.dir); end
N=p.cam.Npx; pxu=p.cam.px*1e6;
coords=((1:N)-(N+1)/2)*pxu; fov=N*pxu/2;
green=[0.05 0.75 0.25]; red=[0.95 0.12 0.12]; blue=[0.12 0.42 0.70];
fig=figure('Name',['Validation diagnostics ' tag],'Color','w','Position',[130 130 1140 490]);
ex=[]; truth=[]; video_id=nan;
if isstruct(sample) && isscalar(sample)
    if isfield(sample,'ex'), ex=sample.ex; end
    if isfield(sample,'video') && isscalar(sample.video), video_id=sample.video; end
    if isfield(sample,'traj') && isstruct(sample.traj)
        if isfield(sample.traj,'pos_camera_wrap')
            truth=sample.traj.pos_camera_wrap;
        elseif isfield(sample.traj,'pos_wrap')
            truth=sample.traj.pos_wrap;
        end
    end
end

%% 对比度、曝光平均真位置与同帧检测：全部换算为同一个µm坐标系。
subplot(1,2,1); hold on;
has_frame=isstruct(ex) && isfield(ex,'frame') && isscalar(ex.frame) && ...
    isfinite(ex.frame) && isfield(ex,'C') && ~isempty(ex.C);
if has_frame
    imagesc(coords,coords,double(ex.C));
    set(gca,'YDir','normal'); axis image; colormap(gca,gray); colorbar;
    cl=max(abs(ex.C(:)));
    if isfinite(cl) && cl>0, caxis([-cl cl]); end
    visible_truth=zeros(0,2);
    f=ex.frame;
    if ~isempty(truth) && f>=1 && f==round(f) && f<=size(truth,1)
        xy=reshape(truth(f,:,1:2),[],2)*1e6;
        inside=all(isfinite(xy),2) & abs(xy(:,1))<=fov & abs(xy(:,2))<=fov;
        visible_truth=xy(inside,:);
    end
    detected=zeros(0,2);
    if isfield(ex,'det') && ~isempty(ex.det)
        detected=(ex.det(:,1:2)-(N+1)/2)*pxu;
    end
    h1=plot(visible_truth(:,1),visible_truth(:,2),'+','Color',green,'MarkerSize',11,'LineWidth',1.8);
    h2=plot(detected(:,1),detected(:,2),'o','Color',red,'MarkerSize',12,'LineWidth',1.3);
    legend([h1 h2],{'曝光平均真位置（仅诊断）','图像检测位置'},'Location','southoutside');
    xlim([-fov fov]); ylim([-fov fov]);
    title(sprintf('%s：视频%d，第%d帧',tag,video_id,f),'Interpreter','none');
else
    text(0.5,0.5,'没有可展示的相机示例帧','Units','normalized','HorizontalAlignment','center');
    xlim([-fov fov]); ylim([-fov fov]); axis square;
    title('相机帧诊断不可用');
end
xlabel('x (µm)'); ylabel('y (µm)'); box on;

%% 同一视频中入选且可匹配的轨迹，计算逐定位点径向残差，不用于再筛选。
subplot(1,2,2); hold on;
zall=zeros(0,1); errorall=zeros(0,1); matched_tracks=0;
if ~isempty(truth) && isfinite(video_id) && ~isempty(tracks)
    if isfield(tracks,'video') && isfield(tracks,'truth_id')
        candidates=find(selected & [tracks.video]==video_id & [tracks.truth_id]>0);
        for k=candidates
            t=tracks(k); j=t.truth_id; fr=t.frames(:);
            if ~isscalar(j) || j~=round(j) || j>size(truth,2), continue; end
            xy=[t.x_um(:),t.y_um(:)];
            if size(xy,1)~=numel(fr), continue; end
            ok=fr>=1 & fr<=size(truth,1) & fr==round(fr) & all(isfinite(xy),2);
            if ~any(ok), continue; end
            actual=reshape(truth(fr(ok),j,:),[],3);
            err=sqrt(sum((xy(ok,:)-actual(:,1:2)*1e6).^2,2))*1e3;
            good=isfinite(err) & isfinite(actual(:,3));
            if ~any(good), continue; end
            zall=[zall;actual(good,3)*1e6]; %#ok<AGROW>
            errorall=[errorall;err(good)]; %#ok<AGROW>
            matched_tracks=matched_tracks+1;
        end
    end
end
if ~isempty(errorall)
    plot(zall,errorall,'.','Color',blue,'MarkerSize',5);
    rms_error=sqrt(mean(errorall.^2));
    xl=[0 p.sim.H_z*1e6];
    plot(xl,[rms_error rms_error],'k--','LineWidth',1.3);
    xlim(xl);
    legend('入选且可匹配的定位点','径向 RMS','Location','northwest');
    title(sprintf('%d条轨迹，%d点；径向 RMS = %.2f nm',matched_tracks,numel(errorall),rms_error));
else
    text(0.5,0.5,'该视频没有入选且可匹配的定位点','Units','normalized','HorizontalAlignment','center');
    title('定位残差诊断不可用');
end
xlabel('曝光平均真高度 z (µm)'); ylabel('定位径向误差 (nm)'); grid on; box on;
print(fig,fullfile(p.out.dir,['fig6_validation_' tag '.png']),'-dpng','-r150');
end
