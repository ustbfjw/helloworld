function plot_spatial_illumination_diagnostics(sample,p,tag)
%PLOT_SPATIAL_ILLUMINATION_DIAGNOSTICS 紧凑空间模型及固定留出门控诊断。
%   读取sample.ex.illumination_model及illumination.score，不重新拟合。
%   每频率sqrt(cos系数^2+sin系数^2)是谐波幅度，不是时间RMS。
%   三个位置在看结果前固定；实际增益由生产candidate_gain_block生成。
%   不读取生成器真值，不把空间平均增益或拟合图当作像素误差/性能证据。
if nargin<3, tag='sample'; end
if ~isfield(p,'out') || ~isfield(p.out,'dir')
    fprintf('空间光照诊断未保存：缺少p.out.dir。\n'); return;
end
if ~exist(p.out.dir,'dir'), mkdir(p.out.dir); end
info=[]; model=[]; video_id=nan;
if isstruct(sample) && isscalar(sample)
    if isfield(sample,'ex') && isstruct(sample.ex) && isscalar(sample.ex)
        if isfield(sample.ex,'illumination'), info=sample.ex.illumination; end
        if isfield(sample.ex,'illumination_model'), model=sample.ex.illumination_model; end
    end
    if isfield(sample,'video') && isscalar(sample.video), video_id=sample.video; end
end
valid=isstruct(model) && isscalar(model) && isfield(model,'is_spatial') && ...
    isequal(model.is_spatial,true);
if valid
    try, candidate_gain_block(model,1); catch, valid=false; end
end
if ~valid
    fig=figure('Name',['Spatial illumination diagnostics ' tag],'Color','w','Position',[100 100 1430 470]);
    axis off;
    text(0.5,0.5,{'空间模型诊断不可用：未保存有效紧凑系数/时间基', ...
        '不能用空间平均增益代替像素增益，也不能从生成器真值补造模型'}, ...
        'Units','normalized','HorizontalAlignment','center','Interpreter','none');
    sgtitle([tag '：空间模型字段缺失'],'Interpreter','none');
    save_figure(fig,p,tag); return;
end
ny=size(model.spatial_coefficients,1); nx=size(model.spatial_coefficients,2);
K=size(model.spatial_coefficients,3)/2; F=numel(model.global_gain);
frequencies=nan(1,K);
if isstruct(info) && isscalar(info)
    if isfield(info,'accepted_frequencies') && numel(info.accepted_frequencies)==K
        frequencies=double(info.accepted_frequencies(:).');
    elseif isfield(info,'local_frequencies') && numel(info.local_frequencies)==K
        frequencies=double(info.local_frequencies(:).');
    end
end
% 固定三个相对位置：四分之一、中心、四分之三，不按结果选择亮点。
points=max(1,round([0.25 0.25;0.5 0.5;0.75 0.75].*[ny,nx]));
points(:,1)=min(points(:,1),ny); points(:,2)=min(points(:,2),nx);
point_ids=sub2ind([ny,nx],points(:,1),points(:,2));
colors=[0.06 0.36 0.75;0.86 0.28 0.10;0.12 0.55 0.28];
ncols=max(3,K+1);
fig=figure('Name',['Spatial illumination diagnostics ' tag],'Color','w','Position',[80 80 1650 830]);
layout=tiledlayout(fig,2,ncols,'TileSpacing','compact','Padding','compact');
for k=1:K
    ax=nexttile(layout,k);
    amplitude=hypot(double(model.spatial_coefficients(:,:,2*k-1)), ...
        double(model.spatial_coefficients(:,:,2*k)));
    imagesc(ax,amplitude*100); axis(ax,'image'); colorbar(ax); hold(ax,'on');
    for j=1:3
        plot(ax,points(j,2),points(j,1),'o','Color',colors(j,:), ...
            'MarkerFaceColor','w','LineWidth',1.2,'MarkerSize',7);
        text(ax,points(j,2)+1,points(j,1),sprintf('%d',j),'Color',colors(j,:), ...
            'FontWeight','bold','Interpreter','none');
    end
    if isfinite(frequencies(k)), label=sprintf('频率 %.3f Hz',frequencies(k));
    else, label=sprintf('频率项 %d（Hz未保存）',k); end
    title(ax,{label,'谐波系数幅度（%，非时间RMS）'},'Interpreter','none');
    xlabel(ax,'x / 列像素'); ylabel(ax,'y / 行像素');
end

ax=nexttile(layout,K+1,[1,ncols-K]); hold(ax,'on');
valid_fps=isfield(p,'sim') && isfield(p.sim,'fps') && isnumeric(p.sim.fps) && ...
    isscalar(p.sim.fps) && isfinite(p.sim.fps) && p.sim.fps>0;
if valid_fps
    t_ms=(0:F-1)'/p.sim.fps*1e3; idx=find(t_ms<=100);
    point_gain=zeros(numel(idx),3,'single');
    for first=1:100:numel(idx)
        rows=first:min(first+99,numel(idx));
        actual=reshape(candidate_gain_block(model,idx(rows)),ny*nx,[]);
        point_gain(rows,:)=actual(point_ids,:).';
    end
    for j=1:3
        plot(ax,t_ms(idx),double(point_gain(:,j)),'Color',colors(j,:), ...
            'LineWidth',1.1,'DisplayName',sprintf('点%d：行%d，列%d',j,points(j,1),points(j,2)));
    end
    xlabel(ax,'相对首帧时间 (ms)'); ylabel(ax,'实际 single 增益');
    title(ax,'前100 ms：三个预选像素的实际增益','Interpreter','none');
    legend(ax,'Location','southoutside','Interpreter','none'); grid(ax,'on'); box(ax,'on');
else
    axis(ax,'off'); text(ax,0.5,0.5,'缺少有效帧率，不能绘制前100 ms增益', ...
        'Units','normalized','HorizontalAlignment','center','Interpreter','none');
end

ax=nexttile(layout,ncols+1,[1,ncols-1]);
score=[];
if isstruct(info) && isscalar(info) && isfield(info,'score'), score=info.score; end
score_fields={'local_residual_amplitude','global_residual_amplitude','uncorrected_amplitude'};
valid_score=isstruct(score) && isscalar(score) && all(isfield(score,score_fields));
if valid_score
    for k=1:numel(score_fields)
        value=score.(score_fields{k});
        valid_score=valid_score && isnumeric(value) && isreal(value) && numel(value)==16;
    end
end
if valid_score
    values=[double(score.uncorrected_amplitude(:)), ...
        double(score.global_residual_amplitude(:)),double(score.local_residual_amplitude(:))];
    valid_score=isequal(size(values),[16,3]) && all(isfinite(values(:))) && all(values(:)>=0);
end
if valid_score
    bars=bar(ax,1:16,100*values,'grouped');
    bars(1).FaceColor=[0.65 0.65 0.65]; bars(2).FaceColor=[0.85 0.50 0.10]; bars(3).FaceColor=[0.08 0.42 0.70];
    xticks(ax,1:16); xlim(ax,[0.3,16.7]); grid(ax,'on'); box(ax,'on');
    xlabel(ax,'固定留出验证区域编号'); ylabel(ax,'多频系数向量范数（%，非时间RMS）');
    title(ax,{'16个留出区域：周期分量幅度诊断', ...
        '全局基准为同训练频率/训练像素拟合；不是PN或原full总噪声'},'Interpreter','none');
    legend(ax,{'未校正','训练拟合的全局基准','空间候选残差'}, ...
        'Location','southoutside','Orientation','horizontal','Interpreter','none');
else
    axis(ax,'off'); text(ax,0.5,0.5,'未保存完整的16区域留出评分，不能补造比较', ...
        'Units','normalized','HorizontalAlignment','center','Interpreter','none');
end

ax=nexttile(layout,2*ncols); axis(ax,'off');
lines={'实际使用：空间增益模型','未叠加global_gain',''};
if isstruct(info) && isscalar(info) && isfield(info,'reason')
    lines=[lines,wrap_lines(['门控结论：' char(info.reason)],22)]; %#ok<AGROW>
else
    lines{end+1}='门控结论字段未保存';
end
if valid_score && isfield(score,'relative_energy_reduction') && ...
        isscalar(score.relative_energy_reduction) && isfinite(score.relative_energy_reduction)
    lines{end+1}=sprintf('验证能量相对全局降低：%.2f%%',100*score.relative_energy_reduction);
end
if valid_score && isfield(score,'active_validation_tiles')
    lines{end+1}=sprintf('活跃验证区域：%d / 16',numel(score.active_validation_tiles));
end
lines=[lines,{''},wrap_lines('固定开发门槛，不是统计显著性检验。',22), ...
    {''},wrap_lines('幅度场与增益曲线来自拟合模型，不是照明真值误差或粒径精度证明。',22)];
text(ax,0,0.97,lines,'Units','normalized','VerticalAlignment','top','Interpreter','none','FontSize',10);
sgtitle(layout,sprintf('%s，视频%d：空间照明模型与留出门控（诊断图，非性能证明）',tag,video_id), ...
    'Interpreter','none');
save_figure(fig,p,tag);
end

function lines=wrap_lines(value,width)
lines=cell(1,ceil(numel(value)/width));
for k=1:numel(lines), lines{k}=value((k-1)*width+1:min(k*width,numel(value))); end
end

function save_figure(fig,p,tag)
set(findall(fig,'-property','FontName'),'FontName','Microsoft YaHei');
print(fig,fullfile(p.out.dir,['fig7_illumination_' tag '.png']),'-dpng','-r150');
end
