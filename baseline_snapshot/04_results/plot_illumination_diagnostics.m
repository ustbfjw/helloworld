function plot_illumination_diagnostics(sample,p,tag)
%PLOT_ILLUMINATION_DIAGNOSTICS 共同周期照明的独立诊断图，不重新估计背景。
%   sample.ex.illumination保存估计器诊断；sample.meta.illumination_gain仅作
%   生成器真值对照，绝不用于追踪、入选或粒径计算。前后RMS只指空间共同
%   亮度序列的相对波动，不能解释为全部背景或散粒噪声被消除。
if nargin<3, tag='sample'; end
if ~isfield(p,'out') || ~isfield(p.out,'dir')
    fprintf('光照诊断未保存：缺少p.out.dir。\n'); return;
end
if ~exist(p.out.dir,'dir'), mkdir(p.out.dir); end
info=[]; truth=[]; video_id=nan;
if isstruct(sample) && isscalar(sample)
    if isfield(sample,'ex') && isstruct(sample.ex) && isfield(sample.ex,'illumination')
        info=sample.ex.illumination;
    end
    if isfield(sample,'meta') && isstruct(sample.meta) && isfield(sample.meta,'illumination_gain')
        truth=sample.meta.illumination_gain(:);
    end
    if isfield(sample,'video') && isscalar(sample.video), video_id=sample.video; end
end
fig=figure('Name',['Illumination diagnostics ' tag],'Color','w','Position',[100 100 1430 470]);
valid=isstruct(info) && isscalar(info) && isfield(info,'raw_common') && ...
    isfield(info,'gain') && ~isempty(info.raw_common) && numel(info.raw_common)==numel(info.gain);
if valid
    raw=double(info.raw_common(:)); gain=double(info.gain(:));
    valid=all(isfinite(raw)) && all(isfinite(gain)) && all(gain>0) && mean(raw)>0;
end
if ~valid
    for k=1:3
        subplot(1,3,k); axis off;
        text(0.5,0.5,{'未保存频谱共同光照诊断','median模式或诊断字段缺失'}, ...
            'Units','normalized','HorizontalAlignment','center');
    end
    sgtitle([tag '：光照诊断不可用'],'Interpreter','none');
    print(fig,fullfile(p.out.dir,['fig7_illumination_' tag '.png']),'-dpng','-r150');
    fprintf('光照诊断 %s：缺少频谱估计字段，已保存提示图。\n',tag); return;
end
F=numel(raw); t=[];
if isstruct(sample) && isfield(sample,'traj') && isstruct(sample.traj) && ...
        isfield(sample.traj,'t') && numel(sample.traj.t)==F
    t=double(sample.traj.t(:));
end
if isempty(t)
    if isfield(p,'sim') && isfield(p.sim,'fps') && isscalar(p.sim.fps) && p.sim.fps>0
        t=(0:F-1)'/p.sim.fps;
    else
        t=(0:F-1)';
    end
end
t_ms=(t-t(1))*1e3; first=t_ms<=100;
if ~any(first), first(1)=true; end
raw_norm=raw/mean(raw); gain_norm=gain/mean(gain);
accepted=[];
if isfield(info,'accepted_frequencies'), accepted=info.accepted_frequencies(:); end
applied=~isempty(accepted);
if isfield(info,'applied'), applied=logical(info.applied); end
if applied, status='已启动共同周期光照校正'; else, status='未启动校正：无可信共同周期分量'; end

%% 前100ms共同亮度与估计增益，真值仅作可选参考。
subplot(1,3,1); hold on;
plot(t_ms(first),raw_norm(first),'-','Color',[0.65 0.65 0.65],'DisplayName','观测共同亮度');
plot(t_ms(first),gain_norm(first),'-','Color',[0.05 0.35 0.78],'LineWidth',1.1,'DisplayName','估计增益');
if numel(truth)==F && all(isfinite(truth)) && mean(truth)>0
    truth=truth/mean(truth);
    plot(t_ms(first),truth(first),'--','Color',[0.86 0.24 0.12],'LineWidth',0.9, ...
        'DisplayName','生成器真增益（仅诊断）');
end
xlabel('相对首帧时间 (ms)'); ylabel('归一化共同亮度 / 增益');
title('前100 ms：观测与估计'); legend('Location','southoutside'); grid on; box on;

%% 共同亮度候选FFT谱；标记通过所有门控的频率。
subplot(1,3,2); hold on;
if isfield(info,'frequency_axis') && isfield(info,'spectrum') && ...
        ~isempty(info.frequency_axis) && numel(info.frequency_axis)==numel(info.spectrum)
    freq=double(info.frequency_axis(:)); spec=double(info.spectrum(:));
    ok=isfinite(freq) & isfinite(spec) & freq>=0 & spec>=0;
    freq=freq(ok); spec=spec(ok);
    if ~isempty(freq)
        floor_y=max(max(spec)*1e-6,1e-12); top_y=max(max(spec)*1.8,floor_y*10);
        semilogy(freq,max(spec,floor_y),'-','Color',[0.05 0.35 0.78],'DisplayName','共同亮度FFT幅值');
        for k=1:numel(accepted)
            plot([accepted(k) accepted(k)],[floor_y top_y],'--','Color',[0.86 0.24 0.12], ...
                'HandleVisibility','off');
        end
        set(gca,'YScale','log'); ylim([floor_y top_y]);
        xlim([0 max(freq)+eps]);
        if ~isempty(accepted)
            label=strjoin(arrayfun(@(x) sprintf('%.3f',x),accepted,'UniformOutput',false),', ');
            text(0.02,0.96,['采纳频率：' label ' Hz'],'Units','normalized','VerticalAlignment','top');
        else
            text(0.02,0.96,'未采纳任何频率','Units','normalized','VerticalAlignment','top');
        end
    else
        text(0.5,0.5,'频谱数据无有效数值','Units','normalized','HorizontalAlignment','center');
    end
else
    text(0.5,0.5,'未保存候选频谱','Units','normalized','HorizontalAlignment','center');
end
xlabel('频率 (Hz)'); ylabel('共同分量相对幅值'); title('候选频谱与采纳频率'); grid on; box on;

%% 全视频共同亮度波动；这里的RMS不是像素噪声或定位误差。
post=raw./gain; post=post/mean(post)-1; pre=raw_norm-1;
rms_pre=sqrt(mean(pre.^2))*100; rms_post=sqrt(mean(post.^2))*100;
subplot(1,3,3); hold on;
plot(t_ms,pre*100,'-','Color',[0.65 0.65 0.65],'DisplayName','校正前共同波动');
plot(t_ms,post*100,'-','Color',[0.05 0.35 0.78],'DisplayName','校正后共同波动');
xlabel('相对首帧时间 (ms)'); ylabel('共同亮度相对波动 (%)');
title(sprintf('共模RMS：%.4f%% → %.4f%%',rms_pre,rms_post));
legend('Location','southoutside'); grid on; box on;
sgtitle(sprintf('%s，视频%d：%s（不代表全部背景或散粒噪声）',tag,video_id,status),'Interpreter','none');
print(fig,fullfile(p.out.dir,['fig7_illumination_' tag '.png']),'-dpng','-r150');
end
