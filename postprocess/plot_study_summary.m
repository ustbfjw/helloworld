function record=plot_study_summary(summary_dir,image_dir)
%PLOT_STUDY_SUMMARY 只读完整CSV的独立后处理，禁止用于运行中部分结果。
%   record=plot_study_summary(summary_dir,image_dir)
%   输入summary.csv、pairwise.csv、protocol.json，不读取生成器或挑示例。
%   schema1:108个汇总/90个配对；schema2:72个汇总/54个配对。
%   PNG为300 dpi，PDF为矢量。合成fixture强制在文件名和标题标注。
%   不根据统计status自动声称验收或优越；缺CI不能当作零宽CI。
assert(nargin==2,'studyplot:Arguments','必须给出汇总目录和图片目录。');
summary_dir=char(string(summary_dir)); image_dir=char(string(image_dir));
assert(isfolder(summary_dir),'studyplot:MissingDirectory','汇总目录不存在。');
maxNumCompThreads(2);
names={'summary.csv','pairwise.csv','protocol.json'};
for k=1:numel(names)
    assert(isfile(fullfile(summary_dir,names{k})),'studyplot:MissingFile', ...
        '缺少完整输入：%s。',names{k});
end
protocol=jsondecode(fileread(fullfile(summary_dir,'protocol.json')));
S=readtable(fullfile(summary_dir,'summary.csv'),'TextType','string');
P=readtable(fullfile(summary_dir,'pairwise.csv'),'TextType','string');
[groups,conditions,methods]=validate_inputs(S,P,protocol);
protocol_hash=file_hash(fullfile(summary_dir,'protocol.json'));
hash_path=fullfile(summary_dir,'protocol.sha256');
if isfile(hash_path)
    assert(strcmp(protocol_hash,strtrim(fileread(hash_path))), ...
        'studyplot:ProtocolHash','protocol.json与保存的hash不符。');
end
fixture=isfield(protocol,'synthetic_statistics_fixture') && ...
    isequal(protocol.synthetic_statistics_fixture,true);
if ~fixture
    assert(~contains(lower(summary_dir),'fixture'),'studyplot:FixtureLabel', ...
        '输入目录含fixture，但协议没有明确合成标记。');
end
if ~isfolder(image_dir), mkdir(image_dir); end
prefix=sprintf('schema%d_',protocol.schema_version);
if fixture, prefix=['fixture_' prefix]; end
context=sprintf('schema %d；探索性视频 bootstrap B=%d；点态区间，未校正多重比较', ...
    protocol.schema_version,protocol.bootstrap_replicates);
if fixture, context=['FIXTURE 合成统计自测，非正式结果；' context]; end
colors=zeros(numel(methods),3); labels=strings(size(methods));
for m=1:numel(methods), [labels(m),colors(m,:)]=method_style(methods(m)); end
condition_labels=strings(size(conditions));
for c=1:numel(conditions), condition_labels(c)=condition_label(conditions(c)); end
artifacts=struct('stem',{},'png',{},'pdf',{},'log_axis',{},'display_floor_pct',{});

% 一、粒径偏差。CI端点单独画线，允许点估计在百分位区间之外。
fig=new_figure(); layout=tiledlayout(fig,3,4,'TileSpacing','compact','Padding','compact');
bias=100*(S.diameter_nm./S.diameter_true_nm-1);
lo=100*(S.d_ci95_low./S.diameter_true_nm-1);
hi=100*(S.d_ci95_high./S.diameter_true_nm-1);
all_values=[bias;lo;hi;-5;5]; all_values=all_values(isfinite(all_values));
padding=max(0.5,0.08*(max(all_values)-min(all_values)));
limits=[min(all_values)-padding,max(all_values)+padding];
for g=1:3
    ax=nexttile(layout,(g-1)*4+1,[1 3]); hold(ax,'on');
    first_group=g==1;
    for m=1:numel(methods)
        for c=1:numel(conditions)
            r=pick(S,groups(g),conditions(c),methods(m));
            x=c+method_offset(m,numel(methods));
            y=100*(r.diameter_nm/r.diameter_true_nm-1);
            if ~isfinite(y), continue; end
            has_ci=isfinite(r.d_ci95_low) && isfinite(r.d_ci95_high);
            if has_ci
                a=100*(r.d_ci95_low/r.diameter_true_nm-1);
                b=100*(r.d_ci95_high/r.diameter_true_nm-1);
                plot(ax,[x x],[a b],'-','Color',colors(m,:),'LineWidth',1.05,'HandleVisibility','off');
                plot(ax,x+[-.025 .025],[a a],'-','Color',colors(m,:),'HandleVisibility','off');
                plot(ax,x+[-.025 .025],[b b],'-','Color',colors(m,:),'HandleVisibility','off');
                plot(ax,x,y,'o','Color',colors(m,:),'MarkerFaceColor',colors(m,:), ...
                    'MarkerSize',4.5,'HandleVisibility','off');
            else
                plot(ax,x,y,'x','Color',colors(m,:),'LineWidth',1.8, ...
                    'MarkerSize',8,'HandleVisibility','off');
            end
        end
    end
    yline(ax,0,'-','Color',[.75 .75 .75],'HandleVisibility','off');
    yline(ax,-5,'--','Color',[.4 .4 .4],'HandleVisibility','off');
    yline(ax,5,'--','Color',[.4 .4 .4],'HandleVisibility','off');
    configure_axis(ax,condition_labels); ylim(ax,limits);
    ylabel(ax,'粒径相对误差（%）'); title(ax,group_title(S,groups(g)),'Interpreter','none');
    if first_group, add_method_legend(ax,labels,colors,true); end
    tx=nexttile(layout,g*4); draw_counts(tx,S,groups(g),conditions,methods,true);
end
title_and_note(layout,'粒径相对误差及95%区间',context, ...
    '虚线为±5%参考线，不自动判定验收；×为有限点估计但无CI；缺粒径不画成0，见右栏。C编号按横轴，M编号按图例。');
artifacts(end+1)=export_figure(fig,image_dir,[prefix 'diameter_error'],false,NaN);

% 二、逐像素照明增益残余RMS，不是相机图像总噪声。
fig=new_figure(); layout=tiledlayout(fig,3,4,'TileSpacing','compact','Padding','compact');
scale=rms_scale(S.gain_rms_pct);
for g=1:3
    ax=nexttile(layout,(g-1)*4+1,[1 3]); hold(ax,'on');
    for m=1:numel(methods)
        for c=1:numel(conditions)
            r=pick(S,groups(g),conditions(c),methods(m));
            rms_point(ax,c+method_offset(m,numel(methods)),r.gain_rms_pct,colors(m,:),scale);
        end
    end
    configure_axis(ax,condition_labels); apply_rms_scale(ax,scale);
    ylabel(ax,'全图照明增益残余 RMS（%）'); title(ax,group_title(S,groups(g)),'Interpreter','none');
    if g==1, add_method_legend(ax,labels,colors,false); end
    tx=nexttile(layout,g*4); draw_counts(tx,S,groups(g),conditions,methods,false);
end
title_and_note(layout,'全图照明误差',context, ...
    ['逐像素时间归一化残差的空间RMS，再跨视频RMS；不是总噪声或粒径误差。' scale.note ' C编号按横轴；M编号按图例。']);
artifacts(end+1)=export_figure(fig,image_dir,[prefix 'illumination_global'],scale.use_log,scale.floor);

% 三、全部预设局部条件的左右半区，不选择效果有利的区域或样本。
if protocol.schema_version==1, local_conditions="local_left_half";
else, local_conditions=["gaussian_patch","tilted_smooth","oblique_step"]; end
local_rows=ismember(S.condition,local_conditions);
scale=rms_scale([S.gain_rms_left_pct(local_rows);S.gain_rms_right_pct(local_rows)]);
fig=new_figure(); layout=tiledlayout(fig,3,4,'TileSpacing','compact','Padding','compact');
region_labels=strings(1,2*numel(local_conditions));
for c=1:numel(local_conditions)
    region_labels(2*c-1)=condition_label(local_conditions(c))+"（左半区）";
    region_labels(2*c)=condition_label(local_conditions(c))+"（右半区）";
end
for g=1:3
    ax=nexttile(layout,(g-1)*4+1,[1 3]); hold(ax,'on');
    for m=1:numel(methods)
        for c=1:numel(local_conditions)
            r=pick(S,groups(g),local_conditions(c),methods(m));
            y=[r.gain_rms_left_pct,r.gain_rms_right_pct];
            for side=1:2
                rms_point(ax,2*c-2+side+method_offset(m,numel(methods)),y(side),colors(m,:),scale);
            end
        end
    end
    configure_axis(ax,region_labels); apply_rms_scale(ax,scale);
    ylabel(ax,'区域照明增益残余 RMS（%）'); title(ax,group_title(S,groups(g)),'Interpreter','none');
    if g==1, add_method_legend(ax,labels,colors,false); end
    tx=nexttile(layout,g*4); draw_counts(tx,S,groups(g),local_conditions,methods,false);
end
if protocol.schema_version==1
    local_note='左半图是真值受扰区域，右半图未注入动态照明；右部残差可显示人工闪烁。';
else
    local_note='左右为固定几何半区，并非各掩膜的纯受扰/未受扰区；包含全部三个局部条件。';
end
title_and_note(layout,'局部照明条件：左右区域对照',context, ...
    [local_note ' 此量不是相机总噪声。' scale.note ' 右栏C编号按局部条件顺序，左右共用一次计数。']);
artifacts(end+1)=export_figure(fig,image_dir,[prefix 'illumination_regions'],scale.use_log,scale.floor);

% 四、接受次数只表示模型选择，不等于性能或通过验收。
if protocol.schema_version==2
    fig=new_figure(); layout=tiledlayout(fig,3,4,'TileSpacing','compact','Padding','compact');
    [~,candidate_color]=method_style("spatial_validated");
    for g=1:3
        ax=nexttile(layout,(g-1)*4+1,[1 3]); hold(ax,'on');
        a=zeros(1,numel(conditions)); n=a; known=a;
        for c=1:numel(conditions)
            r=pick(S,groups(g),conditions(c),"spatial_validated");
            a(c)=r.spatial_accepted_videos; n(c)=r.n_videos; known(c)=r.spatial_decision_videos;
        end
        bar(ax,1:numel(conditions),a,.55,'FaceColor',candidate_color);
        for c=1:numel(conditions)
            text(ax,c,a(c)+.035*n(c),sprintf('%d/%d',a(c),n(c)), ...
                'HorizontalAlignment','center','FontSize',10,'Interpreter','none');
        end
        configure_axis(ax,condition_labels); ylim(ax,[0,1.2*max(n)]);
        ylabel(ax,'空间模型接受视频数'); title(ax,group_title(S,groups(g)),'Interpreter','none');
        tx=nexttile(layout,g*4); draw_counts(tx,S,groups(g),conditions,"spatial_validated",false);
        decisions=string(known)+"/"+string(n);
        text(tx,0,.20,{'门控记录数（C1至C6）：', ...
            char(strjoin(decisions(1:3),', ')),char(strjoin(decisions(4:6),', '))}, ...
            'Units','normalized','VerticalAlignment','top','Interpreter','none','FontSize',9);
    end
    title_and_note(layout,'空间模型接受次数 / 全部视频数',context, ...
        '读取is_spatial汇总计数，不把applied全局校正当空间接受；包含无定位视频。次数不是性能结论；右栏M1仅指空间留出候选。');
    artifacts(end+1)=export_figure(fig,image_dir,[prefix 'spatial_acceptance'],false,NaN);
end

record=struct('schema_version',protocol.schema_version,'fixture',fixture, ...
    'summary_directory',summary_dir,'image_directory',image_dir,'summary_rows',height(S), ...
    'pairwise_rows',height(P),'protocol_sha256',protocol_hash, ...
    'summary_sha256',file_hash(fullfile(summary_dir,'summary.csv')), ...
    'pairwise_sha256',file_hash(fullfile(summary_dir,'pairwise.csv')), ...
    'helper_sha256',file_hash([mfilename('fullpath') '.m']), ...
    'matlab_version',version,'generated_at',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss Z')), ...
    'bootstrap_replicates',protocol.bootstrap_replicates, ...
    'missing_diameter_points',sum(~isfinite(S.diameter_nm)), ...
    'finite_diameter_without_ci',sum(isfinite(S.diameter_nm) & ~isfinite(S.d_ci95_low)), ...
    'zero_global_illumination_rows',sum(S.gain_rms_pct==0),'artifacts',artifacts, ...
    'interpretation','CSV-only exploratory plots; no automatic performance or acceptance conclusion');
fid=fopen(fullfile(image_dir,[prefix 'plot_provenance.json']),'w','n','UTF-8');
assert(fid>=0,'studyplot:WriteFailed','无法写入绘图来源记录。'); cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(record,'PrettyPrint',true));
fprintf('绘图完成：schema%d，fixture=%d，%d张PNG/PDF；CSV保持只读。\n', ...
    protocol.schema_version,fixture,numel(artifacts));
end

function [groups,conditions,methods]=validate_inputs(S,P,protocol)
required={'schema_version','groups','conditions','methods','bootstrap_replicates'};
assert(isstruct(protocol) && isscalar(protocol) && all(isfield(protocol,required)), ...
    'studyplot:Protocol','协议字段缺失。');
assert(ismember(protocol.schema_version,[1 2]) && numel(protocol.groups)==3 && ...
    numel(protocol.conditions)==6 && protocol.bootstrap_replicates==2000, ...
    'studyplot:Protocol','只支持已冻结的两种完整协议。');
groups=string({protocol.groups.name}); conditions=string({protocol.conditions.name});
methods=reshape(string(protocol.methods),1,[]);
if protocol.schema_version==1
    expected_methods=["median","pn","full","no_spatial_gate","no_magnetic_guard","no_refinement"];
    expected_conditions=["static","global_offbin","near_inside","near_outside","same_frequency","local_left_half"];
    primary="full";
else
    expected_methods=["median","pn","full","spatial_validated"];
    expected_conditions=["static","global_offbin","same_frequency","gaussian_patch","tilted_smooth","oblique_step"];
    primary="spatial_validated";
end
assert(isequal(methods,expected_methods) && isequal(conditions,expected_conditions) && ...
    isequal(groups,["mono30","mono100","mono200"]) && ...
    isequal([protocol.groups.diameter_nm],[30 100 200]) && ...
    isequal([protocol.groups.n_videos],[20 60 60]), ...
    'studyplot:Protocol','方法、条件、组别或样本数与固定协议不符。');
requiredS={'group','condition','method','diameter_true_nm','n_videos','active_videos', ...
    'diameter_nm','diameter_bias_pp','d_ci95_low','d_ci95_high', ...
    'gain_rms_pct','gain_rms_left_pct','gain_rms_right_pct','bootstrap_replicates'};
requiredP={'group','condition','method','reference','n_videos','bootstrap_replicates'};
assert(all(ismember(requiredS,S.Properties.VariableNames)) && ...
    all(ismember(requiredP,P.Properties.VariableNames)),'studyplot:Columns','CSV必需列缺失。');
if protocol.schema_version==2
    assert(all(ismember({'spatial_accepted_videos','spatial_decision_videos'},S.Properties.VariableNames)), ...
        'studyplot:Columns','缺少空间接受/决策次数。');
end
assert(height(S)==3*6*numel(methods) && height(P)==3*6*(numel(methods)-1), ...
    'studyplot:RowCount','CSV行数不完整。');
skey=S.group+"|"+S.condition+"|"+S.method;
pkey=P.group+"|"+P.condition+"|"+P.method+"|"+P.reference;
assert(numel(unique(skey))==height(S) && numel(unique(pkey))==height(P), ...
    'studyplot:Duplicate','CSV存在重复组合。');
assert(all(S.bootstrap_replicates==protocol.bootstrap_replicates) && ...
    all(P.bootstrap_replicates==protocol.bootstrap_replicates), ...
    'studyplot:Bootstrap','重抽样次数与协议不符。');
for g=1:3
    for c=1:6
        for m=1:numel(methods)
            r=pick(S,groups(g),conditions(c),methods(m));
            assert(r.n_videos==protocol.groups(g).n_videos && ...
                r.diameter_true_nm==protocol.groups(g).diameter_nm && ...
                isfinite(r.active_videos) && r.active_videos>=0 && r.active_videos<=r.n_videos && ...
                r.active_videos==round(r.active_videos),'studyplot:Counts','组标尺或视频计数不符。');
            assert(~isinf(r.diameter_nm) && (isnan(r.diameter_nm) || r.diameter_nm>0), ...
                'studyplot:Size','尺寸必须是正数或明确NaN。');
            if isfinite(r.diameter_nm)
                assert(isfinite(r.diameter_bias_pp) && abs(r.diameter_bias_pp- ...
                    100*(r.diameter_nm/r.diameter_true_nm-1))<1e-8, ...
                    'studyplot:Bias','CSV尺寸与偏差列不一致。');
            end
            interval=[r.d_ci95_low r.d_ci95_high];
            assert(all(isnan(interval)) || (all(isfinite(interval)) && all(interval>0) && ...
                interval(1)<=interval(2) && isfinite(r.diameter_nm)), ...
                'studyplot:Interval','CI须为完整有序正端点或双NaN。');
            y=[r.gain_rms_pct,r.gain_rms_left_pct,r.gain_rms_right_pct];
            assert(all(isfinite(y) & y>=0),'studyplot:RMS','照明误差必须有限非负。');
            if protocol.schema_version==2
                a=r.spatial_accepted_videos; k=r.spatial_decision_videos;
                assert(all(isfinite([a k])) && 0<=a && a<=k && k<=r.n_videos && ...
                    all([a k]==round([a k])),'studyplot:DecisionCount','空间次数无效。');
            end
            if methods(m)~=primary
                q=P(P.group==groups(g) & P.condition==conditions(c) & ...
                    P.method==primary & P.reference==methods(m),:);
                assert(height(q)==1 && q.n_videos==r.n_videos, ...
                    'studyplot:Pairing','完整配对组合或视频数不符。');
            end
        end
    end
end
end

function r=pick(S,g,c,m)
r=S(S.group==g & S.condition==c & S.method==m,:);
assert(height(r)==1,'studyplot:Combination','汇总组合缺失或重复。');
end

function [name,color]=method_style(method)
switch char(method)
    case 'median', name="中值背景"; color=[.35 .35 .35];
    case 'pn', name="全帧功率归一化 PN"; color=[.90 .50 .05];
    case 'full', name="完整全局校正 full"; color=[.00 .35 .70];
    case 'no_spatial_gate', name="去空间一致性门控"; color=[.60 .35 .70];
    case 'no_magnetic_guard', name="去磁频保护"; color=[.80 .20 .15];
    case 'no_refinement', name="去格点内细化"; color=[.15 .55 .25];
    case 'spatial_validated', name="空间留出候选"; color=[.00 .60 .60];
    otherwise, error('studyplot:Method','未知方法。');
end
end

function label=condition_label(condition)
switch char(condition)
    case 'static', label="静态";
    case 'global_offbin', label="全局423.7 Hz";
    case 'near_inside', label="保护带内11.4 Hz";
    case 'near_outside', label="保护带外13.4 Hz";
    case 'same_frequency', label="同频10 Hz";
    case 'local_left_half', label="左半图423.7 Hz";
    case 'gaussian_patch', label="高斯局部斑";
    case 'tilted_smooth', label="平滑非均匀场";
    case 'oblique_step', label="斜向阶跃";
    otherwise, error('studyplot:Condition','未知条件。');
end
end

function fig=new_figure()
fig=figure('Visible','off','Color','w','Position',[50 50 1700 1250], ...
    'PaperPositionMode','auto','Renderer','painters');
end

function x=method_offset(m,n)
if n==1, x=0; else, x=-.30+.60*(m-1)/(n-1); end
end

function configure_axis(ax,labels)
xlim(ax,[.5 numel(labels)+.5]); xticks(ax,1:numel(labels)); xticklabels(ax,labels);
ax.TickLabelInterpreter='none'; ax.FontSize=10;
grid(ax,'on'); box(ax,'on'); ax.GridAlpha=.16;
end

function value=group_title(S,group)
rows=S(S.group==group,:);
value=sprintf('%g nm组；全部独立视频 n=%d；active见右侧矩阵', ...
    rows.diameter_true_nm(1),rows.n_videos(1));
end

function add_method_legend(ax,names,colors,show_missing)
h=gobjects(1,numel(names)); labels=strings(1,numel(names));
for m=1:numel(names)
    h(m)=plot(ax,NaN,NaN,'o','Color',colors(m,:),'MarkerFaceColor',colors(m,:));
    labels(m)=sprintf('M%d ',m)+names(m);
end
if show_missing
    h(end+1)=plot(ax,NaN,NaN,'kx','LineWidth',1.5); labels(end+1)="无CI";
end
legend(ax,h,labels,'Location','northoutside','Orientation','horizontal', ...
    'NumColumns',min(4,numel(h)),'Interpreter','none','FontSize',10);
end

function draw_counts(ax,S,group,conditions,methods,show_missing)
axis(ax,'off');
text(ax,0,.96,'active / n（有有效位移对的视频）','Units','normalized', ...
    'VerticalAlignment','top','FontSize',9,'Interpreter','none');
% 每个表格单元单独定位，避免六列计数因比例字体或位数变化溢出图幅。
xs=linspace(.15,.78,max(2,numel(conditions))); xs=xs(1:numel(conditions));
for c=1:numel(conditions)
    text(ax,xs(c),.85,sprintf('C%d',c),'Units','normalized', ...
        'HorizontalAlignment','center','FontSize',8.5,'Interpreter','none');
end
for m=1:numel(methods)
    y=.76-(m-1)*.09;
    text(ax,0,y,sprintf('M%d',m),'Units','normalized','FontSize',8.5,'Interpreter','none');
    for c=1:numel(conditions)
        r=pick(S,group,conditions(c),methods(m));
        value=sprintf('%d/%d',r.active_videos,r.n_videos);
        if show_missing && ~isfinite(r.diameter_nm), value=[value '?'];
        elseif show_missing && ~isfinite(r.d_ci95_low), value=[value '*']; end
        text(ax,xs(c),y,value,'Units','normalized','HorizontalAlignment','center', ...
            'FontSize',8.5,'Interpreter','none');
    end
end
if show_missing
    text(ax,0,.18,{'*：有限粒径但无CI（图中×）','?：缺粒径（不按0绘图）'}, ...
        'Units','normalized','VerticalAlignment','top','FontSize',9,'Interpreter','none');
end
end

function scale=rms_scale(values)
positive=values(values>0); scale=struct('use_log',false,'floor',NaN,'limits',[0 1],'note','');
if ~isempty(positive) && max(positive)/min(positive)>=100
    scale.use_log=true; scale.floor=10^(floor(log10(min(positive)))-1);
    if scale.floor==0, scale.floor=min(positive)/2; end
    scale.limits=[scale.floor/2,2*max(positive)];
    scale.note=sprintf('对数轴显示下限%.3g%%；真实0用▽放在该下限并标0，未改写数值。',scale.floor);
else
    scale.limits=[0,max(.01,1.18*max(values))];
    scale.note='线性轴；数字0明确表示真实零值。';
end
end

function rms_point(ax,x,value,color,scale)
if value==0
    y=0; marker='o'; if scale.use_log, y=scale.floor; marker='v'; end
    plot(ax,x,y,marker,'Color',color,'MarkerFaceColor','w','MarkerSize',5,'HandleVisibility','off');
    text(ax,x,y,'0','Color',color,'VerticalAlignment','bottom','HorizontalAlignment','center', ...
        'FontSize',8,'Interpreter','none','Clipping','off');
else
    plot(ax,x,value,'o','Color',color,'MarkerFaceColor',color,'MarkerSize',5,'HandleVisibility','off');
end
end

function apply_rms_scale(ax,scale)
if scale.use_log, ax.YScale='log'; end
ylim(ax,scale.limits);
end

function title_and_note(layout,title_text,context,note)
sgtitle(layout,{title_text,context},'Interpreter','none','FontSize',14);
xlabel(layout,note,'Interpreter','none','FontSize',10);
end

function artifact=export_figure(fig,outdir,stem,use_log,display_floor)
cleanup=onCleanup(@()close(fig)); %#ok<NASGU>
set(findall(fig,'-property','FontName'),'FontName','Microsoft YaHei');
png=fullfile(outdir,[stem '.png']); pdf=fullfile(outdir,[stem '.pdf']);
exportgraphics(fig,png,'Resolution',300);
exportgraphics(fig,pdf,'ContentType','vector','BackgroundColor','white');
info=imfinfo(png); assert(info.Width>1000 && info.Height>1000,'studyplot:Image','PNG输出尺寸异常。');
assert(isfile(pdf),'studyplot:PDF','PDF未生成。');
artifact=struct('stem',stem,'png',png,'pdf',pdf,'log_axis',use_log,'display_floor_pct',display_floor);
end

function hash=file_hash(path)
fid=fopen(path,'rb'); assert(fid>=0,'studyplot:ReadFailed','文件无法读取。');
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid), md.update(fread(fid,1024*1024,'*uint8')); end
hash=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
