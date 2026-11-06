function manifest=create_candidate_examples(source_root,results_dir,outdir)
%CREATE_CANDIDATE_EXAMPLES 冻结验收之后的单视频复现实例，不改变正式统计。
% 输入显式指定冻结源码、已完成结果及新的独立输出目录。
% 先验选择：每组最早tilted_smooth合格且局部接纳的视频；否则最早任意
% 条件有合格轨迹的视频（标明无接受的局部示例）；全无则跳过。
% 选择不读取真值、照明误差或粒径结果。全部检查点齐全才允许重建。
% 只跑spatial_validated生产路径，与检查点充分统计量逐字段核验。
% fixture_only=true的冻结协议允许小于5000帧，仅用于接口测试。
source_root=canonical_path(source_root); results_dir=canonical_path(results_dir);
outdir=canonical_path(outdir);
assert(isfolder(source_root) && isfolder(results_dir),'examples:InvalidPath','源码或结果目录不存在。');
assert(~isfolder(outdir) && ~isfile(outdir),'examples:OutputExists','请使用新的独立示例输出目录。');
assert(~inside(outdir,source_root) && ~inside(outdir,results_dir), ...
    'examples:UnsafeOutput','示例目录不能位于冻结源码或正式结果目录内。');
oldpath=path; oldrng=rng; oldvisible=get(groot,'DefaultFigureVisible');
cleanup=onCleanup(@() restore_environment(oldpath,oldrng,oldvisible));
addpath(source_root,'-begin'); setup_paths(source_root); maxNumCompThreads(2);
set(groot,'DefaultFigureVisible','off');
assert(strcmpi(which('inta_track_video'),fullfile(source_root,'03_inta_analysis','inta_track_video.m')), ...
    'examples:WrongSource','生产追踪函数未解析到指定源码。');
protocol_path=fullfile(results_dir,'protocol.json');
protocol=jsondecode(fileread(protocol_path));
protocol_hash=strtrim(fileread(fullfile(results_dir,'protocol.sha256')));
assert(strcmp(targeted_sha256(protocol_path),protocol_hash),'examples:ProtocolMismatch','协议哈希不一致。');
assert(protocol.schema_version==2 && numel(protocol.groups)==3 && numel(protocol.conditions)==6, ...
    'examples:InvalidProtocol','本helper只处理候选schema 2的三组六条件。');
fixture=isfield(protocol,'fixture_only') && isequal(protocol.fixture_only,true);
assert(fixture || protocol.frames==5000,'examples:InvalidProtocol','正式示例必须来自5000帧协议。');
assert(strcmp(targeted_sha256(fullfile(results_dir,'parameters.mat')),protocol.parameters_sha256), ...
    'examples:ParameterMismatch','冻结参数哈希不一致。');
for k=1:numel(protocol.source_files)
    assert(strcmp(targeted_sha256(fullfile(source_root,protocol.source_files{k})),protocol.source_hashes{k}), ...
        'examples:SourceMismatch','冻结源码不一致：%s',protocol.source_files{k});
end
frozen=load(fullfile(results_dir,'parameters.mat'),'parameters');
method=find(strcmp(protocol.methods,'spatial_validated'));
tilted=find(strcmp({protocol.conditions.name},'tilted_smooth'));
assert(isscalar(method) && isscalar(tilted),'examples:InvalidProtocol','候选方法或选择条件缺失/重复。');

% 全部组的完整性检查先于任何重建或输出，不读取metric或粒径数值作选择。
chosen=nan(1,3); selection=cell(1,3);
for g=1:3
    first_active=nan; first_local=nan;
    for v=1:protocol.groups(g).n_videos
        checkpoint=fullfile(results_dir,protocol.groups(g).name,sprintf('video_%03d.mat',v));
        assert(isfile(checkpoint),'examples:IncompleteResults','结果尚未齐全：%s',checkpoint);
        old=load(checkpoint,'items','p','seed','video_index','protocol_hash');
        expected_seed=protocol.groups(g).seed_base+protocol.seed_stride*(v-1);
        assert(strcmp(old.protocol_hash,protocol_hash) && old.seed==expected_seed && old.video_index==v && ...
            isequaln(old.p,frozen.parameters{g}) && old.p.sim.n_frames==protocol.frames && ...
            isequal(size(old.items),[6,numel(protocol.methods)]), ...
            'examples:CheckpointMismatch','检查点与冻结协议不一致：%s',checkpoint);
        active=arrayfun(@(x) x.stats.n_tracks>0,old.items(:,method));
        if isnan(first_active) && any(active), first_active=v; end
        if isnan(first_local) && active(tilted) && ...
                isequal(old.items(tilted,method).illumination_info.is_spatial,true)
            first_local=v;
        end
    end
    if isfinite(first_local)
        chosen(g)=first_local; selection{g}='最早tilted_smooth合格且局部接纳的视频';
    elseif isfinite(first_active)
        chosen(g)=first_active; selection{g}='无接受的局部示例；使用最早任一条件有合格轨迹的视频';
    else
        selection{g}='跳过：六条件均无合格轨迹';
    end
end

mkdir(outdir);
manifest=struct('schema_version',1,'status','in_progress','fixture_only',fixture, ...
    'source_root',source_root,'results_dir',results_dir,'protocol_sha256',protocol_hash, ...
    'helper_sha256',targeted_sha256([mfilename('fullpath') '.m']), ...
    'full_frame_count',protocol.frames,'max_num_threads',2,'created_at', ...
    char(datetime('now','Format','yyyy-MM-dd HH:mm:ss Z')), ...
    'selection_rule','Per group, lowest video_index with tilted_smooth candidate n_tracks>0 and is_spatial=true; otherwise lowest video_index with candidate n_tracks>0 in any condition; otherwise skip. No truth or sizing-error selection.', ...
    'scope','Single-video visualization only; not population performance or acceptance statistics.', ...
    'floating_comparison','Each finite numeric field: max(abs(actual-reference)) <= 1e-12 * max(abs(actual),abs(reference)) across that field; integer fields exact; NaN/Inf locations and signs exact.', ...
    'truth_use','Generator truth is reference/diagnostic only; production analysis uses it only for existing post-inference matching diagnostics.', ...
    'groups',[]);
group_template=struct('name','','status','','selection_reason','','video_index',nan,'seed',nan, ...
    'checkpoint_sha256','','accepted_local_example',false,'conditions',[]);
manifest.groups=repmat(group_template,1,3);
for g=1:3
    name=protocol.groups(g).name; manifest.groups(g).name=name;
    manifest.groups(g).selection_reason=selection{g}; manifest.groups(g).video_index=chosen(g);
    if isnan(chosen(g))
        manifest.groups(g).status='skipped_no_qualified_tracks'; continue;
    end
    v=chosen(g); checkpoint=fullfile(results_dir,name,sprintf('video_%03d.mat',v));
    old=load(checkpoint,'items','p','seed'); p=old.p; seed=old.seed;
    manifest.groups(g).seed=seed;
    manifest.groups(g).checkpoint_sha256=targeted_sha256(checkpoint);
    manifest.groups(g).accepted_local_example=old.items(tilted,method).stats.n_tracks>0 && ...
        old.items(tilted,method).illumination_info.is_spatial;
    fprintf('事后示例 %s video=%d seed=%d：%s\n',name,v,seed,selection{g});
    tab=ipsf_lookup_table(p); traj=simulate_trajectories(p,seed);
    [videos,metas]=render_candidate_conditions(traj,p,tab,protocol.conditions);
    pm=p; pm.ana.background_mode='spatial_validated'; pm.ana.bootstrap_n=0;
    samples=cell(1,6); records=cell(1,6);
    for c=1:6
        [tracks,ex]=inta_track_video(videos{c},pm);
        tracks=analyze_tracks(tracks,traj,pm,v);
        res=inta_population_analysis(tracks,pm);
        stats=targeted_sufficient_stats(tracks,res,pm);
        verification=compare_stats(stats,old.items(c,method).stats);
        assert(isequal(ex.illumination.is_spatial,old.items(c,method).illumination_info.is_spatial), ...
            'examples:ModelMismatch','%s/%s的实际空间接纳路径不一致。',name,protocol.conditions(c).name);
        assert(isempty(ex.frame) || isequal(ex.raw,videos{c}(:,:,ex.frame)), ...
            'examples:RawMismatch','示例raw不是原始相机帧。');
        samples{c}=struct('ex',ex,'tracks',tracks,'res',res,'stats',stats,'meta',metas{c});
        records{c}=struct('name',protocol.conditions(c).name,'qualified_tracks',stats.n_tracks, ...
            'selected_localizations',stats.n_loc,'linked_tracks',numel(tracks), ...
            'is_spatial',ex.illumination.is_spatial,'illumination_reason',ex.illumination.reason, ...
            'statistics_verification',verification,'example_frame',ex.frame, ...
            'excerpt_frame_indices',[],'frame_note','','files',{{}});
    end

    % 六条件全部重算一致后才保存本组示例；不以局部通过冒充全组复现成功。
    groupdir=fullfile(outdir,name); mkdir(groupdir);
    pp=p; pp.out.dir=groupdir;
    plot_trajectories_3d(traj,pp,[name '_truth_reference']);
    fig=gcf; sgtitle(sprintf('%s 视频%d：生成器三维轨迹参考；单视频实例，非测量结果',name,v),'Interpreter','none');
    print(fig,fullfile(groupdir,['fig1_trajectories_3d_' name '_truth_reference.png']),'-dpng','-r150'); close(fig);
    save(fullfile(groupdir,'trajectory_reference.mat'),'traj','p','seed','v','protocol_hash','-v7.3');
    for c=1:6
        condition=protocol.conditions(c); conditiondir=fullfile(groupdir,condition.name); mkdir(conditiondir);
        sample=samples{c}; ex=sample.ex; tracks=sample.tracks; res=sample.res; stats=sample.stats; meta=sample.meta;
        model=ex.illumination_model; info=ex.illumination;
        pp=pm; pp.out.dir=conditiondir; tag=[name '_' condition.name];
        if isempty(ex.frame)
            center=ceil(protocol.frames/2);
            records{c}.frame_note='无链接轨迹，生产ex为空；仅展示固定中间raw帧，不补造C/RVT/检测位置。';
            plot_empty_frame(videos{c}(:,:,center),pp,tag,center);
        else
            center=ex.frame;
            if stats.n_tracks==0
                records{c}.frame_note='无合格轨迹；C/RVT及检测圈来自生产示例，不能解释为合格粒径轨迹。';
            else
                records{c}.frame_note='生产路径最长链接轨迹的中间帧；检测圈不表示该轨迹必然入选粒径统计。';
            end
            plot_example_frames(ex,pp,tag); fig=gcf;
            sgtitle(sprintf('%s：单视频实例，合格轨迹%d；检测圈为本帧检测位置',tag,stats.n_tracks),'Interpreter','none');
            print(fig,fullfile(conditiondir,['fig2_example_frames_' tag '.png']),'-dpng','-r150'); close(fig);
        end
        plot_sample=struct('ex',ex,'meta',meta,'traj',traj,'video',v);
        plot_illumination_diagnostics(plot_sample,pp,tag); close(gcf);
        plot_positions(tracks,res,traj,pp,tag);
        frame_indices=max(1,center-100):min(protocol.frames,center+100);
        raw_excerpt=videos{c}(:,:,frame_indices); full_frame_count=protocol.frames; is_excerpt=true;
        display_only=true; fixture_only=fixture;
        assert(isa(raw_excerpt,'uint16') && size(raw_excerpt,3)<=201);
        save(fullfile(conditiondir,'example.mat'),'ex','model','info','meta','tracks','res','stats', ...
            'pp','seed','v','condition','protocol_hash','fixture_only','-v7.3');
        save(fullfile(conditiondir,'raw_excerpt_display_only.mat'),'raw_excerpt','frame_indices', ...
            'full_frame_count','is_excerpt','display_only','fixture_only','seed','v','condition','protocol_hash','-v7.3');
        records{c}.excerpt_frame_indices=frame_indices;
        files=dir(conditiondir); files=files(~[files.isdir]); records{c}.files={files.name};
    end
    manifest.groups(g).conditions=[records{:}]; manifest.groups(g).status='six_conditions_reproduced';
    clear videos metas traj samples records tab
end
manifest.status='complete';
manifest.reproduction_note='All generated examples passed checkpoint-statistics comparison. Skipped groups contain no fabricated trajectories. These examples do not replace full-study statistics.';
write_text(fullfile(outdir,'example_manifest.json'),jsonencode(manifest,PrettyPrint=true));
write_text(fullfile(outdir,'README.txt'),sprintf([ ...
    '这些图和MAT是按预先说明规则选择的单视频实例，不是群体性能证据。\n' ...
    '三组独立选择视频，每组六条件共用同一seed和真实轨迹。\n' ...
    'example.mat保留生产ex、实际model/info及meta原诊断；真值仅为参考。\n' ...
    'raw_excerpt_display_only.mat最多201帧，is_excerpt=true，不能替代完整视频验收统计。\n' ...
    '无轨迹/回退/跳过状态与逐字段数值核验见example_manifest.json。\n' ...
    'fixture_only=%d；完整帧数=%d；protocol SHA-256=%s\n'],fixture,protocol.frames,protocol_hash));
fprintf('示例生成结束（仅展示）：%s\n',outdir);
end

function verification=compare_stats(actual,reference)
assert(isequal(sort(fieldnames(actual)),sort(fieldnames(reference))), ...
    'examples:StatisticsMismatch','统计字段集合不一致。');
integer_fields={'counts','n_tracks','n_loc','unmatched_selected','matched_localizations'};
fields=fieldnames(reference); checks=repmat(struct('field','','exact',false,'max_abs_difference',0,'relative_difference',0),1,numel(fields));
for k=1:numel(fields)
    name=fields{k}; a=actual.(name); b=reference.(name);
    assert(isnumeric(a) && isnumeric(b) && isequal(size(a),size(b)), ...
        'examples:StatisticsMismatch','统计字段类型/形状不一致：%s',name);
    assert(isequal(isnan(a),isnan(b)) && isequal(isinf(a),isinf(b)) && ...
        isequal(a(isinf(a)),b(isinf(b))), 'examples:StatisticsMismatch','统计非有限值不一致：%s',name);
    exact=isequaln(a,b); keep=isfinite(a); delta=0; relative=0;
    if any(keep(:))
        av=double(a(keep)); bv=double(b(keep)); av=av(:); bv=bv(:); delta=max(abs(av-bv));
        scale=max([abs(av);abs(bv)]); relative=delta/max(scale,realmin);
    end
    if ismember(name,integer_fields), ok=exact; else, ok=relative<=1e-12; end
    assert(ok,'examples:StatisticsMismatch','统计不一致：%s，maxdiff=%.17g，relative=%.17g',name,delta,relative);
    checks(k)=struct('field',name,'exact',exact,'max_abs_difference',delta,'relative_difference',relative);
end
verification=struct('passed',true,'integer_fields_exact',true,'floating_relative_tolerance',1e-12, ...
    'all_fields_exact',all([checks.exact]),'fields',checks);
end

function plot_empty_frame(raw,p,tag,frame)
fig=figure('Color','w','Position',[80 80 1250 400]);
subplot(1,3,1); ax=((1:p.cam.Npx)-(p.cam.Npx+1)/2)*p.cam.px*1e6;
imagesc(ax,ax,double(raw)); axis image; colormap(gray); colorbar;
title(sprintf('原始相机帧%d（ADU）',frame)); xlabel('x (\mum)'); ylabel('y (\mum)');
for k=2:3
    subplot(1,3,k); axis off; text(.5,.5,{'无链接轨迹，生产示例为空', ...
        '未补造对比度/RVT/检测位置'},'Units','normalized','HorizontalAlignment','center');
end
sgtitle([tag '：无轨迹单视频实例'],'Interpreter','none');
print(fig,fullfile(p.out.dir,['fig2_example_frames_' tag '.png']),'-dpng','-r150'); close(fig);
end

function plot_positions(tracks,res,traj,p,tag)
fig=figure('Color','w','Position',[80 80 1100 500]);
subplot(1,2,1); hold on;
if any(res.sel)
    for k=find(res.sel), plot(tracks(k).x_um,tracks(k).y_um,'-'); end
    title(sprintf('推断中入选的二维轨迹：%d条',sum(res.sel)));
else
    text(.5,.5,'无合格轨迹；不补造轨迹','Units','normalized','HorizontalAlignment','center');
    title('无合格二维轨迹');
end
axis equal; grid on; xlabel('x (\mum)'); ylabel('y (\mum)');
subplot(1,2,2); hold on;
for k=1:numel(traj.d)
    xy=reshape(traj.pos(:,k,1:2),[],2)*1e6; plot(xy(:,1),xy(:,2),'--');
end
axis equal; grid on; xlabel('x (\mum)'); ylabel('y (\mum)'); title('生成器真实轨迹xy投影（仅参考）');
sgtitle([tag '：单视频轨迹实例，非群体性能'],'Interpreter','none');
print(fig,fullfile(p.out.dir,'tracks_2d_example.png'),'-dpng','-r150'); close(fig);
end

function value=canonical_path(value)
assert((ischar(value) && isrow(value)) || (isstring(value) && isscalar(value)), ...
    'examples:InvalidPath','路径必须是字符行或标量字符串。');
value=char(java.io.File(char(value)).getCanonicalPath());
end

function yes=inside(child,parent)
yes=strcmpi(child,parent) || startsWith(lower([child filesep]),lower([parent filesep]));
end

function restore_environment(oldpath,oldrng,oldvisible)
path(oldpath); rng(oldrng); set(groot,'DefaultFigureVisible',oldvisible);
end

function write_text(file,value)
fid=fopen(file,'w','n','UTF-8'); assert(fid>=0,'examples:WriteFailure','无法写入%s',file);
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',value);
end
