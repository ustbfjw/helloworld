function report=summarize_candidate_experiment(outdir)
%SUMMARIZE_CANDIDATE_EXPERIMENT 按独立视频配对汇总空间候选，不以效果决定通过。
% 每次重抽样重新合并位移矩再估计D、d，所有条件/方法共用视频索引。
% 空视频保持在抽样总体中；真值仅用于推断完成后的误差评价。
protocol=jsondecode(fileread(fullfile(outdir,'protocol.json')));
protocol_hash=strtrim(fileread(fullfile(outdir,'protocol.sha256')));
assert(strcmp(targeted_sha256(fullfile(outdir,'protocol.json')),protocol_hash), ...
    'candidate:ProtocolHashMismatch','冻结协议内容与SHA256不一致。');
methods=cellstr(string(protocol.methods)); conditions={protocol.conditions.name};
assert(protocol.schema_version==2 && numel(conditions)==6 && numel(methods)==4);
assert(isequal(methods(:)',{'median','pn','full','spatial_validated'}));
assert(isequal(conditions,{'static','global_offbin','same_frequency','gaussian_patch','tilted_smooth','oblique_step'}));
assert(numel(protocol.groups)==3 && protocol.bootstrap_replicates==2000);
B=protocol.bootstrap_replicates; nc=numel(conditions); nm=numel(methods);
candidate_index=find(strcmp(methods,'spatial_validated')); rows=struct([]); comparisons=struct([]);
video_rows=struct([]); bootstrap=cell(1,3);
for ig=1:numel(protocol.groups)
    group=protocol.groups(ig); nv=group.n_videos;
    folder=fullfile(outdir,group.name); files=dir(fullfile(folder,'video_*.mat'));
    assert(numel(files)==nv,'视频文件数不符合冻结协议：%s。',group.name);
    items_all=cell(1,nv); reference_p=[];
    for iv=1:nv
        file=fullfile(folder,sprintf('video_%03d.mat',iv));
        s=load(file,'items','p','video_index','seed','protocol_hash');
        assert(s.video_index==iv && s.seed==group.seed_base+protocol.seed_stride*(iv-1));
        assert(strcmp(s.protocol_hash,protocol_hash),'协议hash不一致：%s。',file);
        assert(isequal(size(s.items),[nc nm]),'条件/方法矩阵尺寸不符：%s。',file);
        if iv==1, reference_p=s.p; else
            assert(isequaln(s.p.ana,reference_p.ana) && isequaln(s.p.phys,reference_p.phys) && ...
                s.p.sim.fps==reference_p.sim.fps,'组内推断参数不一致。');
        end
        items_all{iv}=s.items;
    end
    p=reference_p; assert(strcmpi(p.ana.fit_model,'perp'));
    Dtrue=p.phys.kB*p.phys.T/(3*pi*p.phys.eta*group.diameter_nm*1e-9)*1e12;
    rs=RandStream('mt19937ar','Seed',protocol.bootstrap_seed_base+ig);
    indices=randi(rs,nv,[nv B]);
    draws=struct('D',nan(B,nc,nm),'d',nan(B,nc,nm),'u',nan(B,nc,nm), ...
        'v',nan(B,nc,nm),'gain',nan(B,nc,nm),'left',nan(B,nc,nm),'right',nan(B,nc,nm));
    group_rows=repmat(empty_row(),nc,nm);
    for ic=1:nc
        for im=1:nm
            stats=repmat(items_all{1}(ic,im).stats,1,nv);
            metric=repmat(items_all{1}(ic,im).metric,1,nv);
            spatial_accept=false(1,nv); spatial_known=false(1,nv);
            for iv=1:nv
                stats(iv)=items_all{iv}(ic,im).stats;
                metric(iv)=items_all{iv}(ic,im).metric;
                check_stats(stats(iv),p.ana.n_lag);
                [spatial_accept(iv),spatial_known(iv)]=spatial_decision(items_all{iv}(ic,im).illumination_info);
                vr=struct('group',group.name,'condition',conditions{ic},'method',methods{im}, ...
                    'video_index',iv,'seed',group.seed_base+protocol.seed_stride*(iv-1), ...
                    'active',any(stats(iv).counts>0),'n_tracks',stats(iv).n_tracks, ...
                    'n_loc',stats(iv).n_loc,'gain_rms_pct',100*metric(iv).gain_rms, ...
                    'gain_rms_left_pct',100*metric(iv).gain_rms_left, ...
                    'gain_rms_right_pct',100*metric(iv).gain_rms_right, ...
                    'baseline_gain_rms_pct',100*metric(iv).baseline_gain_rms, ...
                    'applied',logical(metric(iv).applied),'spatial_accepted',spatial_accept(iv), ...
                    'spatial_decision_available',spatial_known(iv), ...
                    'saturated_pixels',metric(iv).saturated_pixels);
                if isempty(video_rows), video_rows=vr; else, video_rows(end+1)=vr; end %#ok<AGROW>
            end
            total=combine(stats,1:nv); est=targeted_estimate_stats(total,p);
            active=sum(arrayfun(@(x)any(x.counts>0),stats));
            row=empty_row(); row.group=group.name; row.condition=conditions{ic}; row.method=methods{im};
            row.diameter_true_nm=group.diameter_nm; row.D_true_um2_s=Dtrue;
            row.n_videos=nv; row.active_videos=active; row.n_tracks=total.n_tracks; row.n_loc=total.n_loc;
            row.D_um2_s=est.D_um2_s; row.diameter_nm=est.diameter_nm;
            row.diameter_bias_pp=100*(est.diameter_nm/group.diameter_nm-1);
            row.absolute_diameter_bias_pp=abs(row.diameter_bias_pp);
            row.D_bias_pp=100*(est.D_um2_s/Dtrue-1);
            row.flow_um_s=est.flow_um_s; row.magnetic_velocity_um_s=est.magnetic_velocity_um_s;
            row.waveform_identifiable=est.waveform_identifiable;
            row.truth_velocity_um_s=mean([metric.truth_velocity_um_s]);
            row.unmatched_selected=total.unmatched_selected;
            row.matched_localizations=total.matched_localizations;
            if total.matched_localizations>0
                row.localization_rms_nm=sqrt(total.localization_squared_sum_nm2/total.matched_localizations);
            end
            row.gain_rms_pct=pooled_rms([metric.gain_rms])*100;
            row.gain_rms_left_pct=pooled_rms([metric.gain_rms_left])*100;
            row.gain_rms_right_pct=pooled_rms([metric.gain_rms_right])*100;
            row.baseline_gain_rms_pct=pooled_rms([metric.baseline_gain_rms])*100;
            if isfield(metric,'baseline_gain_rms_left')
                row.baseline_gain_rms_left_pct=pooled_rms([metric.baseline_gain_rms_left])*100;
                row.baseline_gain_rms_right_pct=pooled_rms([metric.baseline_gain_rms_right])*100;
            end
            if row.baseline_gain_rms_pct>1e-12
                row.gain_rms_reduction_pct=100*(1-row.gain_rms_pct/row.baseline_gain_rms_pct);
            end
            row.applied_videos=sum([metric.applied]); row.saturated_pixels=sum([metric.saturated_pixels]);
            row.spatial_decision_videos=sum(spatial_known);
            row.spatial_accepted_videos=sum(spatial_accept==1);
            row.adequate_evidence=active>=10 && total.n_loc>=5000;
            for ib=1:B
                ix=indices(:,ib); eb=targeted_estimate_stats(combine(stats,ix),p);
                draws.D(ib,ic,im)=eb.D_um2_s; draws.d(ib,ic,im)=eb.diameter_nm;
                draws.u(ib,ic,im)=eb.flow_um_s; draws.v(ib,ic,im)=eb.magnetic_velocity_um_s;
            end
            x=[metric.gain_rms]; draws.gain(:,ic,im)=sqrt(mean(x(indices).^2,1))'*100;
            x=[metric.gain_rms_left]; draws.left(:,ic,im)=sqrt(mean(x(indices).^2,1))'*100;
            x=[metric.gain_rms_right]; draws.right(:,ic,im)=sqrt(mean(x(indices).^2,1))'*100;
            good=isfinite(draws.d(:,ic,im)) & isfinite(draws.D(:,ic,im)) & draws.D(:,ic,im)>0;
            row.bootstrap_valid=sum(good); row.bootstrap_replicates=B;
            velocity_good=isfinite(draws.u(:,ic,im)) & isfinite(draws.v(:,ic,im));
            row.velocity_bootstrap_valid=sum(velocity_good);
            if ~est.waveform_identifiable, row.velocity_status='unidentifiable';
            elseif active<4, row.velocity_status='descriptive_only_fewer_than_4_active_videos';
            elseif sum(velocity_good)<0.95*B, row.velocity_status='insufficient_valid_bootstrap';
            else
                q=percentile(draws.u(velocity_good,ic,im),[0.025 0.975]);
                row.flow_ci95_low=q(1); row.flow_ci95_high=q(2);
                q=percentile(draws.v(velocity_good,ic,im),[0.025 0.975]);
                row.magnetic_velocity_ci95_low=q(1); row.magnetic_velocity_ci95_high=q(2);
                row.velocity_status='exploratory_auxiliary_interval_available';
            end
            if active>=4 && row.bootstrap_valid>=0.95*B
                q=percentile(draws.d(good,ic,im),[0.025 0.975]); row.d_ci95_low=q(1); row.d_ci95_high=q(2);
                q=percentile(draws.D(good,ic,im),[0.025 0.975]); row.D_ci95_low=q(1); row.D_ci95_high=q(2);
                if row.adequate_evidence, row.status='exploratory_interval_available';
                else, row.status='interval_available_low_information'; end
            elseif active<4, row.status='descriptive_only_fewer_than_4_active_videos';
            else, row.status='insufficient_valid_bootstrap'; end
            group_rows(ic,im)=row;
            if isempty(rows), rows=row; else, rows(end+1)=row; end %#ok<AGROW>
        end
        for im=setdiff(1:nm,candidate_index,'stable')
            a=group_rows(ic,candidate_index); b=group_rows(ic,im);
            pair=compare_pair(a,b,draws,ic,candidate_index,im,group.diameter_nm,protocol.ni_margin_pp,B);
            if isempty(comparisons), comparisons=pair; else, comparisons(end+1)=pair; end %#ok<AGROW>
        end
    end
    bootstrap{ig}=struct('group',group.name,'seed',protocol.bootstrap_seed_base+ig, ...
        'indices',indices,'draws',draws);
    fprintf('candidate summary: %s complete, %d videos retained.\n',group.name,nv);
end
report=struct('schema_version',2,'protocol_hash',protocol_hash,'protocol',protocol, ...
    'data_complete',true,'inference_scope','exploratory_pointwise_no_multiplicity_adjustment', ...
    'notes',{{'独立视频为重抽样单位；空视频保留；所有条件与方法共享配对索引。', ...
    '点估计与重抽样先合并位移矩再估计D并换算d，不平均单视频粒径。', ...
    'D和d区间须至少4个有效视频、95%有效重抽样；adequate_evidence另须10个有效视频和5000个定位。', ...
    '非劣仅为预设+1百分点评价界限下的探索性结果；证据不足不能解释为等效。', ...
    '主要比较为gaussian_patch和tilted_smooth的spatial_validated对full，三个粒径分别报告；其余探索/边界；所有区间均为点态，未校正多重比较。', ...
    '照明误差为等长视频RMS平方平均再开方；静态基线为零时不计算相对降低率。', ...
    '真值匹配与定位误差只作事后诊断；数据完整不等于效果合格或专利成立。'}}, ...
    'summary',rows,'pairwise',comparisons);
writetable(struct2table(rows),fullfile(outdir,'summary.csv'));
writetable(struct2table(comparisons),fullfile(outdir,'pairwise.csv'));
writetable(struct2table(video_rows),fullfile(outdir,'video_metrics.csv'));
save(fullfile(outdir,'summary.mat'),'report','bootstrap','video_rows','-v7.3');
fid=fopen(fullfile(outdir,'summary.json'),'w','n','UTF-8'); assert(fid>=0);
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true));
end

function total=combine(stats,indices)
% 显式求和包括空视频；重复索引正是有放回视频bootstrap。
selected=stats(indices); total=selected(1); names=fieldnames(total);
for k=1:numel(names)
    name=names{k}; first=selected(1).(name); dim=ndims(first)+1;
    total.(name)=sum(cat(dim,selected.(name)),dim);
end
end

function check_stats(s,nl)
assert(isequal(size(s.counts),[nl 1]) && isequal(size(s.moments),[nl 6]));
assert(isequal(size(s.driftA),[2 2]) && isequal(size(s.driftB),[2 1]));
assert(all(isfinite(s.driftA(:))) && all(isfinite(s.driftB(:))));
assert(all(isfinite(s.counts)) && all(s.counts>=0) && all(isfinite(s.moments(:))));
assert(s.n_tracks>=0 && s.n_loc>=0 && s.unmatched_selected>=0 && ...
    s.matched_localizations>=0 && s.localization_squared_sum_nm2>=0);
end

function r=empty_row()
r=struct('group','','condition','','method','','diameter_true_nm',nan,'D_true_um2_s',nan, ...
    'n_videos',0,'active_videos',0,'n_tracks',0,'n_loc',0,'D_um2_s',nan,'diameter_nm',nan, ...
    'diameter_bias_pp',nan,'absolute_diameter_bias_pp',nan,'D_bias_pp',nan, ...
    'D_ci95_low',nan,'D_ci95_high',nan,'d_ci95_low',nan,'d_ci95_high',nan, ...
    'flow_um_s',nan,'magnetic_velocity_um_s',nan,'truth_velocity_um_s',nan, ...
    'flow_ci95_low',nan,'flow_ci95_high',nan,'magnetic_velocity_ci95_low',nan, ...
    'magnetic_velocity_ci95_high',nan,'waveform_identifiable',false, ...
    'velocity_bootstrap_valid',0,'velocity_status','', ...
    'unmatched_selected',0,'matched_localizations',0,'localization_rms_nm',nan, ...
    'gain_rms_pct',nan,'gain_rms_left_pct',nan,'gain_rms_right_pct',nan, ...
    'baseline_gain_rms_pct',nan,'baseline_gain_rms_left_pct',nan, ...
    'baseline_gain_rms_right_pct',nan,'gain_rms_reduction_pct',nan, ...
    'applied_videos',0,'spatial_decision_videos',0,'spatial_accepted_videos',0, ...
    'saturated_pixels',0,'bootstrap_valid',0,'bootstrap_replicates',0, ...
    'adequate_evidence',false,'status','');
end

function r=compare_pair(a,b,draws,ic,ia,ib,dtrue,margin,B)
r=struct('group',a.group,'condition',a.condition,'method',a.method,'reference',b.method, ...
    'primary_comparison',ismember(a.condition,{'gaussian_patch','tilted_smooth'}) && strcmp(b.method,'full'), ...
    'n_videos',a.n_videos,'active_candidate',a.active_videos,'active_reference',b.active_videos, ...
    'adequate_evidence',a.adequate_evidence && b.adequate_evidence, ...
    'paired_bootstrap_valid',0,'bootstrap_replicates',B,'ni_margin_pp',margin, ...
    'delta_D_um2_s',a.D_um2_s-b.D_um2_s,'delta_D_ci95_low',nan,'delta_D_ci95_high',nan, ...
    'delta_D_upper95',nan,'delta_d_pp',100*(a.diameter_nm-b.diameter_nm)/dtrue, ...
    'delta_d_ci95_low',nan,'delta_d_ci95_high',nan,'delta_d_upper95',nan, ...
    'delta_abs_bias_pp',a.absolute_diameter_bias_pp-b.absolute_diameter_bias_pp, ...
    'delta_abs_bias_ci95_low',nan,'delta_abs_bias_ci95_high',nan, ...
    'delta_abs_bias_lower95',nan,'delta_abs_bias_upper95',nan, ...
    'delta_gain_rms_pp',a.gain_rms_pct-b.gain_rms_pct, ...
    'delta_gain_ci95_low',nan,'delta_gain_ci95_high',nan,'delta_gain_upper95',nan, ...
    'delta_gain_left_pp',a.gain_rms_left_pct-b.gain_rms_left_pct, ...
    'delta_gain_left_ci95_low',nan,'delta_gain_left_ci95_high',nan, ...
    'delta_gain_right_pp',a.gain_rms_right_pct-b.gain_rms_right_pct, ...
    'delta_gain_right_ci95_low',nan,'delta_gain_right_ci95_high',nan, ...
    'interval_status','','ni_status','insufficient_evidence','bias_improvement_status','insufficient_evidence');
good=isfinite(draws.d(:,ic,ia)) & isfinite(draws.d(:,ic,ib)) & ...
    isfinite(draws.D(:,ic,ia)) & isfinite(draws.D(:,ic,ib)) & ...
    draws.D(:,ic,ia)>0 & draws.D(:,ic,ib)>0;
r.paired_bootstrap_valid=sum(good);
if min(a.active_videos,b.active_videos)<4
    r.interval_status='descriptive_only_fewer_than_4_active_videos';
elseif sum(good)<0.95*B
    r.interval_status='insufficient_valid_paired_bootstrap';
else
    r.interval_status='exploratory_pointwise_interval_available';
    x=draws.D(good,ic,ia)-draws.D(good,ic,ib); q=percentile(x,[0.025 0.975 0.95]);
    r.delta_D_ci95_low=q(1); r.delta_D_ci95_high=q(2); r.delta_D_upper95=q(3);
    x=100*(draws.d(good,ic,ia)-draws.d(good,ic,ib))/dtrue; q=percentile(x,[0.025 0.975 0.95]);
    r.delta_d_ci95_low=q(1); r.delta_d_ci95_high=q(2); r.delta_d_upper95=q(3);
    x=100*(abs(draws.d(good,ic,ia)/dtrue-1)-abs(draws.d(good,ic,ib)/dtrue-1));
    q=percentile(x,[0.025 0.975 0.05 0.95]);
    r.delta_abs_bias_ci95_low=q(1); r.delta_abs_bias_ci95_high=q(2);
    r.delta_abs_bias_lower95=q(3); r.delta_abs_bias_upper95=q(4);
    if r.adequate_evidence
        if q(4)<margin, r.ni_status='supports_exploratory_noninferiority';
        elseif q(3)>margin, r.ni_status='exceeds_prespecified_margin';
        else, r.ni_status='inconclusive'; end
        if q(2)<0, r.bias_improvement_status='supports_exploratory_improvement';
        elseif q(1)>0, r.bias_improvement_status='supports_exploratory_worsening';
        else, r.bias_improvement_status='inconclusive'; end
    end
end
% 照明量对所有视频均有定义，不因未检测到颗粒而删除；至少4个视频才给区间。
if a.n_videos>=4
    q=finite_interval(draws.gain(:,ic,ia)-draws.gain(:,ic,ib),B);
    r.delta_gain_ci95_low=q(1); r.delta_gain_ci95_high=q(2); r.delta_gain_upper95=q(3);
    q=finite_interval(draws.left(:,ic,ia)-draws.left(:,ic,ib),B);
    r.delta_gain_left_ci95_low=q(1); r.delta_gain_left_ci95_high=q(2);
    q=finite_interval(draws.right(:,ic,ia)-draws.right(:,ic,ib),B);
    r.delta_gain_right_ci95_low=q(1); r.delta_gain_right_ci95_high=q(2);
end
end

function q=finite_interval(x,B)
good=isfinite(x); q=nan(1,3);
if sum(good)>=0.95*B, q=percentile(x(good),[0.025 0.975 0.95]); end
end

function y=pooled_rms(x)
y=sqrt(mean(x.^2));
end

function [decision,available]=spatial_decision(info)
% is_spatial才表示空间模型接受；applied还可能仅表示回退的全局校正。
available=isfield(info,'is_spatial'); decision=false;
if ~available, return; end
value=info.is_spatial;
assert(isscalar(value) && (islogical(value) || isnumeric(value)) && ...
    isfinite(value) && ismember(value,[0 1]),'空间门控标记须为有限逻辑标量。');
decision=logical(value);
end

function q=percentile(x,probs)
% 与原群体分析相同的线性插值经验分位数，不需要Statistics Toolbox。
x=sort(x(:)); pos=1+(numel(x)-1)*probs;
lo=floor(pos); hi=ceil(pos); q=x(lo)'.*(hi-pos)+x(hi)'.*(pos-lo);
same=lo==hi; q(same)=x(lo(same));
end
