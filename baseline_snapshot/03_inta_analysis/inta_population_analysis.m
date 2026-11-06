function res = inta_population_analysis(tracks, p)
%INTA_POPULATION_ANALYSIS 单分散样本组联合估计，不用真值标签筛选或分组。
% 主结果先汇总位移矩估计D，再做Stokes–Einstein换算，避免个体倒数平均。
% 主方法perp不受固定方向磁泳影响；waveform利用已知场波形作独立对照。
% 置信区间以完整独立视频为抽样单位，不把同一粒子的碎片视为独立样本。
nl=p.ana.n_lag; nv=p.sim.n_videos;
if ~isempty(tracks), nv=max(nv,max([tracks.video])); end
res.sel=false(1,numel(tracks)); res.X=zeros(0,2); res.w=zeros(0,1);
res.gm=[]; res.cluster_table=[];
minpairs=option(p.ana,'min_pairs',5);
for k=1:numel(tracks)
    np=tracks(k).npairs; nf=min(p.ana.n_fit,numel(np));
    res.sel(k)=tracks(k).n_loc>=p.ana.min_len && ...
        sum(np(1:nf)>=minpairs & isfinite(tracks(k).msd_perp(1:nf)))>=2;
end

% 每个视频保存可相加的充分统计量；各滞后用真实有效位移对数作权重。
counts=zeros(nl,nv); moments=zeros(nl,6,nv);
driftA=zeros(2,2,nv); driftB=zeros(2,nv);
for k=find(res.sel)
    t=tracks(k); iv=t.video;
    M=[t.msd,t.msd_perp,t.wave_q2,t.wave_qdt,t.wave_cross,t.flow_cross];
    ok=all(isfinite(M),2) & t.npairs>0;
    counts(ok,iv)=counts(ok,iv)+t.npairs(ok);
    moments(ok,:,iv)=moments(ok,:,iv)+M(ok,:).*t.npairs(ok);
    driftA(:,:,iv)=driftA(:,:,iv)+t.drift_A;
    driftB(:,iv)=driftB(:,iv)+t.drift_b;
end
g=estimate(sum(counts,2),sum(moments,3),sum(driftA,3),sum(driftB,2),p);
g.n_tracks=sum(res.sel); g.n_loc=sum([tracks(res.sel).n_loc]);
g.n_videos=nv; g.n_videos_with_tracks=sum(any(counts>0,1));
g.n_pairs=sum(counts,2); g.fit_model=p.ana.fit_model;
g.d_ci95_interval=[nan nan]; g.D_ci95_interval=[nan nan];
g.d_ci95=nan; g.bootstrap_valid=0;
nb=option(p.ana,'bootstrap_n',500);
if nv>=4 && g.n_videos_with_tracks>=4 && nb>0
    % 独立随机流保证调用统计函数不会改变后续仿真/成像的随机序列。
    rs=RandStream('mt19937ar','Seed',option(p.ana,'bootstrap_seed',31415));
    db=nan(nb,1); Db=nan(nb,1);
    for ib=1:nb
        ix=randi(rs,nv,[1,nv]);
        gb=estimate(sum(counts(:,ix),2),sum(moments(:,:,ix),3), ...
            sum(driftA(:,:,ix),3),sum(driftB(:,ix),2),p);
        db(ib)=gb.d_avgmsd; Db(ib)=gb.D_avgmsd;
    end
    good=isfinite(db) & isfinite(Db) & Db>0;
    g.bootstrap_valid=sum(good);
    % 大量重抽样不能估计时，不剔除坏样本后给出虚假的窄置信区间。
    if sum(good)>=0.95*nb
        g.d_ci95_interval=percentile(db(good),[0.025 0.975]);
        g.D_ci95_interval=percentile(Db(good),[0.025 0.975]);
        g.d_ci95=diff(g.d_ci95_interval)/2;
    end
end

% 个体粒径仅作分布描述；非正D不换算成直径，但其MSD仍进入组估计。
g.d_mean=nan; g.d_std=nan; g.D_mean=nan; g.n_positive_individual=0;
g.cbrtC_mean=nan; g.cbrtC_std=nan; g.v_std=nan;
if any(res.sel)
    tt=tracks(res.sel); w=[tt.n_loc]; w=w/sum(w);
    good=isfinite([tt.d_mes]); g.n_positive_individual=sum(good);
    if any(good)
        ww=w(good)/sum(w(good)); dd=[tt(good).d_mes];
        g.d_mean=sum(ww.*dd); g.d_std=sqrt(sum(ww.*(dd-g.d_mean).^2));
    end
    dd=[tt.D]; good=isfinite(dd);
    if any(good), g.D_mean=sum(w(good).*dd(good))/sum(w(good)); end
    cc=[tt.cbrtC]; good=isfinite(cc);
    if any(good)
        ww=w(good)/sum(w(good)); g.cbrtC_mean=sum(ww.*cc(good));
        g.cbrtC_std=sqrt(sum(ww.*(cc(good)-g.cbrtC_mean).^2));
    end
    good=isfinite([tt.d_mes]) & isfinite([tt.cbrtC]);
    res.X=[[tt(good).d_mes]',[tt(good).cbrtC]']; res.w=[tt(good).n_loc]';
end

% 验证参考值在所有推断完成后才写入；不改变sel、D或模型选择。
active=find(p.part.n_per_size>0);
g.d_true=nan; g.D_true=nan; g.v_true=nan; g.rel_err=nan;
if numel(active)==1
    g.d_true=(p.part.d_list(active)+2*p.part.t_shell)*1e9;
    g.D_true=p.phys.kB*p.phys.T/(3*pi*p.phys.eta*g.d_true*1e-9)*1e12;
    g.rel_err=g.d_avgmsd/g.d_true-1;
end
if any(res.sel)
    vt=[tracks(res.sel).v_true]; vt=vt(isfinite(vt));
    if ~isempty(vt), g.v_true=mean(vt); end
end
res.group=g;
% 旧绘图接口兼容；组内估计没有按真实粒子标签分流。
res.by_size=repmat(g,1,numel(p.part.d_list));
for j=1:numel(p.part.d_list)
    if ~ismember(j,active)
        res.by_size(j).n_tracks=0; res.by_size(j).n_loc=0;
        res.by_size(j).d_true=(p.part.d_list(j)+2*p.part.t_shell)*1e9;
    end
end
end

function g = estimate(counts,sums,A,b,p)
% 固定坐标系的位移矩，直接生成三种可比较估计。
M=sums./max(counts,1); M(counts==0,:)=nan;
g.tau=(1:p.ana.n_lag)'/p.sim.fps;
g.msd_avg=M(:,1); g.msd_perp_avg=M(:,2);
g.flow_mean=nan; g.v_mean=nan; g.waveform_identifiable=false;
if A(1,1)>0
    if A(2,2)==0
        g.flow_mean=b(1)/A(1,1); g.v_mean=0; g.waveform_identifiable=true;
    elseif rcond(A)>1e-8
        coef=A\b; g.flow_mean=coef(1); g.v_mean=coef(2); g.waveform_identifiable=true;
    end
end
u=g.flow_mean; v=g.v_mean;
g.msd_corrected_avg=M(:,1)-2*v*M(:,5)-2*u*M(:,6)+ ...
    v^2*M(:,3)+2*u*v*M(:,4)+u^2*g.tau.^2;
g.model_estimates=struct();
for name={'linear','perp','waveform','drift'}
    model=name{1}; pp=p;
    if strcmp(model,'waveform'), pp.ana.fit_model='linear'; xy=g.msd_corrected_avg;
    else, pp.ana.fit_model=model; xy=g.msd_avg; end
    [D,intercept,vfit]=fit_msd(g.tau,xy,g.msd_perp_avg,pp);
    d=nan;
    if isfinite(D) && D>0, d=p.phys.kB*p.phys.T/(3*pi*p.phys.eta*D*1e-12)*1e9; end
    g.model_estimates.(model)=struct('D',D,'d',d,'b',intercept,'v',vfit);
end
chosen=g.model_estimates.(lower(p.ana.fit_model));
g.D_avgmsd=chosen.D; g.d_avgmsd=chosen.d; g.b_avgmsd=chosen.b;
g.D_linear=g.model_estimates.linear.D; g.D_waveform=g.model_estimates.waveform.D;
end

function x=option(s,name,default)
if isfield(s,name), x=s.(name); else, x=default; end
end

function q=percentile(x,probs)
% 线性插值经验分位数，不依赖Statistics Toolbox。
x=sort(x(:)); pos=1+(numel(x)-1)*probs;
lo=floor(pos); hi=ceil(pos); q=x(lo)'.*(hi-pos)+x(hi)'.*(pos-lo);
same=lo==hi; q(same)=x(lo(same));
end
