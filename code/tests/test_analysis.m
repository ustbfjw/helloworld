function test_analysis()
%TEST_ANALYSIS 磁修正、真值隔离、统一位移对和独立视频置信区间回归检查。
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root); setup_paths(root);
p=sim_params('quick'); p.ana.fit_model='perp'; p.ana.bootstrap_n=0;
p.ana.min_len=25; p.ana.min_pairs=5;
p.cam.exposure=0; p.cam.Npx=1024; p.cam.px=1e-6;
p.part.n_per_size=[0 0 1]; p.sim.n_videos=1;

% 精确解析MSD：群体perp必须使用真实垂直矩，不能用xy/2代替。
r=make_raw(zeros(50,2),p); tr=analyze_tracks(r,[],p,1);
D0=2; v0=100; t=tr.tau;
tr.msd=4*D0*t+v0^2*t.^2+0.001; tr.msd_perp=2*D0*t+0.0005;
tr.npairs=100*ones(size(t)); tr.npairs_perp=tr.npairs;
tr.n_loc=200; tr.D=-1; tr.d_mes=nan;
g=inta_population_analysis(tr,p);
assert(g.sel && abs(g.group.D_avgmsd-D0)<1e-10);
assert(abs(g.group.model_estimates.linear.D-D0-2.5)<1e-8);
assert(abs(g.group.model_estimates.drift.v-v0)<1e-6);
% 已知粒径/磁真值不能改变估计或选择，非正单轨迹D仍进入组MSD。
tr.d_true=1e8; tr.D_true=1e8; tr.v_true=1e8;
p.part.d_list=[31 111 177]*1e-9;
h=inta_population_analysis(tr,p);
assert(isequal(g.sel,h.sel) && g.group.D_avgmsd==h.group.D_avgmsd);

% 两帧之间缺一帧：链接门限应按sqrt(2)扩大，而非强制切断。
p.ana.link_r=6; p.ana.link_gap=1;
det={ [10 10 0.1 1]; zeros(0,4); [17 10 0.1 1] };
linked=link_trajectories(det,p);
assert(numel(linked)==1 && isequal(linked.frames,[1;3]));

% 丢帧与旋转方向：perp和xy共享同一二维边界筛选。
p.mag.dir=[1 1 0]/sqrt(2); p.cam.Npx=32; p.cam.px=0.25e-6;
r=make_raw([linspace(-2,2,80)',0.4*sin((1:80)')],p);
r.frames([10 11 37])=[]; r.x([10 11 37])=[]; r.y([10 11 37])=[];
r.C([10 11 37])=[]; r.n_loc=numel(r.frames);
tt=analyze_tracks(r,[],p,1);
assert(isequal(tt.npairs,tt.npairs_perp));
pcol=p; pcol.mag.dir=p.mag.dir(:);
ttcol=analyze_tracks(r,[],pcol,1);
assert(isequaln(tt.msd_perp,ttcol.msd_perp) && isequaln(tt.wave_cross,ttcol.wave_cross));
pos=[tt.x_um tt.y_um]; fr=tt.frames-tt.frames(1)+1;
P=nan(fr(end),2); P(fr,:)=pos; dp=diff(P); ok=all(isfinite(dp),2);
scale=mean(sum(dp(ok,:).^2,2)); bound=(p.cam.Npx/2-(p.ana.loc_win+2)-0.5)*p.cam.px*1e6;
ep=[-1 1]/sqrt(2);
for lag=1:p.ana.n_lag
    d=P(lag+1:end,:)-P(1:end-lag,:);
    ok=all(isfinite(d),2) & all(abs(P(1:end-lag,:))<bound-p.ana.edge_nsig*sqrt(scale*lag/2),2);
    assert(sum(ok)==tt.npairs(lag));
    if any(ok), assert(abs(mean((d(ok,:)*ep').^2)-tt.msd_perp(lag))<1e-12); end
end

% 真正的随机轨迹（无成像）：联合估计D及交变速度，不向算法传入真值。
p=sim_params('quick'); p.ana.fit_model='perp'; p.ana.bootstrap_n=200;
p.ana.bootstrap_seed=41; p.ana.min_len=25; p.ana.min_pairs=5;
p.cam.exposure=0; p.cam.Npx=1024; p.cam.px=1e-6;
p.part.n_per_size=[0 0 1]; p.sim.n_videos=30; p.mag.dir=[1 1 0]/sqrt(2);
p.mag.wave='sine'; p.mag.phase0=0.37; p.mag.f=13;
D0=2.2; v0=35; u0=3; F=4000; dt=1/p.sim.fps;
rs=RandStream('mt19937ar','Seed',73); alltracks=empty_tracks();
for iv=1:p.sim.n_videos
    times=(1:F)'*dt; q=field_integral(0,times,p);
    pos=cumsum(sqrt(2*D0*dt)*randn(rs,F,2),1)+(v0*q+u0*times)*p.mag.dir(1:2);
    pos=pos+0.01*randn(rs,F,2);
    alltracks=[alltracks,analyze_tracks(make_raw(pos,p),[],p,iv)]; %#ok<AGROW>
end
state=rng; result=inta_population_analysis(alltracks,p); state_after=rng;
assert(isequal(state,state_after),'统计函数不得污染全局随机数状态。');
g=result.group;
assert(abs(g.D_avgmsd/D0-1)<0.03);
assert(abs(g.D_waveform/D0-1)<0.03);
assert(abs(g.v_mean-v0)<1.5 && abs(g.flow_mean-u0)<1.5);
assert(all(isfinite(g.d_ci95_interval)) && diff(g.d_ci95_interval)>0);
fprintf('test_analysis PASS: Dperp=%.6f Dwaveform=%.6f v=%.6f flow=%.6f bootstrap=%d\n', ...
    g.D_avgmsd,g.D_waveform,g.v_mean,g.flow_mean,g.bootstrap_valid);

% 边界状态复用同一检测坐标，检查投影/可辨识性，不再重复生成随机样本。
pe=p; pe.ana.bootstrap_n=0; pe.sim.n_videos=1; pe.mag.dir=[1 0 0];
raw=make_raw(pos,pe); tx=analyze_tracks(raw,[],pe,1);
pe.mag.dir=[0;0;1]; tz=analyze_tracks(raw,[],pe,1);
assert(all(tz.wave_q2==0) && isequaln(tx.msd_perp,tz.msd_perp), ...
    '纯轴向梯度的面内磁项应为零，perp仍为y分量。');
pe.mag.dir=[1 0 0]; pe.mag.G0=0;
tzero=analyze_tracks(raw,[],pe,1); rzero=inta_population_analysis(tzero,pe);
assert(all(tzero.wave_q2==0) && rzero.group.v_mean==0);
assert(isequaln(tx.msd_perp,tzero.msd_perp));
empty=inta_population_analysis(empty_tracks(),pe);
assert(empty.group.n_tracks==0 && isnan(empty.group.d_avgmsd) && ...
    all(isnan(empty.group.d_ci95_interval)),'空结果不能给出有限粒径或置信区间。');
pe.mag.G0=p.mag.G0; pe.mag.f=0; pe.mag.phase0=pi/2;
tc=analyze_tracks(raw,[],pe,1); rc=inta_population_analysis(tc,pe);
assert(~rc.group.waveform_identifiable && isnan(rc.group.D_waveform) && ...
    isfinite(rc.group.D_avgmsd),'恒定波形与恒流不可区分，不应破坏perp估计。');
fprintf('test_analysis edge cases PASS: zero field, pure-z, empty data, unidentifiable constant waveform\n');
end

function r=make_raw(pos,p)
% 原始检测接口只含帧号、像素坐标及对比度。
r.frames=(1:size(pos,1))'; r.x=pos(:,1)/(p.cam.px*1e6)+(p.cam.Npx+1)/2;
r.y=pos(:,2)/(p.cam.px*1e6)+(p.cam.Npx+1)/2;
r.C=0.1*ones(size(r.x)); r.n_loc=size(pos,1);
end
