function m = track_moments(pos, frames, p)
%TRACK_MOMENTS 从检测坐标构造统一的位移对，不使用粒径、磁矩或轨迹真值。
% 所有方向及磁修正共用原始二维视场中的筛选，避免旋转投影后边界错误。
% pos单位µm；Q是曝光平均后的已知场波形原函数，单位s。
dt = 1/p.sim.fps; nl = p.ana.n_lag;
eg = reshape(p.mag.dir(1:2),1,2);
if norm(eg)>1e-12, eg=eg/norm(eg); else, eg=[1 0]; end
ep = [-eg(2),eg(1)]; frames=frames(:);
fr = frames-frames(1)+1;
P=nan(fr(end),2); P(fr,:)=pos;
t=(frames(1):frames(end))'*dt;
if isfield(p.cam,'exposure') && p.cam.exposure>0
    ns=p.sim.n_sub;
    offset=(((1:ns)-0.5)/ns-0.5)*p.cam.exposure;
else
    offset=0;
end
Q=mean(field_integral(0,t+offset,p),2);
if p.mag.G0==0 || norm(p.mag.dir(1:2))<1e-12, Q(:)=0; end
% 初步步幅仅决定边界余量；不把已知真D作为输入。
dp=diff(P); valid=all(isfinite(dp),2);
if any(valid), m1=mean(sum(dp(valid,:).^2,2)); else, m1=0; end
margin=p.ana.edge_nsig*sqrt(max(m1,0)*(1:nl)'/2);
bound=(p.cam.Npx/2-(p.ana.loc_win+2)-0.5)*p.cam.px*1e6;
m.tau=(1:nl)'*dt; m.npairs=zeros(nl,1);
m.msd=nan(nl,1); m.msd_dim=nan(nl,2); m.msd_perp=nan(nl,1);
m.wave_q2=nan(nl,1); m.wave_qdt=nan(nl,1);
m.wave_cross=nan(nl,1); m.flow_cross=nan(nl,1);
m.drift_A=zeros(2); m.drift_b=zeros(2,1); m.drift_n=0;
for lag=1:nl
    if lag>=size(P,1), break; end
    P0=P(1:end-lag,:); delta=P(1+lag:end,:)-P0;
    valid=all(isfinite(delta),2) & all(abs(P0)<bound-margin(lag),2);
    n=sum(valid); m.npairs(lag)=n;
    if n==0, continue; end
    delta=delta(valid,:); q=Q(1+lag:end)-Q(1:end-lag); q=q(valid);
    along=delta*eg'; elapsed=lag*dt;
    m.msd_dim(lag,:)=mean(delta.^2,1); m.msd(lag)=sum(m.msd_dim(lag,:));
    m.msd_perp(lag)=mean((delta*ep').^2);
    m.wave_q2(lag)=mean(q.^2); m.wave_qdt(lag)=mean(q)*elapsed;
    m.wave_cross(lag)=mean(q.*along); m.flow_cross(lag)=mean(along)*elapsed;
    if lag==1
        % 组内估计常量流动与交变漂移；不对每条短轨迹单独去趋势。
        A=[elapsed*ones(n,1),q];
        m.drift_A=A'*A/elapsed; m.drift_b=A'*along/elapsed; m.drift_n=n;
    end
end
end
