function estimate=targeted_estimate_stats(stats,p)
%TARGETED_ESTIMATE_STATS 从相加统计量重建原perp组估计及波形速度对照。
% 此处调用原fit_msd；不平均单视频直径，不用真尺寸修正估计值。
assert(strcmpi(p.ana.fit_model,'perp'),'定向实验主估计器固定为perp。');
M=stats.moments./max(stats.counts,1); M(stats.counts==0,:)=nan;
tau=(1:p.ana.n_lag)'/p.sim.fps;
[D,b]=fit_msd(tau,M(:,1),M(:,2),p);
d=nan;
if isfinite(D) && D>0
    d=p.phys.kB*p.phys.T/(3*pi*p.phys.eta*D*1e-12)*1e9;
end
u=nan; v=nan; identifiable=false;
A=stats.driftA; rhs=stats.driftB;
if A(1,1)>0
    if A(2,2)==0
        u=rhs(1)/A(1,1); v=0; identifiable=true;
    elseif rcond(A)>1e-8
        coefficient=A\rhs; u=coefficient(1); v=coefficient(2); identifiable=true;
    end
end
estimate=struct('D_um2_s',D,'diameter_nm',d,'intercept_um2',b, ...
    'flow_um_s',u,'magnetic_velocity_um_s',v, ...
    'waveform_identifiable',identifiable);
end
