function out = analyze_tracks(tracks, traj, p, video_id)
%ANALYZE_TRACKS 检测轨迹→扩散；真实轨迹只用于最后的独立误差评价。
if nargin<4, video_id=1; end
out=empty_tracks();
N=p.cam.Npx; pxu=p.cam.px*1e6;
fitp=p;
% 波形去漂移在样本组层面拟合；单轨迹只显示垂直方向结果。
if strcmpi(p.ana.fit_model,'waveform'), fitp.ana.fit_model='perp'; end
for k=1:numel(tracks)
    tr=tracks(k); x=(tr.x-(N+1)/2)*pxu; y=(tr.y-(N+1)/2)*pxu;
    m=track_moments([x,y],tr.frames,p);
    [D,b,vfit]=fit_msd(m.tau,m.msd,m.msd_perp,fitp);
    dmes=nan;
    if isfinite(D) && D>0, dmes=p.phys.kB*p.phys.T/(3*pi*p.phys.eta*D*1e-12)*1e9; end
    Cmax=max(tr.C); cbrtC=sign(Cmax)*abs(Cmax)^(1/3);
    vdrift=nan;
    if m.drift_A(2,2)>0, vdrift=m.drift_b(2)/m.drift_A(2,2); end
    % 以下字段仅用于验证、图中参考线和误差报告；不参与D计算或入选规则。
    dtrue=nan; Dtrue=nan; vtrue=nan; sizeidx=0; md=nan; trueid=0; locrms=nan;
    if ~isempty(traj) && isfield(traj,'pos_wrap')
        if isfield(traj,'pos_camera_wrap'), truth=traj.pos_camera_wrap; else, truth=traj.pos_wrap; end
        nf=numel(tr.frames); Pt=reshape(truth(tr.frames,:,1:2)*1e6,[nf,numel(traj.d),2]);
        dist=sqrt((Pt(:,:,1)-x).^2+(Pt(:,:,2)-y).^2);
        [md,j]=min(mean(dist,1));
        if md<0.5
            trueid=j; dtrue=traj.dH(j)*1e9; sizeidx=traj.size_idx(j);
            Dtrue=traj.D(j)*1e12; vtrue=traj.v_mag(j)*1e6*norm(traj.e_g(1:2));
            locrms=sqrt(mean(dist(:,j).^2))*1e3;
        end
    end
    out(k).video=video_id; out(k).frames=tr.frames;
    out(k).x=tr.x; out(k).y=tr.y; out(k).C=tr.C; out(k).n_loc=tr.n_loc;
    out(k).x_um=x; out(k).y_um=y;
    out(k).tau=m.tau; out(k).msd=m.msd; out(k).msd_xy=m.msd_dim;
    out(k).msd_perp=m.msd_perp; out(k).npairs=m.npairs; out(k).npairs_perp=m.npairs;
    out(k).wave_q2=m.wave_q2; out(k).wave_qdt=m.wave_qdt;
    out(k).wave_cross=m.wave_cross; out(k).flow_cross=m.flow_cross;
    out(k).drift_A=m.drift_A; out(k).drift_b=m.drift_b; out(k).drift_n=m.drift_n;
    out(k).D=D; out(k).b=b; out(k).v_fit=vfit; out(k).d_mes=dmes;
    out(k).C_max=Cmax; out(k).cbrtC=cbrtC; out(k).v_drift=vdrift;
    out(k).d_true=dtrue; out(k).size_idx=sizeidx; out(k).D_true=Dtrue;
    out(k).v_true=vtrue; out(k).match_dist=md; out(k).truth_id=trueid; out(k).loc_rms_nm=locrms;
end
end
