function metrics = test_numerical_convergence(outdir)
%TEST_NUMERICAL_CONVERGENCE 同一输入下检验光学离散与曝光积分收敛。
% 此测试只检验数值光学和已知中心附近定位，不代替独立盲追踪验收。
% 相机噪声关闭、背景6000 ADU，减少量化掩盖数值误差的可能。
% 阈值预先规定：LUT相对图像RMS<=5%，中心差p95<=0.25 px；
% 初始查表半径4/6um应由自动扩展产生相同图像（至多1 ADU舍入差），
% 此项检验查表覆盖，不再宣称物理窗函数收敛，因为生产模型没有径向窗。
% 曝光D的4->8变化<=3%、8->16<=2%，且8点图像比4点更接近16点。
if nargin<1, outdir=''; end
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root); setup_paths(root);
p=sim_params('full'); p.sim.L_xy=20e-6; p.sim.H_z=4e-6;
p.cam.noise=false; p.cam.bg_amp=0; p.cam.illum_w=inf;
p.cam.bits=16; p.cam.gain=6000/p.cam.photons;
p.cam.exposure=0; p.sim.n_sub=1; p.opt.R_max=4e-6;
p.opt.drho=20e-9; p.opt.du=20e-9;
tab20=ipsf_lookup_table(p);
fine=p; fine.opt.drho=10e-9; fine.opt.du=10e-9;
tab10=ipsf_lookup_table(fine);
wide=fine; wide.opt.R_max=6e-6; tab6=ipsf_lookup_table(wide);
metrics=struct(); metrics.scope='Noiseless LUT resolution, auto-extension invariance and shared-path exposure convergence; no radial truncation in production';
metrics.thresholds=struct('image_relative_rms',0.05,'localization_p95_px',0.25, ...
    'exposure_D_4_to_8',0.03,'exposure_D_8_to_16',0.02);
sizes=[30 100 200]*1e-9;
names={'mono30','mono100','mono200'}; loc_windows=zeros(1,3);
for k=1:3
    configured=apply_scenario(sim_params('full'),names{k});
    loc_windows(k)=configured.ana.loc_win;
end
metrics.localization_half_windows_px=loc_windows;
z=[0.137 0.283 0.527 0.811 1.013 1.267 1.691 2.139 2.773 3.373 3.817]*1e-6;
offset=[0.13 0.37 0.61 0.89];
for k=1:3
    p.ana.loc_win=loc_windows(k); % 与各组生产细定位窗口一致
    F=numel(z)*numel(offset); P=zeros(F,1,3); f=0;
    for iz=1:numel(z)
        for io=1:numel(offset)
            f=f+1; P(f,1,:)=[offset(io)*p.cam.px,-offset(io)*p.cam.px,z(iz)];
        end
    end
    tr=static_trace(sizes(k),P,p);
    v20=render_iscat_video(tr,p,tab20); v10=render_iscat_video(tr,fine,tab10);
    loc20=centers(v20,p,P); loc10=centers(v10,p,P);
    row=struct('diameter_nm',sizes(k)*1e9,'image_relative_rms',contrast_error(v20,v10,6000), ...
        'localization_p95_px',percentile(sqrt(sum((loc20-loc10).^2,2)),0.95));
    metrics.lookup(k)=row;
    fprintf('LUT d=%g nm relative_RMS=%.6f localization_p95_px=%.6f\n',row.diameter_nm,row.image_relative_rms,row.localization_p95_px);
    % 不同初始表尺度都必须覆盖相同完整光场；含视场外的离焦粒子。
    for mode=1:2
        tw=tr;
        if mode==2
            P2=P; P2(:,1,1)=4.6e-6; P2(:,1,2)=0.3e-6;
            P2(:,1,3)=3.5e-6+0.2e-6*sin((1:F)');
            tw=static_trace([sizes(k);sizes(k)],cat(2,P,P2),p);
        end
        v4=render_iscat_video(tw,fine,tab10); v6=render_iscat_video(tw,wide,tab6);
        l4=centers(v4,p,P); l6=centers(v6,p,P);
        row=struct('diameter_nm',sizes(k)*1e9,'particles',mode, ...
            'image_relative_rms',contrast_error(v4,v6,6000), ...
            'localization_p95_px',percentile(sqrt(sum((l4-l6).^2,2)),0.95), ...
            'max_ADU_difference',max(abs(double(v4(:))-double(v6(:)))));
        metrics.window((k-1)*2+mode)=row;
        fprintf('WINDOW d=%g nm particles=%d relative_RMS=%.6f localization_p95_px=%.6f\n', ...
            row.diameter_nm,mode,row.image_relative_rms,row.localization_p95_px);
    end
end

% 一条共同细网格轨迹包含4/8/16中点位置；禁止用同seed重跑不同n_sub作伪配对。
F=1200; p.cam.exposure=80e-6; p.sim.gravity=false;
pf=p; fine_dt=p.cam.exposure/32; ratio=round((1/p.sim.fps)/fine_dt);
assert(abs(ratio*fine_dt-1/p.sim.fps)<1e-15);
pf.sim.fps=1/fine_dt; pf.sim.n_frames=F*ratio+16;
pf.cam.exposure=0; pf.sim.n_sub=1; pf.part.n_per_size=[1 1 1];
pf.sim.z_init=[0.8e-6 1.2e-6];
allfine=simulate_trajectories(pf,1973);
% 每种粒子独立渲染，初始xy放在视场中心；相同平移应用于全部采样设置。
for k=1:3
    p.ana.loc_win=loc_windows(k);
    positions=allfine.pos(:,k,:);
    positions(:,:,1:2)=positions(:,:,1:2)-positions(1,:,1:2);
    images=cell(1,3); loc=cell(1,3); DD=zeros(1,3); nqs=[4 8 16];
    traces=cell(1,3);
    for iq=1:3
        nq=nqs(iq); offsets=(((1:nq)-0.5)/nq-0.5)*p.cam.exposure;
        index=(1:F)'*ratio+round(offsets/fine_dt);
        assert(max(abs(offsets-round(offsets/fine_dt)*fine_dt))<1e-15);
        tx=static_trace(sizes(k),positions((1:F)*ratio,:,:),p);
        tx.pos_exposure=zeros(F,1,3,nq);
        for j=1:nq, tx.pos_exposure(:,:,:,j)=positions(index(:,j),:,:); end
        tx.pos_camera=mean(tx.pos_exposure,4);
        tx.exposure=p.cam.exposure;
        images{iq}=render_iscat_video(tx,fine,tab10);
        traces{iq}=tx;
    end
    % 所有设置使用同一16点参考中心窗口，排除不同裁剪位置引入的差别。
    reference=traces{3}.pos_camera;
    for iq=1:3
        loc{iq}=centers(images{iq},p,reference);
        DD(iq)=perp_diffusion(loc{iq},p);
    end
    e4=contrast_error(images{1},images{3},6000);
    e8=contrast_error(images{2},images{3},6000);
    row=struct('diameter_nm',sizes(k)*1e9,'image_4_vs_16',e4,'image_8_vs_16',e8, ...
        'D4_um2_s',DD(1),'D8_um2_s',DD(2),'D16_um2_s',DD(3), ...
        'D_change_4_to_8',abs(DD(1)/DD(2)-1),'D_change_8_to_16',abs(DD(2)/DD(3)-1));
    metrics.exposure(k)=row;
    fprintf('EXPOSURE d=%g nm RMS4=%.6f RMS8=%.6f D4/8/16=%.6f %.6f %.6f change4to8=%.6f change8to16=%.6f\n', ...
        row.diameter_nm,e4,e8,DD,row.D_change_4_to_8,row.D_change_8_to_16);
end
metrics.lookup_pass=all([metrics.lookup.image_relative_rms]<=0.05) && ...
    all([metrics.lookup.localization_p95_px]<=0.25);
metrics.window_pass=all([metrics.window.image_relative_rms]<=0.05) && ...
    all([metrics.window.localization_p95_px]<=0.25) && all([metrics.window.max_ADU_difference]<=1);
metrics.exposure_pass=all([metrics.exposure.image_8_vs_16]<[metrics.exposure.image_4_vs_16]) && ...
    all([metrics.exposure.D_change_4_to_8]<=0.03) && all([metrics.exposure.D_change_8_to_16]<=0.02);
metrics.pass=metrics.lookup_pass && metrics.window_pass && metrics.exposure_pass;
if ~isempty(outdir)
    if ~exist(outdir,'dir'), mkdir(outdir); end
    fid=fopen(fullfile(outdir,'numerical_convergence.json'),'w','n','UTF-8');
    fprintf(fid,'%s',jsonencode(metrics,PrettyPrint=true)); fclose(fid);
    save(fullfile(outdir,'numerical_convergence.mat'),'metrics');
end
assert(metrics.pass,'数值收敛未通过；查看指标，升级默认离散精度或处理窗口误差，不应放宽阈值。');
fprintf('PASS test_numerical_convergence\n');
end

function tr=static_trace(d,P,p)
tr=struct('d',d(:),'n_frames',size(P,1),'pos_wrap',P);
end

function xy=centers(video,p,P)
% 已知中心仅用来固定裁剪窗口，是光学数值单元测试，不是生产追踪算法。
N=p.cam.Npx; L=p.ana.loc_win; F=size(video,3); xy=nan(F,2);
for f=1:F
    ci=round(P(f,1,1)/p.cam.px+(N+1)/2);
    ri=round(P(f,1,2)/p.cam.px+(N+1)/2);
    if ci<=L || ci>N-L || ri<=L || ri>N-L, continue; end
    patch=double(video(ri-L:ri+L,ci-L:ci+L,f));
    [dx,dy]=radial_center(patch);
    xy(f,:)=[ci+dx,ri+dy];
end
end

function e=contrast_error(a,b,background)
% 全帧相对对比度L2误差，避免参考场常数项掩盖粒子信号误差。
delta=double(a)-double(b); signal=double(b)-background;
e=sqrt(sum(delta(:).^2)/max(sum(signal(:).^2),eps));
end

function D=perp_diffusion(xy,p)
% 同一时间轴上的y向MSD；保留缺失帧，不把间隔压缩为连续帧。
y=xy(:,2)*p.cam.px*1e6; msd=zeros(4,1);
for lag=1:4
    dy=y(lag+1:end)-y(1:end-lag); dy=dy(isfinite(dy));
    assert(numel(dy)>=100,'共同定位帧不足，不能进行收敛判断。');
    msd(lag)=mean(dy.^2);
end
c=polyfit((1:4)'/p.sim.fps,msd,1); D=c(1)/2;
assert(isfinite(D) && D>0,'数值定位没有产生有效正扩散系数。');
end

function q=percentile(x,p)
x=sort(x(isfinite(x))); assert(~isempty(x),'没有可比较的中心定位。');
t=1+(numel(x)-1)*p; lo=floor(t); hi=ceil(t);
if lo==hi, q=x(lo); else, q=x(lo)*(hi-t)+x(hi)*(t-lo); end
end
