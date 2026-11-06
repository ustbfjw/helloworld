function [video, meta] = render_iscat_video(traj, p, tab)
%RENDER_ISCAT_VIDEO  反射式iSCAT相机：曝光内先求瞬时光强，再进行时间积分。
%   E=1+sum(eps*exp(i*phi)*H)，相机测量 mean(|E|^2)，不是|mean(E)|^2。
%   每像素背景光电子数p.cam.photons指整次曝光总量；积分后仅添加一次噪声。
%   散射场为标量Debye模型，含离焦/干涉相位；未实现界面矢量光学及像素面积积分。
%   横向开域：有限初始颗粒自然扩散离开，不重新注入。每颗粒渲染完整相机矩形，
%   无径向窗函数；查找表不足时自动扩展，R_max不决定物理成像支持。
%   数组第一维是y/行，第二维是x/列。输出uint16 Npx*Npx*F。
Npx=p.cam.Npx; px=p.cam.px; F=traj.n_frames; Np=numel(traj.d);
if isfield(p.sim,'boundary_xy')
    assert(strcmpi(p.sim.boundary_xy,'open'),'当前完整PSF模型只支持open横向边界。');
end
assert(tab.u(1)<=-p.opt.z_f && tab.u(end)>=p.sim.H_z-p.opt.z_f-1e-15, ...
    'iPSF查找表必须覆盖整个仿真盒的离焦范围。');
coords=((1:Npx)-(Npx+1)/2)*px;
[Xg,Yg]=meshgrid(coords,coords);
% 确保所有曝光位置到相机角点的距离均被查表覆盖；禁止边界clamp伪造尾场。
if isfield(traj,'pos_exposure')
    positions=traj.pos_exposure;
else
    positions=traj.pos_wrap;
end
xx=positions(:,:,1,:); yy=positions(:,:,2,:);
rho_needed=hypot(max(abs(xx(:)))+max(abs(coords)),max(abs(yy(:)))+max(abs(coords)));
extended=rho_needed>tab.rho(end);
if extended, tab=ipsf_lookup_table(p,rho_needed); end
[eps_s,phi_s,C_peak]=particle_contrast(traj.d,p.part.n_p,p);
amp=eps_s(:).*exp(1i*phi_s(:));
illum=exp(-2*(Xg.^2+Yg.^2)/p.cam.illum_w^2);
sig_px=p.cam.bg_corr/px;
kw=ceil(3*sig_px); [kx,ky]=meshgrid(-kw:kw);
gk=exp(-(kx.^2+ky.^2)/(2*sig_px^2)); gk=gk/sum(gk(:));
bg=conv2(randn(Npx),gk,'same'); bg=bg/std(bg(:));
bg_static=1+p.cam.bg_amp*bg; base=illum.*bg_static;
assert(all(base(:)>0),'静态背景出现非正强度，请减小bg_amp。');
if isfield(traj,'pos_exposure')
    nq=size(traj.pos_exposure,4);
else
    nq=1;  % 兼容外部提供的瞬时真轨迹
end
adu_max=2^p.cam.bits-1;
% 照明扰动在噪声及ADC之前、曝光积分内部注入，不能事后乘含噪图像冒充光子模型。
fa=0; ff=0; fp=0;
if isfield(p.cam,'flicker_amplitude')
    fa=p.cam.flicker_amplitude; ff=p.cam.flicker_frequency; fp=p.cam.flicker_phase;
end
assert(isnumeric(fa) && isscalar(fa) && isreal(fa) && isfinite(fa) && abs(fa)<1 && ...
    isnumeric(ff) && isscalar(ff) && isreal(ff) && isfinite(ff) && ff>=0 && ...
    isnumeric(fp) && isscalar(fp) && isreal(fp) && isfinite(fp),'周期照明参数无效。');
if isfield(traj,'t_exposure')
    illumination_gain=1+fa*sin(2*pi*ff*traj.t_exposure+fp);
else
    illumination_gain=1+fa*sin(2*pi*ff*(1:F)'/p.sim.fps+fp);
    illumination_gain=repmat(illumination_gain,1,nq);
end
video=zeros(Npx,Npx,F,'uint16'); n_sat=0;
n_rendered=zeros(Np,1); n_render_samples=zeros(Np,1);
for f=1:F
    Imean=zeros(Npx); rendered=false(Np,1);
    for iq=1:nq
        E=complex(ones(Npx),zeros(Npx));
        if isfield(traj,'pos_exposure')
            r=reshape(traj.pos_exposure(f,:,:,iq),[Np,3]);
        else
            r=reshape(traj.pos_wrap(f,:,:),[Np,3]);
        end
        for j=1:Np
            u=r(j,3)-p.opt.z_f;
            assert(u>=tab.u(1)-1e-15 && u<=tab.u(end)+1e-15,'粒子离焦量超出查找表。');
            % 插值区间端点最大值给出幅值上界；含散射强度项，不再硬删u>2.5um。
            fu=(u-tab.u(1))/tab.du+1;
            iu=min(max(floor(fu),1),numel(tab.u)-1);
            amax=abs(amp(j))*max(tab.Hmax_u(iu:iu+1));
            if 2*amax+amax^2<p.opt.C_min, continue; end
            rho=hypot(Xg-r(j,1),Yg-r(j,2));
            H=ipsf_eval(tab,rho,u);
            E=E+amp(j)*H;
            rendered(j)=true; n_render_samples(j)=n_render_samples(j)+1;
        end
        Imean=Imean+abs(E).^2*illumination_gain(f,iq)/nq;
    end
    n_rendered=n_rendered+rendered;
    electrons=p.cam.photons*Imean.*base;
    if ~isfield(p.cam,'noise') || p.cam.noise
        electrons=electrons+sqrt(electrons).*randn(Npx)+p.cam.read_e*randn(Npx);
    end
    adu=round(electrons*p.cam.gain);
    n_sat=n_sat+sum(adu(:)>adu_max);
    video(:,:,f)=uint16(min(max(adu,0),adu_max));
end
meta.eps=eps_s; meta.phi=phi_s; meta.C_peak=C_peak;
meta.illum=illum; meta.bg_static=bg_static;
meta.n_sat=n_sat; meta.n_rendered=n_rendered; meta.n_render_samples=n_render_samples;
meta.exposure_samples=nq;
meta.lookup_rho_max=tab.rho(end); meta.lookup_extended=extended;
meta.window_model='full_camera_open_domain';
meta.illumination_gain=mean(illumination_gain,2); % 仅诊断真值，不传给背景估计器
if isfield(traj,'exposure'), meta.exposure=traj.exposure; else, meta.exposure=0; end
if n_sat>0
    fprintf('  [提示] 共%d个像素饱和，可降低相机增益。\n',n_sat);
end
end
