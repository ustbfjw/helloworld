function [videos, metas] = render_targeted_conditions(traj, p, tab, conditions)
%RENDER_TARGETED_CONDITIONS 定向对照专用：共享轨迹、曝光光场及标准正态噪声。
%   conditions字段：name、amplitude、frequency、phase、spatial_mode。
%   spatial_mode为global或left_half；后者仅调制物面x<0的相机左半边。
%   每次曝光先计算瞬时|E|^2，再乘照明增益并积分，最后添加噪声和ADC。
%   各条件共享zshot/zread，但散粒噪声尺度sqrt(electrons)按各自光子数计算。
%   global保留生产render_iscat_video的算序；一次调用的RNG消耗与单条件相同。
%   真增益图为1+meta.gain_mask*meta.gain_offset(f)，仅供事后评价。
%   局部照明下meta.illumination_gain仅为空间平均诊断，不能替代全像素误差。
%   此helper不修改生产渲染器，也不将照明真值传入任何估计器。
required={'name','amplitude','frequency','phase','spatial_mode'};
assert(isstruct(conditions) && ~isempty(conditions) && isvector(conditions) && ...
    all(isfield(conditions,required)), 'targeted:InvalidConditions', ...
    'conditions必须为非空结构体向量，含全部指定字段。');
K=numel(conditions);
for k=1:K
    c=conditions(k);
    valid_name=(ischar(c.name) && isrow(c.name) && ~isempty(c.name)) || ...
        (isstring(c.name) && isscalar(c.name) && ~ismissing(c.name) && strlength(c.name)>0);
    valid_mode=(ischar(c.spatial_mode) && isrow(c.spatial_mode)) || ...
        (isstring(c.spatial_mode) && isscalar(c.spatial_mode) && ~ismissing(c.spatial_mode));
    assert(valid_name && valid_mode && any(strcmp(c.spatial_mode,{'global','left_half'})) && ...
        finite_scalar(c.amplitude) && abs(c.amplitude)<1 && ...
        finite_scalar(c.frequency) && c.frequency>=0 && finite_scalar(c.phase), ...
        'targeted:InvalidConditions','条件名称、照明参数或空间模式无效。');
end

Npx=p.cam.Npx; px=p.cam.px; F=traj.n_frames; Np=numel(traj.d);
if isfield(p.sim,'boundary_xy')
    assert(strcmpi(p.sim.boundary_xy,'open'),'当前完整PSF模型只支持open横向边界。');
end
assert(tab.u(1)<=-p.opt.z_f && tab.u(end)>=p.sim.H_z-p.opt.z_f-1e-15, ...
    'iPSF查找表必须覆盖整个仿真盒的离焦范围。');
coords=((1:Npx)-(Npx+1)/2)*px;
[Xg,Yg]=meshgrid(coords,coords);
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
if isfield(traj,'pos_exposure'), nq=size(traj.pos_exposure,4); else, nq=1; end
adu_max=2^p.cam.bits-1;
videos=cell(1,K); metas=cell(1,K); illumination_gain=cell(1,K);
gain_mask=cell(1,K); is_global=false(1,K); n_sat=zeros(1,K);
for k=1:K
    fa=conditions(k).amplitude; ff=conditions(k).frequency; fp=conditions(k).phase;
    if isfield(traj,'t_exposure')
        illumination_gain{k}=1+fa*sin(2*pi*ff*traj.t_exposure+fp);
    else
        gain=1+fa*sin(2*pi*ff*(1:F)'/p.sim.fps+fp);
        illumination_gain{k}=repmat(gain,1,nq);
    end
    is_global(k)=strcmp(conditions(k).spatial_mode,'global');
    if is_global(k), gain_mask{k}=ones(Npx); else, gain_mask{k}=double(Xg<0); end
    videos{k}=zeros(Npx,Npx,F,'uint16');
end
n_rendered=zeros(Np,1); n_render_samples=zeros(Np,1);
has_noise=~isfield(p.cam,'noise') || p.cam.noise;
for f=1:F
    Imean=cell(1,K);
    for k=1:K, Imean{k}=zeros(Npx); end
    rendered=false(Np,1);
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
            fu=(u-tab.u(1))/tab.du+1;
            iu=min(max(floor(fu),1),numel(tab.u)-1);
            amax=abs(amp(j))*max(tab.Hmax_u(iu:iu+1));
            if 2*amax+amax^2<p.opt.C_min, continue; end
            rho=hypot(Xg-r(j,1),Yg-r(j,2));
            H=ipsf_eval(tab,rho,u);
            E=E+amp(j)*H;
            rendered(j)=true; n_render_samples(j)=n_render_samples(j)+1;
        end
        instantaneous_intensity=abs(E).^2;
        for k=1:K
            if is_global(k)
                % 保留原标量乘法/除法/累加顺序，避免改变量化边界附近像素。
                Imean{k}=Imean{k}+instantaneous_intensity*illumination_gain{k}(f,iq)/nq;
            else
                gain_map=ones(Npx);
                gain_map(gain_mask{k}~=0)=illumination_gain{k}(f,iq);
                Imean{k}=Imean{k}+instantaneous_intensity.*gain_map/nq;
            end
        end
    end
    n_rendered=n_rendered+rendered;
    if has_noise
        % 仅取一次随机数：同一像素同一帧的共同随机数用于配对对照。
        zshot=randn(Npx); zread=randn(Npx);
    end
    for k=1:K
        electrons=p.cam.photons*Imean{k}.*base;
        if has_noise
            electrons=electrons+sqrt(electrons).*zshot+p.cam.read_e*zread;
        end
        adu=round(electrons*p.cam.gain);
        n_sat(k)=n_sat(k)+sum(adu(:)>adu_max);
        videos{k}(:,:,f)=uint16(min(max(adu,0),adu_max));
    end
end
for k=1:K
    meta=struct;
    meta.eps=eps_s; meta.phi=phi_s; meta.C_peak=C_peak;
    meta.illum=illum; meta.bg_static=bg_static;
    meta.n_sat=n_sat(k); meta.n_rendered=n_rendered; meta.n_render_samples=n_render_samples;
    meta.exposure_samples=nq;
    meta.lookup_rho_max=tab.rho(end); meta.lookup_extended=extended;
    meta.window_model='full_camera_open_domain';
    meta.gain_mask=gain_mask{k};
    meta.gain_offset=mean(illumination_gain{k}-1,2);
    if is_global(k)
        meta.illumination_gain=mean(illumination_gain{k},2);
    else
        meta.illumination_gain=1+mean(gain_mask{k}(:))*meta.gain_offset;
    end
    if isfield(traj,'exposure'), meta.exposure=traj.exposure; else, meta.exposure=0; end
    meta.illumination_condition=conditions(k);
    meta.gain_truth_note=['曝光平均真增益图=1+gain_mask*gain_offset(f)。' ...
        '局部条件的illumination_gain仅为空间平均诊断，不能作为全像素误差真值。'];
    metas{k}=meta;
    if n_sat(k)>0
        fprintf('  [提示] 条件%s共%d个像素饱和，可降低相机增益。\n',conditions(k).name,n_sat(k));
    end
end
end

function ok=finite_scalar(value)
ok=isnumeric(value) && isscalar(value) && isreal(value) && isfinite(value);
end
