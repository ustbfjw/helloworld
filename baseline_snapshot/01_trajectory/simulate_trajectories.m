function traj = simulate_trajectories(p, seed)
%SIMULATE_TRAJECTORIES  纯水中球形磁颗粒的布朗运动、交变磁泳及曝光采样。
%   过阻尼模型：dr=(F/gamma)dt+sqrt(2Ddt)N(0,1)，D=kBT/(3*pi*eta*dH)。
%   恒偏置场下的磁矩固定为m(B0)，磁位移使用精确波形积分。
%   x/y为开域，L_xy只指定有限颗粒的初始区域，不循环包装或重新注入。
%   球心z在[dH/2,H_z-dH/2]反射。采用体相D，不含壁面阻滞。
%   traj.t=k/fps为曝光中心；pos是中心瞬时位置（未包装）。
%   pos_exposure是曝光内均匀中点位置，pos_camera是其平均值，供定位诊断。
%   相机必须平均各样本的光强，而非平均复场或仅渲染pos_camera。
if nargin < 2, seed = p.sim.seed; end
set_seed(seed);
if isfield(p.sim,'boundary_xy')
    assert(strcmpi(p.sim.boundary_xy,'open'),'当前完整PSF模型只支持open横向边界。');
end
d = []; size_idx = [];
for k = 1:numel(p.part.d_list)
    d = [d;repmat(p.part.d_list(k),p.part.n_per_size(k),1)]; %#ok<AGROW>
    size_idx = [size_idx;repmat(k,p.part.n_per_size(k),1)]; %#ok<AGROW>
end
Np = numel(d); dH = d+2*p.part.t_shell;
assert(Np>0 && all(dH>0) && all(dH<p.sim.H_z),'粒子直径或仿真盒高度无效。');
gamma = 3*pi*p.phys.eta*dH; D = p.phys.kB*p.phys.T./gamma;
m = magnetic_moment(d,p.mag.B0,p); v_mag = m*p.mag.G0./gamma;
if p.sim.gravity
    v_grav = -(p.part.rho_p-p.phys.rho_f)*(pi/6*d.^3)*p.phys.g./gamma;
else
    v_grav = zeros(Np,1);
end
e_g = p.mag.dir(:).'; assert(norm(e_g)>0,'磁梯度方向不能为零。');
e_g = e_g/norm(e_g);
F = p.sim.n_frames; dt = 1/p.sim.fps; te = 0;
if isfield(p.cam,'exposure'), te = p.cam.exposure; end
assert(te>=0 && te<=dt,'曝光时间必须在[0,1/fps]范围内。');
nq = p.sim.n_sub;
assert(nq>=1 && nq==round(nq),'n_sub必须是正整数曝光采样数。');
if te==0, nq=1; end
offsets = (((1:nq)-0.5)/nq-0.5)*te;
% 将曝光中心额外加入采样时间：即使nq为偶数也保存准确中心真位置。
[nodes,~,node_id] = unique([offsets,0]);
tframes = (1:F)'*dt;
zlo = dH/2; zhi = p.sim.H_z-dH/2;
if isfield(p.sim,'z_init') && ~isempty(p.sim.z_init)
    zlo_init = max(zlo,p.sim.z_init(1));
    zhi_init = min(zhi,p.sim.z_init(2));
else
    zlo_init=zlo; zhi_init=zhi;
end
assert(all(zhi_init>=zlo_init),'z_init与球心可达区间没有交集。');
L=p.sim.L_xy;
r=[L*(rand(Np,1)-0.5),L*(rand(Np,1)-0.5),zlo_init+(zhi_init-zlo_init).*rand(Np,1)];
pos=zeros(F,Np,3); pos_exposure=zeros(F,Np,3,nq);
previous_t=0;
for f=1:F
    for k=1:numel(nodes)
        next_t=tframes(f)+nodes(k); step=next_t-previous_t;
        q=field_integral(previous_t,next_t,p);
        r=r+(v_mag*q)*e_g+[zeros(Np,2),v_grav*step] ...
            +sqrt(2*D*step).*randn(Np,3);
        % 三角波折返支持一次跨越多个壁面，球心永远留在有效区间。
        width=zhi-zlo; z=mod(r(:,3)-zlo,2*width);
        r(:,3)=zlo+width-abs(z-width);
        if k==node_id(end), pos(f,:,:)=reshape(r,[1,Np,3]); end
        sample_ids=find(node_id(1:nq)==k);
        for iq=sample_ids(:).'
            pos_exposure(f,:,:,iq)=reshape(r,[1,Np,3]);
        end
        previous_t=next_t;
    end
end
pos_camera=mean(pos_exposure,4);
% 保留旧接口字段名；开域下包装字段与原坐标相同。
pos_wrap=pos; pos_camera_wrap=pos_camera;
traj.t=tframes; traj.dt=dt; traj.pos=pos; traj.pos_wrap=pos_wrap;
traj.pos_exposure=pos_exposure; traj.t_exposure=tframes+offsets;
traj.pos_camera=pos_camera; traj.pos_camera_wrap=pos_camera_wrap;
traj.exposure=te; traj.exposure_samples=nq;
traj.d=d; traj.dH=dH; traj.D=D; traj.gamma=gamma; traj.m=m;
traj.v_mag=v_mag; traj.v_grav=v_grav; traj.size_idx=size_idx;
traj.s_step=field_integral(tframes-dt,tframes,p)/dt;
traj.e_g=e_g; traj.n_frames=F;
end
