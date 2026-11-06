function test_physics_imaging()
%TEST_PHYSICS_IMAGING  无工具箱物理/相机回归检查；失败即抛出assert。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'00_config'),fullfile(root,'01_trajectory'), ...
    fullfile(root,'02_iscat_imaging'),fullfile(root,'utils'));
p=sim_params('quick'); p.sim.L_xy=16e-6; p.sim.H_z=4e-6;
p.cam.exposure=80e-6; p.sim.n_sub=4;
p.mag.wave='square'; p.mag.phase0=0; p.mag.f=10;
assert(abs(field_integral(0,0.05,p)-0.05)<1e-12);
assert(abs(field_integral(0.04,0.06,p))<1e-12);
assert(abs(field_integral(-0.03,0.17,p))<1e-12);
q=field_integral([0;0.05],[0.01,0.02],p);
assert(isequal(size(q),[2,2]));
p.mag.wave='sine';
assert(abs(field_integral(0,0.05,p)-1/(pi*10))<1e-12);
p.mag.f=0; p.mag.phase0=pi/2;
assert(abs(field_integral(0.2,0.5,p)-0.3)<1e-12);

% 参考样品Cref的完整强度校准，而非仅用弱散射近似2*eps。
[er,pr,cr]=particle_contrast(p.opt.d_ref,p.opt.n_ref,p);
assert(abs(2*er+er^2-p.opt.C_ref)<1e-12 && abs(cr-p.opt.C_ref)<1e-12);
assert(abs(pr-p.opt.phi0)<1e-12);

% 帧中心、曝光采样与有限半径球心边界。
p.mag.f=10; p.mag.phase0=0; p.mag.wave='square'; p.sim.gravity=false;
p.part.n_per_size=[12,12,12]; p.sim.n_frames=3000; p.sim.z_init=[];
tr=simulate_trajectories(p,91);
assert(isequal(size(tr.pos_exposure),[3000,36,3,4]));
assert(max(abs(mean(tr.t_exposure,2)-tr.t))<1e-12);
assert(max(abs(tr.pos_camera(:)-reshape(mean(tr.pos_exposure,4),[],1)))<1e-18);
z=reshape(tr.pos_exposure(:,:,3,:),[],36,1,4);
zlo=reshape(tr.dH/2,[1,36,1,1]); zhi=p.sim.H_z-zlo;
inside=(z>=zlo-1e-15) & (z<=zhi+1e-15);
assert(all(inside(:)));
% 均匀中点离散曝光的理论模糊项te*(1-1/nq^2)/3；只检查无磁力的y。
for k=1:3
    y=tr.pos_camera(:,tr.size_idx==k,2); dy=diff(y); dy2=y(3:end,:)-y(1:end-2,:);
    D=tr.D(find(tr.size_idx==k,1));
    expected=2*D*(tr.dt-p.cam.exposure*(1-1/p.sim.n_sub^2)/3);
    assert(abs(mean(dy(:).^2)/expected-1)<0.035,'曝光后的布朗方差不正确。');
    Dfit=(mean(dy2(:).^2)-mean(dy(:).^2))/(2*tr.dt);
    assert(abs(Dfit/D-1)<0.05,'曝光后MSD斜率应保持正确。');
end
p.cam.exposure=0; p.sim.n_frames=8;
ti=simulate_trajectories(p,3);
assert(isequal(ti.pos,ti.pos_camera) && ti.exposure_samples==1);

% 整盒轴向覆盖，不允许在u_max处无物理依据地硬关信号。
p.cam.Npx=33; p.cam.noise=false; p.cam.bg_amp=0; p.cam.illum_w=inf;
p.cam.bits=16; p.part.n_p=p.opt.n_ref;
tab=ipsf_lookup_table(p);
assert(tab.u(end)>=p.sim.H_z-p.opt.z_f-1e-15);
one=struct('d',p.opt.d_ref,'n_frames',1,'pos_wrap',reshape([0,0,p.opt.z_f],[1,1,3]));
[image,~]=render_iscat_video(one,p,tab);
expected=round(p.cam.photons*p.cam.gain*(1+p.opt.C_ref));
assert(abs(double(image(17,17))-expected)<=1,'相机实测参考对比度不匹配。');
p.part.n_p=2.40+0.50i; one.d=100e-9;
one.pos_wrap(1,1,3)=p.opt.z_f+2.7e-6;
[~,meta]=render_iscat_video(one,p,tab);
assert(meta.n_rendered==1,'u>2.5um处不应硬删仍可见的粒子。');

% 验证曝光积分是两幅瞬时强度的平均，不能变成复场平均。
a=one; b=one;
a.pos_wrap(1,1,3)=p.opt.z_f-0.1e-6;
b.pos_wrap(1,1,3)=p.opt.z_f+0.1e-6;
va=render_iscat_video(a,p,tab); vb=render_iscat_video(b,p,tab);
both=one; both.pos_exposure=cat(4,a.pos_wrap,b.pos_wrap); both.exposure=80e-6;
vi=render_iscat_video(both,p,tab);
delta=double(vi)-(double(va)+double(vb))/2;
assert(max(abs(delta(:)))<=1,'曝光平均光强与直接两帧平均不一致。');

% 非法照明参数必须在渲染前报错，不能由NaN光强静默生成全黑uint16图。
keys={'flicker_amplitude','flicker_frequency','flicker_phase'};
invalid={nan,1i,[0 1]};
for ik=1:numel(keys)
    for iv=1:numel(invalid)
        bad=p; bad.cam.(keys{ik})=invalid{iv}; rejected=false;
        try, render_iscat_video(one,bad,tab);
        catch err, rejected=contains(err.message,'周期照明参数无效'); end
        assert(rejected,'周期照明参数必须为有限实标量。');
    end
end
bad=p; bad.cam.flicker_amplitude=1; rejected=false;
try, render_iscat_video(one,bad,tab);
catch err, rejected=contains(err.message,'周期照明参数无效'); end
assert(rejected,'照明幅度绝对值必须小于1。');
bad=p; bad.cam.flicker_frequency=-1; rejected=false;
try, render_iscat_video(one,bad,tab);
catch err, rejected=contains(err.message,'周期照明参数无效'); end
assert(rejected,'照明频率不能为负数。');

% 完整PSF开域不允许悄悄套用周期包装；旧参数必须显式报错。
p.sim.boundary_xy='periodic';
rejected=false;
try
    render_iscat_video(one,p,tab);
catch err
    rejected=contains(err.message,'只支持open');
end
assert(rejected,'完整PSF必须拒绝未经实现的周期边界。');
fprintf('PASS test_physics_imaging: field integral, calibration, exposure, walls, PSF coverage, illumination input validation.\n');
end
