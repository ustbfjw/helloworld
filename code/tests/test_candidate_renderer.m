function metrics = test_candidate_renderer(outdir)
%TEST_CANDIDATE_RENDERER 保留形状生成与冻结渲染器的独立数值回归。
%   三组各100帧验证global/left_half逐位一致；新形状用静止颗粒的
%   已知曝光平均增益及独立逐像素/逐曝光参考检查，不评价估计器效果。
if nargin<1, outdir=''; end
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root); setup_paths(root);
saved_rng=rng; restore_rng=onCleanup(@() rng(saved_rng)); %#ok<NASGU>
timer=tic; p0=sim_params('full'); p0.sim.n_frames=100;
protocol=targeted_experiment_protocol(); conditions=protocol.conditions;
new_modes={'gaussian_patch','tilted_smooth','oblique_step'};
for k=1:3
    conditions(6+k)=struct('name',new_modes{k},'amplitude',0.08,'frequency',423.7, ...
        'phase',0.37,'spatial_mode',new_modes{k});
end
tab=ipsf_lookup_table(p0); groups={'mono30','mono100','mono200'};
metrics=struct('scope','Renderer correctness only; held-out shapes are not estimator acceptance', ...
    'frames_per_group',100,'camera_Npx',p0.cam.Npx,'exposure_samples',p0.sim.n_sub, ...
    'baseline_exact_cases',0,'stationary_gain_cases',0,'pixel_reference_cases',0, ...
    'wrong_exposure_order_different_pixels',0,'wrong_noise_order_different_pixels',0, ...
    'maximum_stationary_adu_error',0,'maximum_mask_error',0);
for ig=1:3
    p=apply_scenario(p0,groups{ig}); traj=simulate_trajectories(p,820000+ig);
    rng(920000+ig,'twister'); initial_rng=rng;
    [videos,metas]=render_candidate_conditions(traj,p,tab,conditions); candidate_rng=rng;
    rng(initial_rng); [old_videos,old_metas]=render_targeted_conditions(traj,p,tab,conditions(1:6));
    assert(isequal(candidate_rng,rng),'candidate:RngMismatch','候选渲染改变随机数消耗。');
    for k=1:6
        assert(isequal(videos{k},old_videos{k}) && isequaln(metas{k},old_metas{k}), ...
            'candidate:BaselineMismatch','旧条件的影像或完整元数据不同。');
        metrics.baseline_exact_cases=metrics.baseline_exact_cases+1;
    end
    for k=7:9
        mask=reference_mask(p.cam.Npx,conditions(k).spatial_mode);
        mask_error=max(abs(mask(:)-metas{k}.gain_mask(:)));
        assert(mask_error<2e-15 && all(mask(:)>=0 & mask(:)<=1));
        metrics.maximum_mask_error=max(metrics.maximum_mask_error,mask_error);
        truth_offset=mean(conditions(k).amplitude*sin(2*pi*conditions(k).frequency*traj.t_exposure+conditions(k).phase),2);
        assert(max(abs(truth_offset-metas{k}.gain_offset))<5e-16);
    end
    clear videos old_videos metas old_metas

    % 将所有曝光位置固定，光场不再随曝光变化，此时每像素真强度严格
    % 等于静态强度乘曝光平均增益。只允许两个独立ADC量化引入的误差。
    stationary=traj;
    stationary.pos_exposure=repmat(traj.pos_exposure(1,:,:,1),[traj.n_frames,1,1,p.sim.n_sub]);
    ps=p; ps.cam.noise=false; ps.cam.bits=16; ps.cam.gain=100/ps.cam.photons;
    rng(initial_rng); [static_video,static_meta]=render_iscat_video(stationary,ps,tab);
    rng(initial_rng); [shape_videos,shape_metas]=render_candidate_conditions(stationary,ps,tab,conditions(7:9));
    assert(static_meta.n_sat==0);
    for k=1:3
        assert(shape_metas{k}.n_sat==0);
        mask=reference_mask(p.cam.Npx,new_modes{k});
        for f=1:p.sim.n_frames
            delta=0;
            for q=1:p.sim.n_sub
                delta=delta+0.08*sin(2*pi*423.7*traj.t_exposure(f,q)+0.37)/p.sim.n_sub;
            end
            known_gain=1+mask*delta;
            error=abs(double(shape_videos{k}(:,:,f))-double(static_video(:,:,f)).*known_gain);
            assert(all(error(:)<=0.5*(1+known_gain(:))+2e-12), ...
                'candidate:KnownGainMismatch','静止颗粒的逐像素已知增益检验失败。');
            metrics.maximum_stationary_adu_error=max(metrics.maximum_stationary_adu_error,max(error(:)));
        end
        metrics.stationary_gain_cases=metrics.stationary_gain_cases+1;
    end
    fprintf('CANDIDATE_RENDERER %s: six frozen cases exact; three stationary gain maps verified.\n',groups{ig});

    % 独立小图参考：改变曝光内z使光强与照明相关，防止检验只覆盖
    % mean(I)*mean(g)也能成立的平凡情况。逐像素插值/求和不调用候选helper。
    pr=p; pr.cam.Npx=12; pr.cam.bits=16; pr.cam.gain=1000/pr.cam.photons;
    pr.cam.exposure=160e-6; pr.sim.n_frames=6;
    one=reference_trajectory(pr,ig);
    test_conditions=conditions(7:9);
    for k=1:3, test_conditions(k).amplitude=0.35; test_conditions(k).frequency=1700; end
    rng(930000+ig,'twister'); initial_rng=rng;
    [observed,~]=render_candidate_conditions(one,pr,tab,test_conditions); candidate_rng=rng;
    [reference,wrong_exposure,wrong_noise,reference_rng]=pixel_reference(one,pr,tab,test_conditions,initial_rng);
    assert(isequal(candidate_rng,reference_rng),'独立像素参考随机流不一致。');
    for k=1:3
        assert(isequal(observed{k},reference{k}), ...
            'candidate:PixelReferenceMismatch','独立逐像素噪声/曝光/ADC参考不一致。');
        metrics.pixel_reference_cases=metrics.pixel_reference_cases+1;
        metrics.wrong_exposure_order_different_pixels=metrics.wrong_exposure_order_different_pixels+ ...
            nnz(observed{k}~=wrong_exposure{k});
        metrics.wrong_noise_order_different_pixels=metrics.wrong_noise_order_different_pixels+ ...
            nnz(observed{k}~=wrong_noise{k});
    end
end
assert(metrics.wrong_exposure_order_different_pixels>0 && metrics.wrong_noise_order_different_pixels>0, ...
    '独立参考必须能够检出错误曝光或噪声顺序，不能是平凡同一影像。');

% 零幅值时新形状也应逐位退化到静态；非法模式在抽随机数之前拒绝。
for k=1:3, test_conditions(k).amplitude=0; end
rng(initial_rng); zero_videos=render_candidate_conditions(one,pr,tab,test_conditions); zero_rng=rng;
pr.cam.flicker_amplitude=0; rng(initial_rng); static_video=render_iscat_video(one,pr,tab);
assert(isequal(zero_rng,rng));
for k=1:3, assert(isequal(zero_videos{k},static_video)); end
bad=test_conditions(1); bad.spatial_mode='quadpoly'; before=rng; rejected=false;
try, render_candidate_conditions(one,pr,tab,bad); catch err, rejected=strcmp(err.identifier,'targeted:InvalidConditions'); end
assert(rejected && isequal(before,rng),'未知空间模式必须在消耗随机数前拒绝。');
metrics.baseline_metadata_rng_exact=true; metrics.zero_amplitude_identity=true;
metrics.unknown_shape_rejected=true; metrics.pass=true; metrics.elapsed_seconds=toc(timer);
if ~isempty(outdir)
    if ~exist(outdir,'dir'), mkdir(outdir); end
    fid=fopen(fullfile(outdir,'candidate_renderer_tests.json'),'w','n','UTF-8');
    assert(fid>=0,'测试指标JSON无法写入。');
    close_file=onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%s',jsonencode(metrics,PrettyPrint=true));
end
fprintf('PASS test_candidate_renderer: %d frozen cases, %d known-gain cases, %d independent pixel references (%.2f s).\n', ...
    metrics.baseline_exact_cases,metrics.stationary_gain_cases,metrics.pixel_reference_cases,metrics.elapsed_seconds);
end

function mask=reference_mask(N,mode)
% 独立逐像素坐标计算，不使用生成器的meshgrid或元数据掩膜。
mask=zeros(N);
for row=1:N
    for col=1:N
        x=(col-(N+1)/2)/N; y=(row-(N+1)/2)/N;
        switch mode
            case 'gaussian_patch', mask(row,col)=exp(-((x-0.18)^2+(y+0.12)^2)/(2*0.28^2));
            case 'tilted_smooth', mask(row,col)=0.55+0.20*sin(2*pi*(x+0.17))+0.15*cos(2*pi*(y-0.11));
            case 'oblique_step', mask(row,col)=double(x+0.6*y>0.08);
        end
    end
end
end

function tr=reference_trajectory(p,size_index)
F=p.sim.n_frames; nq=p.sim.n_sub;
tr=struct('d',p.part.d_list(size_index),'n_frames',F,'exposure',p.cam.exposure);
tr.t_exposure=(1:F)'/p.sim.fps+(((1:nq)-0.5)/nq-0.5)*p.cam.exposure;
tr.pos_exposure=zeros(F,1,3,nq); z_offset=[-0.35,-0.1,0.1,0.35]*1e-6;
for f=1:F
    for q=1:nq
        tr.pos_exposure(f,1,:,q)=reshape([0.04e-6*sin(f+q),0.03e-6*cos(f-q),p.opt.z_f+z_offset(q)],1,1,3);
    end
end
end

function [videos,wrong_exposure,wrong_noise,ending_rng]=pixel_reference(tr,p,tab,conditions,initial_rng)
% 物理独立参考：调用冻结生产renderer只获取相同静态背景；复场沿u/rho
% 使用MATLAB通用interp1逐像素求值，不调用ipsf_eval或候选渲染器。
% 正确与两个故意错误的模型使用同一zshot/zread，验证测试可辨别运算顺序。
pn=p; pn.cam.noise=false; pn.cam.flicker_amplitude=0;
rng(initial_rng); [~,background]=render_iscat_video(tr,pn,tab);
N=p.cam.Npx; F=tr.n_frames; nq=p.sim.n_sub; K=numel(conditions);
[epsilon,phase]=particle_contrast(tr.d,p.part.n_p,p); amplitude=epsilon*exp(1i*phase);
videos=cell(1,K); wrong_exposure=cell(1,K); wrong_noise=cell(1,K); masks=cell(1,K);
for k=1:K
    videos{k}=zeros(N,N,F,'uint16'); wrong_exposure{k}=videos{k}; wrong_noise{k}=videos{k};
    masks{k}=reference_mask(N,conditions(k).spatial_mode);
end
rng(initial_rng); randn(N); % 冻结生产模型只用此一次随机矩阵生成静态背景。
for f=1:F
    shot=randn(N); read=randn(N); profiles=zeros(numel(tab.rho),nq);
    for q=1:nq
        z=tr.pos_exposure(f,1,3,q)-p.opt.z_f;
        profiles(:,q)=interp1(tab.u,tab.H.',z,'linear').';
    end
    for row=1:N
        for col=1:N
            x=(col-(N+1)/2)*p.cam.px; y=(row-(N+1)/2)*p.cam.px;
            intensity=zeros(1,nq);
            for q=1:nq
                radius=norm([x-tr.pos_exposure(f,1,1,q),y-tr.pos_exposure(f,1,2,q)]);
                field=1+amplitude*interp1(tab.rho,profiles(:,q),radius,'linear');
                intensity(q)=real(field)^2+imag(field)^2;
            end
            base=background.illum(row,col)*background.bg_static(row,col);
            for k=1:K
                c=conditions(k);
                gain=1+masks{k}(row,col)*c.amplitude*sin(2*pi*c.frequency*tr.t_exposure(f,:)+c.phase);
                expected_e=p.cam.photons*sum(intensity.*gain)/nq*base;
                observed_e=expected_e+sqrt(expected_e)*shot(row,col)+p.cam.read_e*read(row,col);
                videos{k}(row,col,f)=adc(observed_e,p);
                % 故意错误一：先积分I，后乘平均gain，丢失曝光内相关项。
                wrong_e=p.cam.photons*mean(intensity)*mean(gain)*base;
                wrong_exposure{k}(row,col,f)=adc(wrong_e+sqrt(wrong_e)*shot(row,col)+p.cam.read_e*read(row,col),p);
                % 故意错误二：保持无噪声均值与正确模型严格相同，仅将
                % 散粒/读出噪声错误地放在等效增益之前，独立检查噪声顺序。
                mean_gain=mean(gain); pre_gain_e=expected_e/mean_gain;
                wrong_noise{k}(row,col,f)=adc((pre_gain_e+sqrt(pre_gain_e)*shot(row,col)+p.cam.read_e*read(row,col))*mean_gain,p);
            end
        end
    end
end
ending_rng=rng;
end

function pixel=adc(electrons,p)
pixel=uint16(min(max(round(electrons*p.cam.gain),0),2^p.cam.bits-1));
end
