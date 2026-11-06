function metrics=test_common_illumination(outdir)
%TEST_COMMON_ILLUMINATION 公共周期光照估计的已知合成输入回归。
% 非均匀静态平场、局部运动前景和噪声均由测试构造，真增益只供误差评价。
% 该测试不替代反射iSCAT图像的完整盲检测/粒径验收。
if nargin<1, outdir=''; end
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root); setup_paths(root);
p=sim_params('full'); p.sim.fps=5000; p.mag.f=10;
N=48; F=5000; t=(0:F-1)'/p.sim.fps;
[X,Y]=meshgrid(linspace(-1,1,N));
B=800*(0.65+0.35*exp(-2*(X.^2+Y.^2))).*(1+0.03*cos(7*X).*sin(5*Y));
rs=RandStream('mt19937ar','Seed',50813); clean=zeros(N,N,F,'single');
for f=1:F
    xc=0.55*sin(2*pi*1.7*t(f)); yc=0.5*cos(2*pi*2.3*t(f));
    foreground=0.12*exp(-((X-xc).^2+(Y-yc).^2)/0.004);
    clean(:,:,f)=single(B.*(1+foreground)+2*randn(rs,N));
end
[g0,i0]=estimate_common_illumination(clean,p);
assert(isequal(g0,ones(F,1)) && ~i0.applied,'无周期照明时必须为严格恒等。');
metrics=struct('scope','Synthetic common multiplicative illumination, not end-to-end particle sizing', ...
    'zero_flicker_identity',true);
for k=1:2
    freqs=[420 423.7]; freq=freqs(k);
    truth=1+0.05*cos(2*pi*freq*t+0.37);
    video=clean.*reshape(single(truth),1,1,F);
    before=rng; [gain,info]=estimate_common_illumination(video,p); after=rng;
    assert(isequal(before,after),'光照估计不能污染仿真随机流。');
    assert(info.applied && numel(info.accepted_frequencies)==1);
    target=truth/mean(truth); error_rms=sqrt(mean((gain-target).^2));
    frequency_error=abs(info.accepted_frequencies(1)-freq);
    assert(error_rms<3e-4 && frequency_error<0.02,'增益或非整周期频率恢复失败。');
    corrected=video./reshape(single(gain),1,1,F);
    rel_before=norm(double(video(:)-clean(:)))/norm(double(clean(:)));
    ref=clean*mean(truth);
    rel_after=norm(double(corrected(:)-ref(:)))/norm(double(ref(:)));
    assert(rel_after<0.02*rel_before,'校正应显著降低公共闪烁误差并保留逐帧前景。');
    % 生成器的真光照参数不得影响估计。
    pp=p; pp.cam.flicker_amplitude=0.9; pp.cam.flicker_frequency=999; pp.cam.flicker_phase=1.2;
    [other,~]=estimate_common_illumination(video,pp);
    assert(isequal(gain,other),'估计器读取了生成器光照真值。');
    row=struct('frequency',freq,'estimated_frequency',info.accepted_frequencies(1), ...
        'gain_RMS_error',error_rms,'relative_image_error_before',rel_before, ...
        'relative_image_error_after',rel_after,'amplitude',info.amplitudes(1));
    metrics.recovery(k)=row;
    fprintf('ILLUMINATION f=%.3f estimated=%.6f gain_RMS=%.8f image_before=%.6f image_after=%.8f\n', ...
        freq,info.accepted_frequencies(1),error_rms,rel_before,rel_after);
end

% 相位一致但幅度不同的四象限，不可冒充全视场共同乘性光照。
video=clean; a=[0.08 0.02 0.02 0.02]; q=0;
for iy=1:2
    rr=(iy-1)*N/2+(1:N/2);
    for ix=1:2
        q=q+1; cc=(ix-1)*N/2+(1:N/2);
        video(rr,cc,:)=clean(rr,cc,:).*reshape(single(1+a(q)*cos(2*pi*420*t)),1,1,F);
    end
end
[g,info]=estimate_common_illumination(video,p);
assert(isequal(g,ones(F,1)) && ~info.applied,'空间幅度不一致应拒绝。');
metrics.inconsistent_tiles_rejected=true;

% 磁驱动基频附近不盲目估计为照明，不能删去真正同步的颗粒信息。
video=clean.*reshape(single(1+0.05*cos(2*pi*10*t)),1,1,F);
[g,info]=estimate_common_illumination(video,p);
assert(isequal(g,ones(F,1)) && ~info.applied);
assert(~isempty(info.rejection_reasons),'应记录受保护频率的拒绝原因。');
metrics.magnetic_frequency_protected=true;

% 相同幅度但一个象限相位反向，同样不能作为共同照明校正。
video=clean.*reshape(single(1+0.05*cos(2*pi*420*t)),1,1,F);
video(1:N/2,1:N/2,:)=clean(1:N/2,1:N/2,:).* ...
    reshape(single(1-0.05*cos(2*pi*420*t)),1,1,F);
[g,info]=estimate_common_illumination(video,p);
assert(isequal(g,ones(F,1)) && ~info.applied);
metrics.inconsistent_phase_rejected=true;

% 两条非整周期频率同时存在，联合回归应恢复其叠加而非强制单频。
truth=1+0.03*cos(2*pi*423.7*t)+0.02*sin(2*pi*731.2*t);
[g,info]=estimate_common_illumination(clean.*reshape(single(truth),1,1,F),p);
assert(info.applied && numel(info.accepted_frequencies)==2);
metrics.multifrequency_gain_RMS=sqrt(mean((g-truth/mean(truth)).^2));
assert(metrics.multifrequency_gain_RMS<5e-4);

% 过大的共同调制虽能识别，仍拒绝可能不稳定的归一化。
[g,info]=estimate_common_illumination(clean.*reshape(single(1+0.6*cos(2*pi*420*t)),1,1,F),p);
assert(isequal(g,ones(F,1)) && ~info.applied && ~isempty(info.rejected_reason));
metrics.extreme_gain_rejected=true;

% 非法输入必须报错；全零参考仅退回恒等并解释原因。
failed=false; try, estimate_common_illumination(zeros(N),p); catch, failed=true; end
assert(failed);
bad=clean; bad(1)=-1; failed=false;
try, estimate_common_illumination(bad,p); catch, failed=true; end
assert(failed);
[g,info]=estimate_common_illumination(zeros(N,N,16,'single'),p);
assert(isequal(g,ones(16,1)) && ~info.applied && ~isempty(info.rejected_reason));
metrics.input_validation=true; metrics.pass=true;
if ~isempty(outdir)
    if ~exist(outdir,'dir'), mkdir(outdir); end
    fid=fopen(fullfile(outdir,'common_illumination_tests.json'),'w','n','UTF-8');
    fprintf(fid,'%s',jsonencode(metrics,PrettyPrint=true)); fclose(fid);
end
fprintf('PASS test_common_illumination: identity, on/off-bin gain recovery, spatial consistency, magnetic guard, no truth leakage.\n');
end
