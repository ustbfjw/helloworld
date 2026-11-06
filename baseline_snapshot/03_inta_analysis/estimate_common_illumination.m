function [gain, info] = estimate_common_illumination(video,p)
%ESTIMATE_COMMON_ILLUMINATION 频谱辅助的共同乘性周期光照估计。
%   gain为F×1正值；调用方将原始每帧除以gain，再做原有静态平场和追踪。
%   仅从四个空间象限的稳健亮度序列识别共同窄带分量，不对粒子图像做
%   单频重建，不平滑轨迹，也不使用粒径/位置真值。无可信频率时严格返回1。
%   FFT与功率归一化均是已有方法；这里实现针对本项目的受限组合改进。
%   适用模型：I(x,y,t)=B(x,y)*g(t)*(1+C(x,y,t))+noise。
%   不校正非共同空间变化、宽带照明、加性暗电流或散粒噪声。
%   三个可选消融开关均默认true；关闭时只移除对应组件，其他阈值不变。
assert(isnumeric(video) && isreal(video) && ndims(video)==3, ...
    '输入必须是实数二维灰度视频，维度为行×列×帧。');
[ny,nx,F]=size(video);
assert(ny>=8 && nx>=8 && F>=16,'至少需要8×8像素和16帧。');
assert(isfield(p,'sim') && isfinite(p.sim.fps) && p.sim.fps>0,'帧率必须为正数。');
if ~isa(video,'uint16') && ~isa(video,'uint8')
    assert(all(isfinite(video(:))) && all(video(:)>=0),'原始光强必须有限且非负。');
end
Fs=p.sim.fps; cfg=struct();
if isfield(p,'ana') && isfield(p.ana,'illum'), cfg=p.ana.illum; end
step=opt(cfg,'sample_step',2); max_ref=opt(cfg,'reference_frames',256);
min_amp=opt(cfg,'min_amplitude',0.002); snr_min=opt(cfg,'noise_factor',8);
coh_min=opt(cfg,'min_coherence',0.95); spread_max=opt(cfg,'max_amplitude_cv',0.3);
max_lines=opt(cfg,'max_frequencies',3); protected_n=opt(cfg,'protect_harmonics',5);
protect_bins=opt(cfg,'protect_bins',2);
enable_spatial_gate=boolopt(cfg,'enable_spatial_gate');
enable_magnetic_guard=boolopt(cfg,'enable_magnetic_guard');
refine_frequency=boolopt(cfg,'refine_frequency');
assert(step>=1 && step==round(step) && max_ref>=1 && max_ref==round(max_ref));
assert(min_amp>0 && snr_min>0 && coh_min>0 && coh_min<=1 && spread_max>=0);
assert(max_lines>=1 && max_lines==round(max_lines) && protected_n>=0 && protected_n==round(protected_n));
assert(isscalar(protect_bins) && isfinite(protect_bins) && protect_bins>=0);
gain=ones(F,1);
info=struct('enabled',true,'applied',false,'accepted_frequencies',[], ...
    'enable_spatial_gate',enable_spatial_gate,'enable_magnetic_guard',enable_magnetic_guard, ...
    'refine_frequency',refine_frequency, ...
    'amplitudes',[],'coherence',[],'amplitude_cv',[],'tile_traces',zeros(F,4), ...
    'gain',gain,'raw_common',ones(F,1),'frequency_axis',[], ...
    'spectrum',[],'rejected_reason','未发现可信共同周期分量', ...
    'rejected_frequencies',[],'rejection_reasons',{{}},'frequency_resolution',Fs/F);

% 最多256个时间位置仅用于估计每个采样像素的参考亮度，非时间滤波输出。
V=single(video(1:step:ny,1:step:nx,:)); [sy,sx,~]=size(V);
assert(sy>=4 && sx>=4,'sample_step过大，无法形成四个空间象限。');
reference=unique(round(linspace(1,F,min(F,max_ref)))); B=median(V(:,:,reference),3);
traces=zeros(F,4); bounds={[1 floor(sy/2)],[floor(sy/2)+1 sy]};
xbounds={[1 floor(sx/2)],[floor(sx/2)+1 sx]}; k=0;
for iy=1:2
    for ix=1:2
        k=k+1; rr=bounds{iy}(1):bounds{iy}(2); cc=xbounds{ix}(1):xbounds{ix}(2);
        b=B(rr,cc); values=reshape(V(rr,cc,:),[],F); b=b(:);
        valid=isfinite(b) & b>0;
        if sum(valid)<4
            info.rejected_reason='某象限缺少正值背景像素'; return;
        end
        ratios=values(valid,:)./b(valid);
        traces(:,k)=double(median(ratios,1)).';
    end
end
clear V values ratios
scale=median(traces,1);
if any(~isfinite(scale) | scale<=0)
    info.rejected_reason='共同亮度参考值无效'; return;
end
traces=traces./scale;
info.tile_traces=traces; info.raw_common=median(traces,2);
t=(0:F-1)'/Fs; tc=t-mean(t); nuisance=[ones(F,1),tc];
Y=traces-1; Yd=Y-nuisance*(nuisance\Y);
common=median(Yd,2);
Z=fft(common); nfreq=floor(F/2)+1;
freq=(0:nfreq-1)'*Fs/F; spectrum=2*abs(Z(1:nfreq))/F;
spectrum(1)=abs(Z(1))/F;
if mod(F,2)==0, spectrum(end)=abs(Z(nfreq))/F; end
info.frequency_axis=freq; info.spectrum=spectrum;
df=Fs/F; fmin=max(opt(cfg,'min_frequency',5),2*df);
fmax=min(opt(cfg,'max_frequency',0.45*Fs),0.49*Fs);
peaks=find(spectrum(2:end-1)>=spectrum(1:end-2) & ...
    spectrum(2:end-1)>spectrum(3:end))+1;
peaks=peaks(spectrum(peaks)>=min_amp);
[~,order]=sort(spectrum(peaks),'descend'); peaks=peaks(order);
accepted=[];
for j=1:numel(peaks)
    bin=peaks(j); f0=freq(bin); reason='';
    if f0<fmin || f0>fmax
        reason='位于保留的低频或Nyquist邻域';
    elseif enable_magnetic_guard && isfield(p,'mag') && isfield(p.mag,'f') && p.mag.f>0 && ...
            any(abs(f0-p.mag.f*(1:protected_n))<=protect_bins*df)
        reason='位于磁驱动基频或受保护谐波邻域';
    else
        nb=max(2,bin-20):min(nfreq-1,bin+20); nb(abs(nb-bin)<=2)=[];
        floor_amp=median(spectrum(nb));
        if spectrum(bin)<snr_min*max(floor_amp,eps)
            reason='局部频谱信噪比不足';
        end
    end
    if ~isempty(reason)
        info.rejected_frequencies(end+1)=f0;
        info.rejection_reasons{end+1}=reason; continue;
    end
    % 非整周期采样会泄漏到相邻DFT格点；只在该峰±半个格点内细化。
    fhat=f0;
    if refine_frequency
        lo=max(fmin,f0-df/2); hi=min(fmax,f0+df/2);
        fun=@(f) residual_energy(f,t,nuisance,common);
        fhat=fminbnd(fun,lo,hi,optimset('Display','off','TolX',1e-7));
    end
    if enable_magnetic_guard && isfield(p,'mag') && isfield(p.mag,'f') && p.mag.f>0 && ...
            any(abs(fhat-p.mag.f*(1:protected_n))<=protect_bins*df)
        info.rejected_frequencies(end+1)=fhat;
        info.rejection_reasons{end+1}='频率复查落入磁驱动保护区'; continue;
    end
    A=[nuisance,cos(2*pi*fhat*t),sin(2*pi*fhat*t)]; coef=A\Y;
    z=coef(end-1,:)-1i*coef(end,:); amp=abs(z);
    coherence=abs(mean(z))/max(mean(amp),eps);
    cv=std(amp)/max(mean(amp),eps);
    reason='';
    if median(amp)<min_amp
        reason='四象限共同振幅不足';
    elseif enable_spatial_gate && coherence<coh_min
        reason='四象限相位一致性不足';
    elseif enable_spatial_gate && cv>spread_max
        reason='四象限幅度一致性不足';
    end
    if ~isempty(reason)
        info.rejected_frequencies(end+1)=fhat;
        info.rejection_reasons{end+1}=reason; continue;
    end
    if any(abs(fhat-accepted)<=protect_bins*df), continue; end
    accepted(end+1)=fhat; %#ok<AGROW>
    info.coherence(end+1)=coherence; info.amplitude_cv(end+1)=cv;
    if numel(accepted)>=max_lines, break; end
end
if isempty(accepted), return; end

% 联合估计周期项；截距和线性趋势仅为拟合干扰项，绝不从数据中移除。
A=nuisance;
for f=accepted, A=[A,cos(2*pi*f*t),sin(2*pi*f*t)]; end %#ok<AGROW>
c=A\(info.raw_common-1); periodic=A(:,3:end)*c(3:end);
periodic=periodic-mean(periodic); proposed=1+periodic;
if any(~isfinite(proposed)) || min(proposed)<=0.5 || max(proposed)>1.5
    info.rejected_reason='估计增益超出安全正值范围，拒绝校正';
    info.rejected_frequencies=[info.rejected_frequencies accepted];
    info.rejection_reasons=[info.rejection_reasons,repmat({info.rejected_reason},1,numel(accepted))];
    info.coherence=[]; info.amplitude_cv=[];
    return;
end
gain=proposed; info.gain=gain; info.applied=true;
info.accepted_frequencies=accepted;
info.amplitudes=hypot(c(3:2:end),c(4:2:end)).';
info.rejected_reason='';
end

function value=opt(s,key,default)
if isfield(s,key), value=s.(key); else, value=default; end
end

function value=boolopt(s,key)
value=opt(s,key,true);
assert(islogical(value) && isscalar(value),'iNTA:illum:InvalidBooleanFlag', ...
    '%s必须是标量logical（true或false）。',key);
end

function rss=residual_energy(f,t,nuisance,y)
A=[nuisance,cos(2*pi*f*t),sin(2*pi*f*t)]; residual=y-A*(A\y);
rss=sum(residual.^2);
end
