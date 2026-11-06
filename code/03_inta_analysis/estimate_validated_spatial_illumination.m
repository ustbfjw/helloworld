function [model,info]=estimate_validated_spatial_illumination(video,p)
%ESTIMATE_VALIDATED_SPATIAL_ILLUMINATION 留出验证的低阶空间周期增益候选。
% 局部模型仅用于磁频保护区外；全局回退严格调用原估计器，不改其结果。
% 训练/验证为每个4x4网格内不重叠的两个对角小块，块间留空隙。
% 空间中位数与多项式基函数均值并非严格等价，由留出评分约束此近似。
% 留出像素不保证无颗粒；本函数不能解决磁同频照明与颗粒运动的可辨识性。
% is_spatial=true时：gain(:,:,f)=1+sum(coeff.*basis(f,:),3)，不要再除global_gain。
% info.applied/accepted_frequencies描述实际校正；local_frequencies另存局部拟合频率。
assert(isnumeric(video) && isreal(video) && ndims(video)==3, ...
    'iNTA:spatial:InvalidVideo','输入必须为实数灰度视频。');
[ny,nx,F]=size(video);
assert(ny>=32 && nx>=32 && F>=16,'iNTA:spatial:InvalidVideo', ...
    '空间候选至少需要32x32像素和16帧。');
assert(all(isfinite(video(:))) && all(video(:)>=0),'iNTA:spatial:InvalidVideo', ...
    '原始光强必须有限且非负。');
assert(isstruct(p) && isfield(p,'sim') && isfield(p.sim,'fps') && ...
    scalar_number(p.sim.fps) && p.sim.fps>0,'iNTA:spatial:InvalidParameters','帧率必须为正标量。');
assert(isfield(p,'mag') && isfield(p.mag,'f') && scalar_number(p.mag.f) && p.mag.f>=0, ...
    'iNTA:spatial:InvalidParameters','必须给出实验磁频率的非负有限标量。');
cfg=struct();
if isfield(p,'ana') && isfield(p.ana,'illum'), cfg=p.ana.illum; end
assert(isstruct(cfg) && isscalar(cfg),'iNTA:spatial:InvalidParameters','illum必须为标量结构。');
guard=opt(cfg,'enable_magnetic_guard',true);
protected_n=opt(cfg,'protect_harmonics',5); protected_bins=opt(cfg,'protect_bins',2);
assert(islogical(guard) && isscalar(guard) && guard && ...
    scalar_number(protected_n) && protected_n==5 && scalar_number(protected_bins) && protected_bins==2, ...
    'iNTA:spatial:ProtectedConfiguration','空间候选禁止关闭磁频保护或改变5倍频/2格点保护范围。');
min_amp=opt(cfg,'min_amplitude',0.002); snr_min=opt(cfg,'noise_factor',8);
max_lines=opt(cfg,'max_frequencies',3); max_ref=opt(cfg,'reference_frames',256);
min_freq=opt(cfg,'min_frequency',5); max_freq=opt(cfg,'max_frequency',0.45*p.sim.fps);
assert(scalar_number(min_amp) && min_amp>0 && scalar_number(snr_min) && snr_min>0 && ...
    scalar_number(max_lines) && max_lines>=1 && max_lines==round(max_lines) && ...
    scalar_number(max_ref) && max_ref>=1 && max_ref==round(max_ref) && ...
    scalar_number(min_freq) && min_freq>=0 && scalar_number(max_freq) && max_freq>min_freq, ...
    'iNTA:spatial:InvalidParameters','照明频谱参数无效。');
max_lines=min(3,max_lines); max_ref=min(256,max_ref);
[fallback_gain,fallback_info]=estimate_common_illumination(video,p);
model=struct('global_gain',fallback_gain,'spatial_coefficients',[], ...
    'temporal_basis',[],'is_spatial',false);
info=struct('method','validated_spatial_polynomial','applied',fallback_info.applied,'is_spatial',false, ...
    'reason','未发现可信保护区外局部分量','fallback_info',fallback_info, ...
    'candidate_frequencies',[],'local_frequencies',[], ...
    'accepted_frequencies',fallback_info.accepted_frequencies,'rejected_frequencies',[], ...
    'rejection_reasons',{{}},'partition',struct(),'training_coefficients',[], ...
    'validation_coefficients',[],'polynomial_coefficients',[],'score',struct(), ...
    'frequency_resolution',p.sim.fps/F,'frequency_axis',[],'training_rms_spectrum',[], ...
    'gain_range',[NaN NaN],'rules',struct('grid',[4 4],'polynomial_order',2, ...
    'minimum_relative_energy_reduction',0.5,'tile_amplitude_tolerance',0.002, ...
    'minimum_global_residual_amplitude_rms',0.002, ...
    'minimum_active_tiles',2,'minimum_valid_pixels_per_partition',4, ...
    'protect_harmonics',5,'protect_bins',2,'max_frequencies',max_lines));

% 固定空间划分与坐标：每tile四个内缩角块，对角块训练、另一对角验证。
[X,Y]=meshgrid(linspace(-1,1,nx),linspace(-1,1,ny));
phi=[ones(ny*nx,1),X(:),Y(:),X(:).^2,X(:).*Y(:),Y(:).^2];
redges=round(linspace(0,ny,5)); cedges=round(linspace(0,nx,5));
train_ids=cell(16,1); validation_ids=cell(16,1);
train=zeros(F,16); validation=zeros(F,16); Gtrain=zeros(16,6); Gval=zeros(16,6);
counts=zeros(16,2); references=unique(round(linspace(1,F,min(F,max_ref))));
flat=reshape(video,ny*nx,F); tile=0;
for iy=1:4
    rr=redges(iy)+(1:redges(iy+1)-redges(iy));
    [rlo,rhi]=inner_blocks(rr);
    for ix=1:4
        tile=tile+1; cc=cedges(ix)+(1:cedges(ix+1)-cedges(ix));
        [clo,chi]=inner_blocks(cc);
        train_ids{tile}=[rectangle_ids(rlo,clo,ny,nx);rectangle_ids(rhi,chi,ny,nx)];
        validation_ids{tile}=[rectangle_ids(rlo,chi,ny,nx);rectangle_ids(rhi,clo,ny,nx)];
        [train(:,tile),Gtrain(tile,:),counts(tile,1)]=tile_trace(flat,train_ids{tile},references,phi);
        [validation(:,tile),Gval(tile,:),counts(tile,2)]=tile_trace(flat,validation_ids{tile},references,phi);
    end
end
info.partition=struct('row_edges',redges,'column_edges',cedges, ...
    'training_linear_indices',{train_ids},'validation_linear_indices',{validation_ids}, ...
    'valid_pixel_counts',counts,'training_basis_means',Gtrain,'validation_basis_means',Gval, ...
    'coordinate_convention','x=column,y=row,both normalized to [-1,1]', ...
    'reference_frames',references,'block_rule','inner 1/8:3/8 and 5/8:7/8; diagonal train, antidiagonal validation');
if any(counts(:)<4) || any(~isfinite(train(:))) || any(~isfinite(validation(:))) || rank(Gtrain)<6
    info.reason='训练或验证背景支持不足，保留原全局结果'; return;
end
Fs=p.sim.fps; t=(0:F-1)'/Fs; nuisance=[ones(F,1),t-mean(t)];
Ytrain=train-1; Yd=Ytrain-nuisance*(nuisance\Ytrain);
Z=fft(Yd); nfreq=floor(F/2)+1; frequency=(0:nfreq-1)'*Fs/F;
spectra=2*abs(Z(1:nfreq,:))/F; spectra(1,:)=abs(Z(1,:))/F;
if mod(F,2)==0, spectra(end,:)=abs(Z(nfreq,:))/F; end
spectrum=sqrt(mean(spectra.^2,2));
info.frequency_axis=frequency; info.training_rms_spectrum=spectrum;
df=Fs/F; fmin=max(min_freq,2*df); fmax=min(max_freq,0.49*Fs);
peaks=find(spectrum(2:end-1)>=spectrum(1:end-2) & spectrum(2:end-1)>spectrum(3:end))+1;
peaks=peaks(spectrum(peaks)>=min_amp);
[~,order]=sort(spectrum(peaks),'descend'); peaks=peaks(order);
accepted=[];
for j=1:numel(peaks)
    bin=peaks(j); f0=frequency(bin); reason='';
    info.candidate_frequencies(end+1)=f0;
    if f0<fmin || f0>fmax
        reason='位于保留的低频或Nyquist邻域';
    elseif protected(f0,p.mag.f,df)
        reason='粗频率位于磁驱动保护区';
    else
        nb=max(2,bin-20):min(nfreq-1,bin+20); nb(abs(nb-bin)<=2)=[];
        if isempty(nb) || spectrum(bin)<snr_min*max(median(spectrum(nb)),eps)
            reason='局部频谱信噪比不足';
        end
    end
    if ~isempty(reason)
        info.rejected_frequencies(end+1)=f0; info.rejection_reasons{end+1}=reason; continue;
    end
    [~,strongest]=max(spectra(bin,:));
    fun=@(f) residual_energy(f,t,nuisance,Yd(:,strongest));
    fhat=fminbnd(fun,max(fmin,f0-df/2),min(fmax,f0+df/2), ...
        optimset('Display','off','TolX',1e-7));
    if protected(fhat,p.mag.f,df)
        info.rejected_frequencies(end+1)=fhat;
        info.rejection_reasons{end+1}='细化频率位于磁驱动保护区'; continue;
    end
    if any(abs(fhat-accepted)<=2*df), continue; end
    accepted(end+1)=fhat; %#ok<AGROW>
    if numel(accepted)>=max_lines, break; end
end
info.local_frequencies=accepted;
if isempty(accepted), return; end

% 频率与所有模型参数只由训练像素决定。验证投影仅用于测量残差，不更新模型。
T=zeros(F,2*numel(accepted));
for k=1:numel(accepted)
    T(:,2*k-1)=cos(2*pi*accepted(k)*t); T(:,2*k)=sin(2*pi*accepted(k)*t);
end
A=[nuisance,T]; fit_train=A\Ytrain; coefficients=fit_train(3:end,:).';
fit_global=A\(median(train,2)-1); common=fit_global(3:end).';
beta=Gtrain\coefficients; predicted=Gval*beta;
fit_validation=A\(validation-1); observed=fit_validation(3:end,:).';
rl=sqrt(sum((observed-predicted).^2,2));
rg=sqrt(sum((observed-common).^2,2)); ru=sqrt(sum(observed.^2,2));
El=mean(rl.^2); Eg=mean(rg.^2); Eu=mean(ru.^2);
active=ru>=min_amp;
info.training_coefficients=coefficients; info.validation_coefficients=observed;
info.polynomial_coefficients=beta;
info.score=struct('local_residual_amplitude',rl,'global_residual_amplitude',rg, ...
    'uncorrected_amplitude',ru,'local_energy',El,'global_energy',Eg, ...
    'uncorrected_energy',Eu,'active_validation_tiles',find(active).', ...
    'global_residual_amplitude_rms',sqrt(Eg), ...
    'predicted_validation_coefficients',predicted,'global_coefficients',common, ...
    'relative_energy_reduction',NaN, ...
    'amplitude_definition','sqrt(mean_tiles(sum_frequencies(cos_coefficient^2+sin_coefficient^2))); not time-domain RMS', ...
    'score_note','固定开发门槛，不是统计显著性检验');
if Eg>0, info.score.relative_energy_reduction=1-El/Eg; end
if sum(active)<2
    info.reason='至少两个验证区的周期幅度支持不足，保留原全局结果'; return;
end
% 这里0.002指各频率cos/sin残差向量的跨区RMS，不是正弦时间RMS。
if sqrt(Eg)<0.002
    info.reason='全局验证残差幅度RMS小于0.002，无需增加局部模型'; return;
end
if ~(Eg>0 && El<=0.5*Eg)
    info.reason='局部模型未达到验证能量降低50%的门槛，保留原全局结果'; return;
end
if any(rl>rg+0.002)
    info.reason='存在验证区比全局模型恶化超过0.002幅度，保留原全局结果'; return;
end
if El>Eu || any(rl>ru+0.002)
    info.reason='局部模型比未校正验证结果恶化，保留原全局结果'; return;
end
T=T-mean(T,1); spatial=phi*beta;
gain_min=inf; gain_max=-inf;
for start=1:128:F
    idx=start:min(F,start+127); gain=1+spatial*T(idx,:).';
    gain_min=min(gain_min,min(gain(:))); gain_max=max(gain_max,max(gain(:)));
    if any(~isfinite(gain(:))) || gain_min<=0.5 || gain_max>1.5
        info.gain_range=[gain_min gain_max];
        info.reason='局部增益超出全图全时安全正值范围，保留原全局结果'; return;
    end
end
model.global_gain=ones(F,1);
model.spatial_coefficients=reshape(spatial,ny,nx,[]); model.temporal_basis=T; model.is_spatial=true;
info.applied=true; info.is_spatial=true; info.gain_range=[gain_min gain_max];
info.accepted_frequencies=accepted;
info.reason='局部模型通过固定空间留出门槛，仅应用局部增益';
end

function [lo,hi]=inner_blocks(indices)
n=numel(indices);
lo=indices(floor(n/8)+1:floor(3*n/8)); hi=indices(floor(5*n/8)+1:floor(7*n/8));
end

function ids=rectangle_ids(rr,cc,ny,nx)
[rows,cols]=ndgrid(rr,cc); ids=sub2ind([ny nx],rows(:),cols(:));
end

function [trace,basis,n]=tile_trace(flat,ids,references,phi)
values=single(flat(ids,:)); B=median(values(:,references),2); valid=isfinite(B) & B>0;
n=sum(valid); trace=nan(size(flat,2),1); basis=nan(1,6);
if n<4, return; end
trace=double(median(values(valid,:)./B(valid),1)).'; scale=median(trace);
if ~isfinite(scale) || scale<=0, return; end
trace=trace/scale; basis=mean(phi(ids(valid),:),1);
end

function yes=protected(f,fmag,df)
yes=fmag>0 && any(abs(f-fmag*(1:5))<=2*df);
end

function rss=residual_energy(f,t,nuisance,y)
A=[nuisance,cos(2*pi*f*t),sin(2*pi*f*t)]; residual=y-A*(A\y); rss=sum(residual.^2);
end

function yes=scalar_number(x)
yes=isnumeric(x) && isreal(x) && isscalar(x) && isfinite(x);
end

function value=opt(s,key,default)
if isfield(s,key), value=s.(key); else, value=default; end
end
