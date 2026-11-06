function metrics=review_fft_replacement(outdir)
% 独立验证test.m的数学核心，不读取用户BMP目录、不执行原脚本。
% 坐标滤波测试只说明单频限制破坏扩散信息，不等于完整iSCAT追踪测试。
if nargin<1, outdir=fileparts(mfilename('fullpath')); end
if ~exist(outdir,'dir'), mkdir(outdir); end
rng(7113,'twister'); Fs=5000; N=5000; D=4.39698; nrep=200;
y=cumsum(sqrt(2*D/Fs)*randn(nrep,N),2);
raw=estimate_D(y,Fs);
freqs=[10 420]; filt=zeros(1,2); retention=zeros(1,2);
for j=1:2
    z=single_bin(y,Fs,freqs(j));
    filt(j)=estimate_D(z,Fs);
    yc=y-mean(y,2);
    retention(j)=sum(z(:).^2)/sum(yc(:).^2);
end
metrics.coordinate_example=struct('D_input',D,'N',N,'nrep',nrep, ...
    'D_raw',raw,'frequencies_Hz',freqs,'D_single_bin',filt, ...
    'D_ratio_to_raw',filt/raw,'centered_signal_energy_retained',retention);
fprintf('COORDINATE ANALOG only: inputD=%.6f rawD=%.6f\n',D,raw);
for j=1:2
    fprintf('  freq=%gHz filteredD=%.9f ratio=%.9g energy=%.9g\n',freqs(j),filt(j),filt(j)/raw,retention(j));
end
% 移动高斯图像（不含iSCAT物理）：演示单频重建的确定性秩2约束。
N=2001; npx=32; t=(0:N-1)/Fs;
x=16+cumsum(sqrt(2*2/Fs)*randn(1,N));
y=16+cumsum(sqrt(2*2/Fs)*randn(1,N));
[X,Y]=meshgrid(1:npx); flat=1+0.5*(X-16)/32+0.3*(Y-16).^2/256;
V=zeros(npx*npx,N);
for k=1:N
    shape=1+0.2*exp(-((X-x(k)).^2+(Y-y(k)).^2)/(2*1.4^2));
    V(:,k)=reshape(flat.*shape,[],1);
end
[R,actual,index]=single_bin(V,Fs,420);
w=2*pi*(index-1)/N;
A=[cos(w*(0:N-1));sin(w*(0:N-1))];
reproj=(R*A')/(A*A')*A;
rank_residual=norm(R-reproj,'fro')/norm(R,'fro');
sv=svd(R(1:17:end,:),'econ');
Vc=V-mean(V,2); retained=sum(R(:).^2)/sum(Vc(:).^2);
metrics.image_example=struct('N',N,'target_Hz',420,'actual_Hz',actual, ...
    'frequency_resolution',Fs/N,'rank2_relative_residual',rank_residual, ...
    'third_to_first_singular_value',sv(3)/sv(1),'temporal_energy_retained',retained);
fprintf('MOVING GAUSSIAN IMAGE: requested420 actual=%.9f deltaf=%.9f rank2Residual=%.3g thirdSingularRatio=%.3g energy=%.9g\n', ...
    actual,Fs/N,rank_residual,sv(3)/sv(1),retained);
assert(rank_residual<1e-12 && sv(3)/sv(1)<1e-12);
% 没有进行空间平场归一化：倍增局部光照会同比放大FFT幅值。
N=5000; tt=(0:N-1)/Fs; harmonic=3*cos(2*pi*420*tt+0.37);
[ra,~,~,amp1]=single_bin(10+harmonic,Fs,420);
[rb,~,~,amp2]=single_bin(2*(10+harmonic),Fs,420);
assert(abs(amp1-3)<1e-12 && abs(amp2/amp1-2)<1e-12);
assert(norm(rb-2*ra)<1e-10);
metrics.flatfield=struct('illumination_ratio',2,'FFT_amplitude_ratio',amp2/amp1);
% 常规FFT正负频率、直流、Nyquist处理核对。
[~,~,~,ampNy]=single_bin(10+2*cos(pi*(0:N-1)),Fs,Fs/2);
[rDC,~,~,ampDC]=single_bin(10+harmonic,Fs,0);
assert(abs(ampNy-2)<1e-12 && abs(ampDC)<1e-12 && norm(rDC)<1e-10);
metrics.boundaries=struct('nyquist_amplitude',ampNy,'DC_after_mean_subtraction',ampDC);
metrics.pass=true;
metrics.scope='Analytic FFT implementation and synthetic counterexamples only; not full iSCAT size validation';
fid=fopen(fullfile(outdir,'fft_replacement_review.json'),'w','n','UTF-8');
fprintf(fid,'%s',jsonencode(metrics,PrettyPrint=true)); fclose(fid);
save(fullfile(outdir,'fft_replacement_review.mat'),'metrics');
fprintf('PASS: mathematical diagnostics; direct replacement is NOT justified.\n');
end

function D=estimate_D(y,Fs)
msd=zeros(4,1);
for lag=1:4
    delta=y(:,lag+1:end)-y(:,1:end-lag);
    msd(lag)=mean(delta(:).^2);
end
c=polyfit((1:4)'/Fs,msd,1); D=c(1)/2;
end

function [R,actual,index,amplitude]=single_bin(V,Fs,target)
% 原test.m第198-210/218-265/404-424行的二维等价数学核心。
N=size(V,2); V=V-mean(V,2); spectrum=fft(V,[],2);
fa=(0:floor(N/2))*Fs/N;
[~,index]=min(abs(fa-target)); actual=fa(index);
isDC=index==1; isNyquist=mod(N,2)==0 && index==N/2+1;
amplitude=abs(spectrum(:,index))/N;
if ~isDC && ~isNyquist, amplitude=2*amplitude; end
filtered=zeros(size(spectrum)); filtered(:,index)=spectrum(:,index);
if ~isDC && ~isNyquist
    other=N-index+2; filtered(:,other)=spectrum(:,other);
end
R=real(ifft(filtered,[],2));
end
