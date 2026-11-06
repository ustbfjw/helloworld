function kern = rvt_kernels(N, radii)
%RVT_KERNELS  预计算径向方差变换 (RVT) 所需的环形核及其 FFT
%   RVT（论文参考文献 26，Kashkanova et al., Opt. Express 2021）：
%   对每个像素，计算以其为圆心、半径 r ∈ radii 的各圆环上的像素平均值，
%   再取这些环平均值关于 r 的方差。对称的 iPSF（无论中心亮或暗、是否离焦）在其对称中心给出峰值。
%   radii：半径向量 [px]，例如 [2:8 10 12 14 16]（较大的 r_max 可覆盖离焦粒子的外环）
%   卷积通过 FFT 实现；为避免周期卷积的边缘回绕，图像先镜像填充 pad = max(radii) 像素，
%   因此核尺寸为 (N + 2·pad)。
pad = max(radii);
Np  = N + 2*pad;
[X, Y] = meshgrid(-floor(Np/2):ceil(Np/2)-1);
R  = round(sqrt(X.^2 + Y.^2));
Nr = numel(radii);
Fk = zeros(Np, Np, Nr, 'single');          % 单精度 FFT 更快，精度足够
for i = 1:Nr
    K = double(R == radii(i));
    K = K/sum(K(:));                      % 环平均 → 核归一化
    Fk(:,:,i) = single(fft2(ifftshift(K)));   % 将核中心移到 (1,1) 后做 FFT
end
kern.Fk = Fk; kern.radii = radii; kern.N = N; kern.pad = pad;
kern.idx = [pad:-1:1, 1:N, N:-1:N-pad+1]; % 镜像填充索引
end
