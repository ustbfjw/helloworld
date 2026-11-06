function R = rvt_apply(C, kern)
%RVT_APPLY  对图像（或图像堆栈）做径向方差变换：RVT = Var_r[ 环平均_r ]
%   C: N×N×F（对比度图）。先镜像填充 pad 像素，再用 FFT 做环平均卷积，最后裁回 N×N。
Nr = size(kern.Fk, 3); N = kern.N; pad = kern.pad;
Cp = single(C(kern.idx, kern.idx, :));                % 镜像填充（单精度）
FC = fft(fft(Cp, [], 1), [], 2);
S1 = zeros(size(Cp), 'single'); S2 = zeros(size(Cp), 'single');
for i = 1:Nr
    M  = real(ifft(ifft(bsxfun(@times, FC, kern.Fk(:,:,i)), [], 1), [], 2));
    S1 = S1 + M;
    S2 = S2 + M.^2;
end
R = S2/Nr - (S1/Nr).^2;
R = double(R(pad+1:pad+N, pad+1:pad+N, :));
R(R < 0) = 0;
end
