function [xc, yc] = radial_center(I)
%RADIAL_CENTER  基于径向对称性的亚像素定位（Parthasarathy, Nat. Methods 9, 724 (2012) 算法）
%   输入 I：以候选峰为中心的 (2w+1)×(2w+1) 小窗口
%   输出 (xc, yc)：对称中心相对窗口中心像素的偏移 [px]（xc 沿列/x，yc 沿行/y）
%   原理：径向对称图样中，各点梯度方向线都通过对称中心；
%   以梯度幅值²（并除以到质心的距离）为权重，求到所有梯度线距离平方和最小的点。
[Ny, Nx] = size(I);
I = double(I);
% 半像素网格（梯度定义在相邻 2×2 像素之间）
xm = (-(Nx-1)/2 + 0.5) : ((Nx-1)/2 - 0.5);
ym = ((-(Ny-1)/2 + 0.5) : ((Ny-1)/2 - 0.5)).';
[XM, YM] = meshgrid(xm, ym);
% 沿 45° 对角方向的差分
dIdu = I(1:Ny-1, 2:Nx) - I(2:Ny, 1:Nx-1);
dIdv = I(1:Ny-1, 1:Nx-1) - I(2:Ny, 2:Nx);
h = ones(3)/9;                                  % 3×3 平滑
fdu = conv2(dIdu, h, 'same'); fdv = conv2(dIdv, h, 'same');
dImag2 = fdu.^2 + fdv.^2;                       % 梯度幅值²
% 梯度线的斜率 m 与截距 b（直线 y = m x + b）
m = -(fdv + fdu)./(fdu - fdv);
bad = ~isfinite(m);
if any(bad(:))
    m_un = -(dIdv + dIdu)./(dIdu - dIdv);
    m(bad) = m_un(bad);
    m(~isfinite(m)) = 1e6;                      % 近竖直的梯度线
end
b = YM - m.*XM;
% 权重：梯度幅值² / 到梯度幅值质心的距离
sdI2 = sum(dImag2(:));
if sdI2 <= 0, xc = 0; yc = 0; return; end
xcen = sum(sum(dImag2.*XM))/sdI2;
ycen = sum(sum(dImag2.*YM))/sdI2;
dist = sqrt((XM - xcen).^2 + (YM - ycen).^2);
dist(dist < 0.5) = 0.5;
wgt = dImag2./dist;
% 最小二乘求所有加权梯度线的"交点"
wm2p1 = wgt./(m.^2 + 1);
sw   = sum(wm2p1(:));
smmw = sum(sum(m.^2.*wm2p1));
smw  = sum(sum(m.*wm2p1));
smbw = sum(sum(m.*b.*wm2p1));
sbw  = sum(sum(b.*wm2p1));
dd = smw^2 - smmw*sw;
if abs(dd) < 1e-12, xc = 0; yc = 0; return; end
xc = (smbw*sw  - smw*sbw)/dd;
yc = (smbw*smw - smmw*sbw)/dd;
% 限制在窗口内，防止离群
xc = max(min(xc, (Nx-1)/2), -(Nx-1)/2);
yc = max(min(yc, (Ny-1)/2), -(Ny-1)/2);
end
