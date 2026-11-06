function [ex, ey] = gauss_ellipse(mu, S, nsig)
%GAUSS_ELLIPSE  二维高斯分布的 nsig 倍标准差等高线（用于绘制 GMM 分量）
[V, E] = eig((S + S.')/2);
E = max(E, 0);
t = linspace(0, 2*pi, 120);
pts = nsig*(V*sqrt(E))*[cos(t); sin(t)];
ex = pts(1,:) + mu(1);
ey = pts(2,:) + mu(2);
end
