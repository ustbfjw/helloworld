function q = field_integral(t0, t1, p)
%FIELD_INTEGRAL  精确计算归一化梯度波形积分 q=积分 s(t)dt，单位秒。
%   t0/t1 支持标量、向量和 MATLAB 隐式扩展；跨方波反转时也适用。
%   生成运动与反演共用此函数，避免把跨反转位移误当作恒速位移。
omega = 2*pi*p.mag.f;
if omega == 0
    q = field_gradient(0,p).*(t1-t0);
    return;
end
a = omega*t0+p.mag.phase0; b = omega*t1+p.mag.phase0;
switch lower(p.mag.wave)
    case 'sine'
        q = (cos(a)-cos(b))/omega;
    case 'square'
        % sign(sin(phi)) 的周期原函数；每个完整周期的积分严格为零。
        A = pi-abs(mod(a,2*pi)-pi);
        B = pi-abs(mod(b,2*pi)-pi);
        q = (B-A)/omega;
    otherwise
        error('未知波形类型: %s（可选 square | sine）',p.mag.wave);
end
end
