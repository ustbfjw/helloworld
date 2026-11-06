function [s, G] = field_gradient(t, p)
%FIELD_GRADIENT  交变梯度磁场的时间波形
%   s(t) ∈ [-1, 1] 为归一化波形；G(t) = G0·s(t) 为瞬时梯度 [T/m]
%   t 可为标量或向量
ph = 2*pi*p.mag.f*t + p.mag.phase0;
switch lower(p.mag.wave)
    case 'square'
        s = sign(sin(ph));
        s(s == 0) = 1;
    case 'sine'
        s = sin(ph);
    otherwise
        error('未知波形类型: %s（可选 square | sine）', p.mag.wave);
end
G = p.mag.G0 * s;
end
