function cols = size_colors(n)
%SIZE_COLORS  为 n 种粒径返回一组区分度高的颜色（n×3，RGB 0~1）
base = [0.85 0.33 0.10;   % 橙
        0.00 0.45 0.74;   % 蓝
        0.47 0.67 0.19;   % 绿
        0.49 0.18 0.56;   % 紫
        0.93 0.69 0.13;   % 黄
        0.30 0.75 0.93];  % 青
cols = base(mod((1:n)-1, size(base,1)) + 1, :);
end
