function bg = masked_median(V, dets_block, R_mask, strong_thr)
%MASKED_MEDIAN  带掩膜的时间中值背景
%   慢速、强散射的粒子会在普通时间中值中留下"鬼影"（其干涉环在同一像素上停留的帧数比例较高），
%   鬼影结构（对比度可达 ±0.03）会导致弱粒子的误检。本函数利用第一遍检测结果，
%   把强检测周围半径 R_mask 内的像素在对应帧中标记为无效（NaN），再逐像素取剩余帧的中值。
%   V：N×N×Fb 单精度图像块；dets_block：cell{Fb}，每帧 [x y C R]；
%   R_mask：掩膜半径 [px]；strong_thr：仅 RVT 值 ≥ strong_thr 的检测参与掩膜
[N1, N2, Fb] = size(V);
mask = false(N1, N2, Fb);
[dy, dx] = ndgrid(-R_mask:R_mask, -R_mask:R_mask);
disc = dx(:).^2 + dy(:).^2 <= R_mask^2; dx = dx(disc); dy = dy(disc);
for f = 1:Fb
    D = dets_block{f};
    if isempty(D), continue; end
    D = D(D(:,4) >= strong_thr, :);
    for i = 1:size(D, 1)
        cc = round(D(i,1)) + dx; rr = round(D(i,2)) + dy;
        ok = cc >= 1 & cc <= N2 & rr >= 1 & rr <= N1;
        mask(sub2ind([N1 N2 Fb], rr(ok), cc(ok), f*ones(sum(ok), 1))) = true;
    end
end
bg0 = median(V, 3);                              % 无掩膜中值（兜底）
if ~any(mask(:)), bg = bg0; return; end
Vm = V; Vm(mask) = NaN;
S = sort(Vm, 3);                                 % NaN 排在末尾
n = sum(~isnan(Vm), 3);
[ii, jj] = ndgrid(1:N1, 1:N2);
klo = max(floor((n + 1)/2), 1); khi = max(ceil((n + 1)/2), 1);
bg = (S(sub2ind(size(S), ii, jj, klo)) + S(sub2ind(size(S), ii, jj, khi)))/2;
bg(n < 0.3*Fb) = bg0(n < 0.3*Fb);                % 有效帧太少的像素退回普通中值
end
