function dets = detect_frame(C, R, p)
%DETECT_FRAME  单帧粒子检测与两段式亚像素定位
%   C：对比度图；R：RVT 图
%   输出 dets = [x y C_val R_val thr]（每行一个粒子；x = 列坐标，y = 行坐标，1-based，亚像素；thr 为本帧阈值）
%   步骤：
%     1) 鲁棒阈值（median + det_k·MAD，在抽样像素上估计）→ (2w+1)² 邻域局部极大 → 剔除边缘；
%     2) 粗定位：RVT 峰的整数像素位置；
%     3) 细定位：在对比度图 C 上以粗位置为中心取 (2L+1)² 窗口，用径向对称中心算法
%        （radial_center）迭代求 iPSF 的对称中心——iPSF 的环结构关于粒子位置径向对称，
%        这一步对离焦粒子（RVT 峰可能落在环上）尤为重要；
%     4) 合并相距过近的重复检测，剔除强粒子周围的环伪影；在亚像素位置双线性读取对比度（iPSF 中心瓣）。
N = size(R, 1); w = p.ana.det_win; L = p.ana.loc_win;
border = L + 2;                                 % RVT 已做镜像填充，只需为细定位窗口留边

% --- 鲁棒噪声阈值（抽样估计中位数与 MAD） ---
v = R(1:3:end, 1:3:end); v = v(:);
med = median(v); mad = median(abs(v - med));
thr = med + p.ana.det_k*1.4826*mad;

% --- 可分离的局部极大滤波 ---
M1 = R;
for dx = -w:w, M1 = max(M1, circshift(R, [0 dx])); end
M = M1;
for dy = -w:w, M = max(M, circshift(M1, [dy 0])); end
pk = (R >= M) & (R > thr);
pk([1:border, N-border+1:N], :) = false;
pk(:, [1:border, N-border+1:N]) = false;
[rr, cc] = find(pk);
if isempty(rr), dets = zeros(0, 5); return; end
vals = R(sub2ind(size(R), rr, cc));
[vals, ord] = sort(vals, 'descend'); rr = rr(ord); cc = cc(ord);

% --- 粗位置预筛：剔除强候选周围的弱环伪峰（减少细定位的计算量） ---
keep = true(size(rr));
for i = 2:numel(rr)
    d2 = (rr(1:i-1) - rr(i)).^2 + (cc(1:i-1) - cc(i)).^2;
    if any(keep(1:i-1) & d2 < p.ana.ring_excl^2 & vals(i) < p.ana.ring_ratio*vals(1:i-1)), keep(i) = false; end
end
rr = rr(keep); cc = cc(keep); vals = vals(keep);

% --- 细定位 + 去重 ---
Cd = double(C);
dets = zeros(0, 5);
for i = 1:numel(rr)
    xc = cc(i); yc = rr(i); step = inf;
    for it = 1:p.ana.loc_iter
        ci = round(xc); ri = round(yc);
        if ci < L+1 || ci > N-L || ri < L+1 || ri > N-L, step = inf; break; end
        [dx, dy] = radial_center(Cd(ri-L:ri+L, ci-L:ci+L));
        step = hypot(ci + dx - xc, ri + dy - yc);
        xc = ci + dx; yc = ri + dy;
    end
    % 收敛性与边界检查：中心位于视场外的粒子，其伸入视场的环会把细定位"卡"在边缘附近，
    % 表现为未收敛或贴边，此类检测不可靠，予以剔除
    if ~(step < 1) || xc < L+2 || xc > N-L-1 || yc < L+2 || yc > N-L-1, continue; end
    if ~isempty(dets)
        d2 = (dets(:,1) - xc).^2 + (dets(:,2) - yc).^2;
        if any(d2 < p.ana.min_sep^2), continue; end                  % 与更强的检测重复（同一粒子）
        % 强散射粒子的离焦环会在其周围产生弱得多的 RVT 次峰：落在更强检测的 ring_excl 范围内、
        % 且 RVT 远小于该检测者视为环伪影剔除（强度相当的则保留，视为另一粒子）
        if any(d2 < p.ana.ring_excl^2 & vals(i) < p.ana.ring_ratio*dets(:,4)), continue; end
    end
    ci = floor(xc); ri = floor(yc); fx = xc - ci; fy = yc - ri;          % 局部双线性插值
    cval = (1-fy)*((1-fx)*Cd(ri,ci) + fx*Cd(ri,ci+1)) + fy*((1-fx)*Cd(ri+1,ci) + fx*Cd(ri+1,ci+1));
    dets(end+1, :) = [xc, yc, cval, vals(i), thr]; %#ok<AGROW>
end
end
