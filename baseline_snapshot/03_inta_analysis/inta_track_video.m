function [tracks, ex] = inta_track_video(video, p)
%INTA_TRACK_VIDEO  对一段 iSCAT 视频做 iNTA 前处理与追踪：
%   时间中值背景扣除 → 对比度图 C = (I - I_bg)/I_bg → RVT → 检测/两段式亚像素定位 → 轨迹链接
%   （对应论文 Methods "Code availability"：median background correction → RVT → tracking）
%   若第一遍检测发现强散射粒子（RVT ≫ 阈值），则用带掩膜的中值重算背景并再检测一遍，
%   以消除慢速强粒子在中值背景中留下的鬼影（p.ana.bg_masked = true 时启用）。
%   为节省内存，按时间块处理，不保存整段对比度/RVT 视频；返回最长轨迹中间一帧作为示例 ex。
N = size(video, 1); F = size(video, 3);
mode='median';
if isfield(p.ana,'background_mode'), mode=p.ana.background_mode; end
if strcmpi(mode,'spectral_common')
    [gain,illumination_info]=estimate_common_illumination(video,p);
elseif strcmpi(mode,'median')
    gain=ones(F,1); illumination_info=struct('mode','median','applied',false,'gain',gain);
else
    error('background_mode必须是median或spectral_common。');
end
kern    = rvt_kernels(N, p.ana.rvt_radii);
det_all = cell(F, 1);
ex = struct('frame', [], 'raw', [], 'C', [], 'R', [], 'det', []);
bg_all = zeros(N, N, 0, 'single'); blk_start = [];
for b0 = 1:p.ana.bg_block:F
    idx = b0:min(b0 + p.ana.bg_block - 1, F);
    V   = single(video(:,:,idx));
    % 只除去估计的共同乘性照明；帧数、时间戳和宽带颗粒信号不作时间投影。
    if any(gain(idx)~=1), V=V./reshape(single(gain(idx)),1,1,[]); end
    bg  = median(V, 3);                                   % 第一遍：普通时间中值背景
    dets_b = detect_block(V, bg, kern, p);
    if p.ana.bg_masked
        % 是否存在强散射粒子：任一检测的 RVT 超过 strong_factor × 本块典型阈值
        rv = cellfun(@(d) max([d(:,4); 0]), dets_b);
        nonempty = dets_b(~cellfun(@isempty, dets_b));
        thr_typ = [];
        if ~isempty(nonempty), thr_typ = median(cellfun(@(d) local_thr(d), nonempty)); end
        if ~isempty(thr_typ) && any(rv > p.ana.strong_factor*thr_typ)
            bg = masked_median(V, dets_b, p.ana.mask_r, p.ana.strong_factor*thr_typ);
            dets_b = detect_block(V, bg, kern, p);               % 第二遍：掩膜背景
        end
    end
    det_all(idx) = dets_b;
    bg_all(:,:,end+1) = bg; blk_start(end+1) = b0; %#ok<AGROW>
end
tracks = link_trajectories(det_all, p);

% ---- 示例帧：取最长轨迹的中间一帧，重算其对比度图与 RVT ----
if ~isempty(tracks)
    [~, k] = max([tracks.n_loc]);
    f_ex = tracks(k).frames(ceil(end/2));
    ib = find(blk_start <= f_ex, 1, 'last');
    bg = bg_all(:,:,ib);
    Vx=single(video(:,:,f_ex));
    if gain(f_ex)~=1, Vx=Vx/single(gain(f_ex)); end
    Cx = (Vx - bg)./bg;
    ex = struct('frame', f_ex, 'raw', video(:,:,f_ex), 'C', Cx, 'R', rvt_apply(Cx, kern), 'det', det_all{f_ex});
end
ex.illumination=illumination_info;
end

%% ---------- 子函数 ----------
function dets_b = detect_block(V, bg, kern, p)
% 对一个时间块：背景扣除 → 对比度 → RVT（分批）→ 逐帧检测
batch = 100; Fb = size(V, 3);
Cb = bsxfun(@rdivide, bsxfun(@minus, V, bg), bg);
dets_b = cell(Fb, 1);
for s0 = 1:batch:Fb
    sidx = s0:min(s0 + batch - 1, Fb);
    Rb = rvt_apply(Cb(:,:,sidx), kern);
    for q = 1:numel(sidx)
        dets_b{sidx(q)} = detect_frame(Cb(:,:,sidx(q)), Rb(:,:,q), p);
    end
end
end

function t = local_thr(d)
% detect_frame 把本帧阈值放在第 5 列（若有）；否则用最弱检测值近似
if size(d, 2) >= 5, t = d(1, 5); else, t = min(d(:, 4)); end
end
