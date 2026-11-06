function tracks = link_trajectories(det_all, p)
%LINK_TRAJECTORIES  扩散尺度归一化的最近邻链接（稀样品，允许短暂丢帧）
%   det_all：cell{F×1}，每帧为 [x y C R] 矩阵（来自 detect_frame）
%   输出 tracks 结构体数组：.frames .x .y .C .n_loc（像素坐标）
%   参数：p.ana.link_r 最大链接位移 [px]；p.ana.link_gap 允许的最大丢帧数
F = numel(det_all);
tracks = struct('frames', {}, 'x', {}, 'y', {}, 'C', {}, 'n_loc', {});
act    = struct('frames', {}, 'x', {}, 'y', {}, 'C', {}, 'last', {});   % 活动轨迹
for f = 1:F
    Dm = det_all{f};
    nd = size(Dm, 1); na = numel(act);
    assigned_d = false(nd, 1); assigned_a = false(na, 1);
    if nd > 0 && na > 0
        ax = arrayfun(@(a) a.x(end), act); ay = arrayfun(@(a) a.y(end), act);
        dist = sqrt(bsxfun(@minus, ax(:), Dm(:,1).').^2 + bsxfun(@minus, ay(:), Dm(:,2).').^2);
        % 布朗位移标准差随间隔平方根增长；丢帧后不能仍用单帧门限。
        elapsed=f-[act.last];
        scaled=dist./sqrt(elapsed(:));
        [dv, idx] = sort(scaled(:));
        for q = 1:numel(dv)                     % 按距离从小到大贪心分配
            if dv(q) > p.ana.link_r, break; end
            [ia, id] = ind2sub([na nd], idx(q));
            if assigned_a(ia) || assigned_d(id), continue; end
            act(ia).frames(end+1) = f;
            act(ia).x(end+1) = Dm(id,1); act(ia).y(end+1) = Dm(id,2); act(ia).C(end+1) = Dm(id,3);
            act(ia).last = f;
            assigned_a(ia) = true; assigned_d(id) = true;
        end
    end
    % 未匹配的检测 → 新轨迹
    for id = find(~assigned_d).'
        act(end+1) = struct('frames', f, 'x', Dm(id,1), 'y', Dm(id,2), 'C', Dm(id,3), 'last', f); %#ok<AGROW>
    end
    % 长时间未更新的活动轨迹 → 关闭
    if ~isempty(act)
        stale = [act.last] < f - p.ana.link_gap;
        for ia = find(stale)
            tracks(end+1) = close_track(act(ia)); %#ok<AGROW>
        end
        act(stale) = [];
    end
end
for ia = 1:numel(act)
    tracks(end+1) = close_track(act(ia)); %#ok<AGROW>
end
end

function t = close_track(a)
t.frames = a.frames(:); t.x = a.x(:); t.y = a.y(:); t.C = a.C(:);
t.n_loc  = numel(a.frames);
end
