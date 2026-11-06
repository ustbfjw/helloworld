function plot_trajectories_3d(traj, p, tag)
%PLOT_TRAJECTORIES_3D  三维轨迹图（含盖玻片/焦平面/视场示意）+ 沿梯度方向位移 + 场波形
if nargin < 3, tag = ''; end
sfx = ''; if ~isempty(tag), sfx = ['_' tag]; end
cols = size_colors(numel(p.part.d_list));
nd   = numel(p.part.d_list);
fig  = figure('Name', ['3D trajectories ' tag], 'Color', 'w', 'Position', [40 40 1300 560]);

%% (a) 三维轨迹（未折叠坐标，单位 µm）
subplot(2,3,[1 2 4 5]); hold on;
Np = numel(traj.d);
for j = 1:Np
    P = squeeze(traj.pos(:,j,:))*1e6;
    c = cols(traj.size_idx(j), :);
    plot3(P(:,1), P(:,2), P(:,3), '-', 'Color', c, 'LineWidth', 0.7);
    plot3(P(1,1), P(1,2), P(1,3), 'o', 'MarkerFaceColor', c, ...
          'MarkerEdgeColor', 'k', 'MarkerSize', 4);
end
% 盖玻片 (z=0) 与焦平面 (z=z_f)，按 7×7 µm 视场大小绘制
fov = p.cam.Npx*p.cam.px*1e6/2;
[Xs, Ys] = meshgrid([-fov fov], [-fov fov]);
surf(Xs, Ys, zeros(2), 'FaceColor', [0.6 0.6 0.6], 'FaceAlpha', 0.35, 'EdgeColor', 'k');
surf(Xs, Ys, p.opt.z_f*1e6*ones(2), 'FaceColor', [0.2 0.6 1], 'FaceAlpha', 0.25, 'EdgeColor', 'b');
xlabel('x (\mum)'); ylabel('y (\mum)'); zlabel('z (\mum)');
grid on; box on; view(38, 24); axis equal;
title(sprintf('三维轨迹：布朗运动 + 交变磁泳  (B_0 = %.0f mT, G_0 = %.0f T/m, f = %.0f Hz, %s)', ...
      p.mag.B0*1e3, p.mag.G0, p.mag.f, p.mag.wave));
h = []; lbl = cell(1, nd); lh = {};
for k = 1:nd
    lbl{k} = sprintf('d = %g nm', p.part.d_list(k)*1e9);
    if any(traj.size_idx == k)
        h(end+1) = plot3(nan, nan, nan, '-', 'Color', cols(k,:), 'LineWidth', 2); %#ok<AGROW>
        lh{end+1} = lbl{k}; %#ok<AGROW>
    end
end
if ~isempty(h), legend(h, lh, 'Location', 'northeast'); end

%% (b) 沿梯度方向 e_g 的位移（每种粒径取第一个粒子）
subplot(2,3,3); hold on;
hb = []; lb = {};
for k = 1:nd
    j = find(traj.size_idx == k, 1);
    if isempty(j), continue; end                         % 该粒径在本样品中不存在
    xg = reshape(traj.pos(:,j,:), [], 3)*traj.e_g.';     % 投影到 e_g [m]
    hb(end+1) = plot(traj.t*1e3, (xg - xg(1))*1e6, '-', 'Color', cols(k,:)); %#ok<AGROW>
    lb{end+1} = lbl{k}; %#ok<AGROW>
end
xlabel('t (ms)'); ylabel('沿 e_g 位移 (\mum)'); grid on; box on;
title('沿梯度方向的位移（布朗 + 往复磁泳）');
if ~isempty(hb), legend(hb, lb, 'Location', 'best'); end

%% (c) 场波形 s(t)
subplot(2,3,6);
plot(traj.t*1e3, traj.s_step, 'k-', 'LineWidth', 1.2);
ylim([-1.3 1.3]); xlabel('t (ms)'); ylabel('s(t) = G(t)/G_0'); grid on; box on;
title(sprintf('交变梯度波形（%s, %.0f Hz）', p.mag.wave, p.mag.f));

if isfield(p.out, 'dir') && exist(p.out.dir, 'dir')
    print(fig, fullfile(p.out.dir, ['fig1_trajectories_3d' sfx '.png']), '-dpng', '-r150');
end
end
