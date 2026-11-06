function plot_example_frames(ex, p, tag)
%PLOT_EXAMPLE_FRAMES  展示一帧示例：原始相机图像 / 对比度图（叠加检测结果）/ RVT 图
if nargin < 3, tag = ''; end
sfx = ''; if ~isempty(tag), sfx = ['_' tag]; end
if isempty(ex) || isempty(ex.frame), return; end
N = p.cam.Npx; ax = ((1:N) - (N+1)/2)*p.cam.px*1e6;
fig = figure('Name', ['Example frames ' tag], 'Color', 'w', 'Position', [80 80 1250 400]);
subplot(1,3,1);
imagesc(ax, ax, double(ex.raw)); axis image; colormap(gray); colorbar;
title(sprintf('%s 原始相机图像 (ADU)，第 %d 帧', tag, ex.frame)); xlabel('x (\mum)'); ylabel('y (\mum)');
subplot(1,3,2);
imagesc(ax, ax, ex.C); axis image; colorbar; hold on;
if ~isempty(ex.det)
    plot((ex.det(:,1) - (N+1)/2)*p.cam.px*1e6, (ex.det(:,2) - (N+1)/2)*p.cam.px*1e6, ...
         'o', 'Color', [1 0.2 0.2], 'MarkerSize', 12, 'LineWidth', 1.2);
end
cl = max(abs(ex.C(:))); if cl > 0, caxis([-1 1]*min(cl, 1)); end
title('背景扣除后的对比度图 + 检测位置'); xlabel('x (\mum)'); ylabel('y (\mum)');
subplot(1,3,3);
imagesc(ax, ax, ex.R); axis image; colorbar;
title('径向方差变换 (RVT)'); xlabel('x (\mum)'); ylabel('y (\mum)');
if isfield(p.out, 'dir') && exist(p.out.dir, 'dir')
    print(fig, fullfile(p.out.dir, ['fig2_example_frames' sfx '.png']), '-dpng', '-r150');
end
end
