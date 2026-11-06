function plot_ipsf_examples(tab, p)
%PLOT_IPSF_EXAMPLES  不同离焦量下的 iPSF 对比度图及其 RVT（论文 Fig.1b 风格）
us   = [-0.6 -0.2 0 0.2 0.6]*1e-6;     % 离焦量 u = z - z_f
eps0 = 0.1;                            % 示例散射幅值比
W    = 22; px = p.cam.px;
[X, Y] = meshgrid((-W:W)*px); rho = sqrt(X.^2 + Y.^2);
kern = rvt_kernels(2*W+1, p.ana.rvt_radii);
fig = figure('Name', 'iPSF examples', 'Color', 'w', 'Position', [60 60 1250 480]);
for i = 1:numel(us)
    H = ipsf_eval(tab, rho, us(i));
    C = abs(1 + eps0*H).^2 - 1;        % 对比度图 (I - I_bg)/I_bg
    R = rvt_apply(C, kern);
    subplot(2, numel(us), i);
    imagesc((-W:W)*px*1e6, (-W:W)*px*1e6, C); axis image; colormap(gray);
    caxis([-1 1]*max(abs(C(:))));
    title(sprintf('iPSF, u = %+.1f \\mum', us(i)*1e6)); xlabel('x (\mum)'); ylabel('y (\mum)');
    subplot(2, numel(us), numel(us)+i);
    imagesc((-W:W)*px*1e6, (-W:W)*px*1e6, R); axis image;
    title('RVT'); xlabel('x (\mum)'); ylabel('y (\mum)');
end
if isfield(p.out, 'dir') && exist(p.out.dir, 'dir')
    print(fig, fullfile(p.out.dir, 'fig0_ipsf_examples.png'), '-dpng', '-r150');
end
end
