function gm = gmm_em2d(X, K, w, opts)
%GMM_EM2D  带样本权重的二维全协方差高斯混合模型（EM 算法），不依赖任何工具箱
%   论文在 (d_mes, ∛C) 平面上用全协方差 2D GMM（scikit-learn）区分不同种群（Fig.3）。
%   X：N×2 数据；K：分量数；w：N×1 权重（默认全 1，论文按轨迹长度加权）
%   输出 gm.mu (K×2)、gm.Sigma (2×2×K)、gm.pi (K×1)、gm.resp (N×K 责任度)、
%        gm.idx (N×1 硬分类)、gm.loglik；分量按 mu(:,1)（粒径）升序排列
if nargin < 3 || isempty(w), w = ones(size(X,1), 1); end
if nargin < 4, opts = struct(); end
max_iter = get_opt(opts, 'max_iter', 300);
tol      = get_opt(opts, 'tol', 1e-8);
n_init   = get_opt(opts, 'n_init', 8);
N = size(X, 1); w = w(:)/mean(w(:));

% 标准化（仅用于数值稳定，结果再变换回原坐标）
mu0 = mean(X, 1); sd0 = std(X, 0, 1); sd0(sd0 == 0) = 1;
Z = bsxfun(@rdivide, bsxfun(@minus, X, mu0), sd0);

best = []; best_ll = -inf;
for it = 1:n_init
    % --- k-means++ 初始化 + 若干步 k-means ---
    mu = kmeanspp_init(Z, K, w);
    for q = 1:20
        [~, lab] = min(pairwise_d2(Z, mu), [], 2);
        for k = 1:K
            sel = lab == k;
            if any(sel), mu(k,:) = sum(bsxfun(@times, Z(sel,:), w(sel)), 1)/sum(w(sel)); end
        end
    end
    Sigma = repmat(0.5*eye(2), [1 1 K]); pk = ones(K, 1)/K;
    ll_old = -inf; resp = zeros(N, K);
    for iter = 1:max_iter
        % E 步：责任度
        logp = zeros(N, K);
        for k = 1:K
            logp(:,k) = log(pk(k)) + log_gauss(Z, mu(k,:), Sigma(:,:,k));
        end
        mx  = max(logp, [], 2);
        lse = mx + log(sum(exp(bsxfun(@minus, logp, mx)), 2));
        resp = exp(bsxfun(@minus, logp, lse));
        ll = sum(w.*lse);
        % M 步：加权更新
        for k = 1:K
            rw = resp(:,k).*w; Nk = sum(rw);
            if Nk < 1e-8
                mu(k,:) = Z(randi(N), :); Sigma(:,:,k) = 0.5*eye(2); pk(k) = 1e-3; continue;
            end
            mu(k,:) = sum(bsxfun(@times, Z, rw), 1)/Nk;
            Zc = bsxfun(@minus, Z, mu(k,:));
            Sigma(:,:,k) = (bsxfun(@times, Zc, rw).'*Zc)/Nk + 1e-4*eye(2);   % 微小正则化
            pk(k) = Nk/sum(w);
        end
        pk = pk/sum(pk);
        if abs(ll - ll_old) < tol*abs(ll), break; end
        ll_old = ll;
    end
    if ll > best_ll
        best_ll = ll; best = struct('mu', mu, 'Sigma', Sigma, 'pi', pk, 'resp', resp);
    end
end

% --- 反标准化并按粒径排序 ---
S = diag(sd0);
gm.mu = bsxfun(@plus, bsxfun(@times, best.mu, sd0), mu0);
gm.Sigma = zeros(2, 2, K);
for k = 1:K, gm.Sigma(:,:,k) = S*best.Sigma(:,:,k)*S; end
gm.pi = best.pi; gm.resp = best.resp; gm.loglik = best_ll;
[~, ord] = sort(gm.mu(:,1));
gm.mu = gm.mu(ord,:); gm.Sigma = gm.Sigma(:,:,ord); gm.pi = gm.pi(ord); gm.resp = gm.resp(:,ord);
[~, gm.idx] = max(gm.resp, [], 2);
end

%% ---------- 子函数 ----------
function lp = log_gauss(Z, mu, S)
Zc = bsxfun(@minus, Z, mu);
[R, pf] = chol(S);
if pf, S = S + 1e-6*eye(2); R = chol(S); end
q  = sum((Zc/R).^2, 2);                       % Mahalanobis 距离²
lp = -0.5*q - log(2*pi) - sum(log(diag(R)));
end

function d2 = pairwise_d2(Z, mu)
d2 = bsxfun(@plus, sum(Z.^2, 2), sum(mu.^2, 2).') - 2*Z*mu.';
d2(d2 < 0) = 0;
end

function mu = kmeanspp_init(Z, K, w)
mu = zeros(K, 2);
mu(1,:) = Z(randsample_w(w), :);
for k = 2:K
    d2 = min(pairwise_d2(Z, mu(1:k-1,:)), [], 2);
    mu(k,:) = Z(randsample_w(w.*d2), :);
end
end

function i = randsample_w(pw)
c = cumsum(pw(:));
if c(end) <= 0, i = randi(numel(pw)); return; end
i = find(rand*c(end) <= c, 1);
end

function v = get_opt(s, f, def)
if isfield(s, f), v = s.(f); else, v = def; end
end
