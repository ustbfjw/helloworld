function set_seed(seed)
%SET_SEED  设置随机数种子（兼容 MATLAB 与 GNU Octave）
if exist('rng', 'file') || exist('rng', 'builtin')
    rng(seed);
else
    rand('state', seed); randn('state', seed);
end
end
