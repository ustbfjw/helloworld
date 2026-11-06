function metrics=test_targeted_ablations(outdir,reference_root)
%TEST_TARGETED_ABLATIONS 三组件开关的定向合成回归，不替代完整粒径验收。
% 六种照明条件只由测试构造；真增益只在估计完成后用于误差评价。
% reference_root可选：与旧交付工程逐条件比较缺省路径的数值结果。
if nargin<1, outdir=''; end
if nargin<2, reference_root=''; end
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root); setup_paths(root);
p=sim_params('full'); p.sim.fps=5000; p.mag.f=10;
N=32; F=5000; t=(0:F-1)'/p.sim.fps;
[X,Y]=meshgrid(linspace(-1,1,N));
B=800*(0.65+0.35*exp(-2*(X.^2+Y.^2))).*(1+0.03*cos(7*X).*sin(5*Y));
rs=RandStream('mt19937ar','Seed',20261105);
clean=repmat(single(B),1,1,F)+single(0.2*randn(rs,N,N,F));
names={'static','global_423p7','near_11p4','near_13p4','same_10','local_left_423p7'};
frequencies=[423.7 423.7 11.4 13.4 10 423.7];
methods={'full','no_spatial_gate','no_magnetic_guard','no_refinement'};
flags={'enable_spatial_gate','enable_magnetic_guard','refine_frequency'};
% 行为断言只针对本测试的已知合成输入，不预设真实iSCAT视频结果。
expected_applied=logical([0 0 0 0; 1 1 1 1; 0 0 1 0; ...
    1 1 1 1; 0 0 1 0; 0 1 0 0]);
rows={}; defaults_identical=true; reference_identical=true;
for k=1:numel(names)
    amplitude=0.08; if k==1, amplitude=0; end
    offsets=((1:p.sim.n_sub)-0.5)*p.cam.exposure/p.sim.n_sub;
    left=1+amplitude*mean(sin(2*pi*frequencies(k)*(t+offsets)+0.37),2);
    right=left; if k==6, right=ones(F,1); end
    video=clean;
    video(:,1:N/2,:)=clean(:,1:N/2,:).*reshape(single(left),1,1,F);
    video(:,N/2+1:end,:)=clean(:,N/2+1:end,:).*reshape(single(right),1,1,F);
    before=rng; [default_gain,default_info]=estimate_common_illumination(video,p);
    pp=p;
    for j=1:numel(flags), pp.ana.illum.(flags{j})=true; end
    [explicit_gain,explicit_info]=estimate_common_illumination(video,pp);
    assert(isequal(default_gain,explicit_gain) && isequaln(default_info,explicit_info), ...
        '缺省设置与显式全开不完全一致。');
    defaults_identical=defaults_identical && isequal(default_gain,explicit_gain);
    if ~isempty(reference_root)
        [reference_gain,reference_info]=reference_estimate(video,p,reference_root);
        assert(isequal(default_gain,reference_gain),'缺省增益与旧交付工程不完全一致。');
        numeric_fields={'applied','accepted_frequencies','amplitudes','coherence', ...
            'amplitude_cv','tile_traces','gain','raw_common','frequency_axis', ...
            'spectrum','rejected_frequencies','frequency_resolution'};
        for j=1:numel(numeric_fields)
            key=numeric_fields{j};
            assert(isequaln(default_info.(key),reference_info.(key)), ...
                '缺省诊断数值与旧交付工程不完全一致：%s。',key);
        end
        reference_identical=reference_identical && isequal(default_gain,reference_gain);
    end
    gains=cell(1,numel(methods)); infos=cell(1,numel(methods));
    for m=1:numel(methods)
        pm=pp;
        if m>1, pm.ana.illum.(flags{m-1})=false; end
        if m==1
            gain=default_gain; info=default_info;
        else
            [gain,info]=estimate_common_illumination(video,pm);
        end
        assert(info.applied==expected_applied(k,m), ...
            '合成条件%s下方法%s的应用状态不符合组件定义。',names{k},methods{m});
        for j=1:numel(flags)
            assert(info.(flags{j})==pm.ana.illum.(flags{j}),'info未准确记录开关。');
        end
        if ~info.applied, assert(isequal(gain,ones(F,1))); end
        if m==4 && info.applied
            grid=info.accepted_frequencies/info.frequency_resolution;
            assert(all(grid==round(grid)),'关闭细化后频率必须落在DFT网格。');
        end
        [all_rms,left_rms,right_rms]=gain_residual(gain,left,right);
        row=struct('condition',names{k},'method',methods{m},'applied',info.applied, ...
            'accepted_frequencies_Hz',info.accepted_frequencies, ...
            'gain_residual_rms',all_rms,'left_residual_rms',left_rms, ...
            'right_residual_rms',right_rms,'coherence',info.coherence, ...
            'amplitude_cv',info.amplitude_cv,'rejection_reasons',{info.rejection_reasons});
        rows{end+1}=row; %#ok<AGROW>
        gains{m}=gain; infos{m}=info;
        fprintf('TARGETED %s / %s applied=%d gainRMS=%.8g f=%s\n', ...
            names{k},methods{m},info.applied,all_rms,mat2str(info.accepted_frequencies,8));
    end
    assert(isequal(before,rng),'开关测试改变了仿真全局随机流。');
    if k==2 || k==4
        refined=gain_residual(gains{1},left,right);
        coarse=gain_residual(gains{4},left,right);
        assert(refined<3e-4 && refined<0.02*coarse,'连续频率细化未恢复非整频照明。');
        assert(abs(infos{1}.accepted_frequencies(1)-frequencies(k))<0.02);
        assert(isequal(gains{1},gains{2}) && isequal(gains{1},gains{3}), ...
            '无对应拒绝条件时，关闭空间或磁频门控不应改变估计。');
    elseif k==3 || k==5
        assert(any(contains(default_info.rejection_reasons,'磁驱动')));
        assert(gain_residual(gains{3},left,right)<3e-4, ...
            '关闭磁频保护后应恢复本测试真实注入的同步公共照明。');
    elseif k==6
        assert(any(contains(default_info.rejection_reasons,'一致性')));
        assert(any(infos{2}.amplitude_cv>p.ana.illum.max_amplitude_cv), ...
            '关闭空间门控后仍需保留幅度不一致诊断。');
        [~,~,quiet_rms]=gain_residual(gains{2},left,right);
        assert(quiet_rms>0.01,'局部照明误校正应在未受扰半区显现，不能被空间平均掩盖。');
    end
    % 仅改变生成器真值字段，估计结果必须完全不变。
    pt=p; pt.cam.flicker_frequency=999; pt.cam.flicker_amplitude=0.9; pt.cam.flicker_phase=1.2;
    [other,~]=estimate_common_illumination(video,pt);
    assert(isequal(default_gain,other),'估计器读取了照明生成器真值字段。');
end

% 单独验证相位一致性消融；幅度一致不能替代相位一致。
wave=reshape(single(0.08*sin(2*pi*423.7*t+0.37)),1,1,F);
video=clean.*(1+wave);
video(1:N/2,1:N/2,:)=clean(1:N/2,1:N/2,:).*(1-wave);
[~,full_info]=estimate_common_illumination(video,p);
pp=p; pp.ana.illum.enable_spatial_gate=false;
[~,off_info]=estimate_common_illumination(video,pp);
assert(~full_info.applied && any(contains(full_info.rejection_reasons,'相位一致性')));
assert(off_info.applied && any(off_info.coherence<p.ana.illum.min_coherence));

% 粗格点13Hz在保护区外，连续细化到12.55Hz后进入10.6±2Hz保护区。
% 此辅助输入专门覆盖第二次磁频复查，不改变六条件或算法默认阈值。
pp=p; pp.mag.f=10.6;
video=clean.*reshape(single(1+0.08*sin(2*pi*12.55*t+0.37)),1,1,F);
[g,guard_info]=estimate_common_illumination(video,pp);
assert(isequal(g,ones(F,1)) && ~guard_info.applied);
assert(any(contains(guard_info.rejection_reasons,'频率复查')));
pp.ana.illum.enable_magnetic_guard=false;
[~,unguarded_info]=estimate_common_illumination(video,pp);
assert(unguarded_info.applied && abs(unguarded_info.accepted_frequencies(1)-12.55)<0.02);

% 三个相位相差120度的象限与一个静态象限：空间中位数可生成三次谐波，
% 但各象限在该频率没有足够振幅。即使关闭一致性门控也必须拒绝此伪分量。
video=clean; phases=[0 2*pi/3 4*pi/3]; q=0;
for iy=1:2
    rr=(iy-1)*N/2+(1:N/2);
    for ix=1:2
        q=q+1; cc=(ix-1)*N/2+(1:N/2);
        if q<=3
            video(rr,cc,:)=clean(rr,cc,:).* ...
                reshape(single(1+0.08*sin(2*pi*420*t+phases(q))),1,1,F);
        end
    end
end
pp=p; pp.ana.illum.enable_spatial_gate=false;
[g,amplitude_info]=estimate_common_illumination(video,pp);
assert(isequal(g,ones(F,1)) && ~amplitude_info.applied);
assert(any(contains(amplitude_info.rejection_reasons,'共同振幅不足')), ...
    '关闭一致性门控不能绕过象限共同振幅下限。');

% 关闭组件不能关闭最小幅度或正增益范围约束。
for j=1:numel(flags), pp.ana.illum.(flags{j})=false; end
[g,small_info]=estimate_common_illumination(clean.* ...
    reshape(single(1+0.0005*sin(2*pi*420*t)),1,1,F),pp);
assert(isequal(g,ones(F,1)) && ~small_info.applied);
[g,large_info]=estimate_common_illumination(clean.* ...
    reshape(single(1+0.6*sin(2*pi*420*t)),1,1,F),pp);
assert(isequal(g,ones(F,1)) && ~large_info.applied);
assert(contains(large_info.rejected_reason,'安全正值范围'));

% 严格标量logical；数值0/1不隐式转换，避免批处理配置歧义。
bad_values={[],[true false],0,1,NaN,'false',"true"}; invalid_count=0;
for j=1:numel(flags)
    for k=1:numel(bad_values)
        pp=p; pp.ana.illum.(flags{j})=bad_values{k}; failed=false;
        try
            estimate_common_illumination(ones(8,8,16,'single'),pp);
        catch err
            assert(strcmp(err.identifier,'iNTA:illum:InvalidBooleanFlag'));
            failed=true;
        end
        assert(failed,'非法开关值未被拒绝。'); invalid_count=invalid_count+1;
    end
end
metrics=struct('scope','Synthetic ablation semantics, not end-to-end sizing acceptance', ...
    'matlab_version',version,'seed',20261105,'n_frames',F,'fps',p.sim.fps, ...
    'default_equals_explicit_true',defaults_identical, ...
    'reference_compared',~isempty(reference_root),'reference_root',reference_root, ...
    'reference_default_identical',reference_identical && ~isempty(reference_root), ...
    'records',[rows{:}],'phase_gate_ablation',true,'refined_magnetic_guard',true, ...
    'median_harmonic_amplitude_guard',true,'amplitude_and_gain_constraints',true, ...
    'invalid_flag_values_rejected',invalid_count,'no_truth_leakage',true,'rng_unchanged',true,'pass',true);
if ~isempty(outdir)
    if ~exist(outdir,'dir'), mkdir(outdir); end
    fid=fopen(fullfile(outdir,'targeted_ablations_tests.json'),'w','n','UTF-8');
    assert(fid>=0,'无法创建测试报告。'); cleanup=onCleanup(@() fclose(fid));
    fprintf(fid,'%s',jsonencode(metrics,PrettyPrint=true));
end
fprintf('PASS test_targeted_ablations: six conditions, component switches, defaults, flags and constraints.\n');
end

function [rms_all,rms_left,rms_right]=gain_residual(gain,left,right)
% 两个等面积区分别按时间均值归一化，最后平均平方；不先平均左右区。
left=left./gain; left=left/mean(left)-1;
right=right./gain; right=right/mean(right)-1;
rms_left=sqrt(mean(left.^2)); rms_right=sqrt(mean(right.^2));
rms_all=sqrt((rms_left^2+rms_right^2)/2);
end

function [gain,info]=reference_estimate(video,p,reference_root)
% 可选旧版数值对照；立即恢复路径，旧版不需要具有新开关字段。
reference_dir=fullfile(reference_root,'03_inta_analysis');
assert(exist(fullfile(reference_dir,'estimate_common_illumination.m'),'file')==2);
old_path=path; cleanup=onCleanup(@() restore_path(old_path));
addpath(reference_dir,'-begin'); clear estimate_common_illumination
assert(strcmp(fileparts(which('estimate_common_illumination')),reference_dir));
[gain,info]=estimate_common_illumination(video,p);
end

function restore_path(old_path)
path(old_path); clear estimate_common_illumination
end
