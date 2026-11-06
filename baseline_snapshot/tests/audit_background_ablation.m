function results=audit_background_ablation(project_root,outdir)
%AUDIT_BACKGROUND_ABLATION 小样本配对消融，非完整精度验收。
% 同一轨迹、同一渲染前RNG状态：clean与8%/423.7Hz曝光内公共照明变化。
% 光强改变散粒噪声幅度，共用随机数不意味着两视频噪声数值完全相同。
% 不用真位置/真D筛选轨迹；真gain与真坐标仅在各方法分析后评价。
if nargin<1
    parent=fileparts(fileparts(mfilename('fullpath')));
    if exist(fullfile(parent,'setup_paths.m'),'file')
        project_root=parent; % 放在交付工程tests目录时可直接运行
    else
        project_root=fullfile(parent,'iNTA_magnetic_fused'); % 本次工作目录布局
    end
end
if nargin<2, outdir=fullfile(fileparts(mfilename('fullpath')),'ablation_output'); end
addpath(project_root); setup_paths(project_root);
if ~exist(outdir,'dir'), mkdir(outdir); end
diary(fullfile(outdir,'ablation_run.log')); cleanup=onCleanup(@() diary('off')); %#ok<NASGU>
names={'mono30','mono100','mono200'};
methods={'clean_median','clean_spectral_common','stress_median','stress_PN','stress_spectral_common'};
rows={}; detail_rows={}; all_results=struct(); counterexamples=struct();
base=sim_params('quick'); tab=ipsf_lookup_table(base);
fprintf('Small-sample development ablation; seed=1, 3 videos x 1500 frames per group.\n');
fprintf('No 5%% size accuracy assertion; PN=full-frame summed-intensity normalization.\n');
for k=1:3
    p=apply_scenario(base,names{k}); p.sim.seed=1;
    p.out.make_plots=false; p.out.save_video=false; p.out.save_diagnostics=false;
    pools=repmat({empty_tracks()},1,numel(methods));
    gain_squared=zeros(1,numel(methods)); gain_count=zeros(1,numel(methods)); applied=zeros(1,numel(methods));
    clean_same=true; fft_done=false;
    for v=1:p.sim.n_videos
        fprintf('\n%s video %d/%d\n',names{k},v,p.sim.n_videos);
        traj=simulate_trajectories(p,p.sim.seed+100*(v-1)); state=rng;
        pc=p; pc.cam.flicker_amplitude=0; pc.cam.flicker_frequency=423.7; pc.cam.flicker_phase=0.37;
        rng(state); [clean,meta_clean]=render_iscat_video(traj,pc,tab);
        ps=pc; ps.cam.flicker_amplitude=0.08;
        rng(state); [stress,meta_stress]=render_iscat_video(traj,ps,tab);
        assert(meta_clean.n_sat==0 && meta_stress.n_sat==0,'Ablation requires unsaturated videos.');
        gains=cell(1,numel(methods)); infos=cell(1,numel(methods)); tracksets=cell(1,numel(methods));
        for m=1:numel(methods)
            pm=p; pm.ana.background_mode='median';
            if m<=2, video=clean; truth_gain=meta_clean.illumination_gain;
            else, video=stress; truth_gain=meta_stress.illumination_gain; end
            if m==2 || m==5
                pm.ana.background_mode='spectral_common';
                [tr,ex]=inta_track_video(video,pm);
                gain=ex.illumination.gain; info=ex.illumination;
            elseif m==4
                % 经典功率归一化基线：不做频谱检测或频率选择。
                gain=double(squeeze(sum(sum(single(video),1),2))); gain=gain/median(gain);
                normalized=single(video)./reshape(single(gain),1,1,[]);
                [tr,ex]=inta_track_video(normalized,pm); %#ok<ASGLU>
                info=struct('applied',true,'accepted_frequencies',[]);
            else
                [tr,ex]=inta_track_video(video,pm); %#ok<ASGLU>
                gain=ones(p.sim.n_frames,1); info=struct('applied',false,'accepted_frequencies',[]);
            end
            tr=analyze_tracks(tr,traj,pm,v);
            tracksets{m}=tr; pools{m}=[pools{m},tr]; %#ok<AGROW>
            gains{m}=gain; infos{m}=info;
            residual=relative_gain_error(gain,truth_gain);
            gain_squared(m)=gain_squared(m)+sum(residual.^2); gain_count(m)=gain_count(m)+numel(residual);
            applied(m)=applied(m)+info.applied;
            rs=inta_population_analysis(tr,pm);
            d=make_row(names{k},methods{m},rs,tr,p,sqrt(mean(residual.^2)),double(info.applied));
            d.video=v; d.seed=p.sim.seed+100*(v-1); d.accepted_frequencies_Hz=info.accepted_frequencies;
            detail_rows{end+1}=d; %#ok<AGROW>
            fprintf('%s: D=%g d=%g tracks=%d loc=%d gainRMS=%.6g applied=%d\n', ...
                methods{m},d.D_um2_s,d.diameter_nm,d.n_tracks,d.n_loc,d.gain_residual_rms,d.applied_videos);
        end
        clean_same=clean_same && isequaln(tracksets{1},tracksets{2});
        % 反例仅在实际有长度合格检测轨迹的一段压力视频运行，不筛真粒子ID。
        if ~fft_done && any([tracksets{3}.n_loc]>=p.ana.min_len)
            [tr_fft,fftinfo]=single_frequency_counterexample(stress,p,420);
            tr_fft=analyze_tracks(tr_fft,traj,p,v); rs_fft=inta_population_analysis(tr_fft,p);
            row=make_row(names{k},'stress_single_frequency_420_diagnostic',rs_fft,tr_fft,p,nan,1);
            row.video=v; row.seed=p.sim.seed+100*(v-1); row.accepted_frequencies_Hz=fftinfo.actual_frequency_Hz;
            detail_rows{end+1}=row; %#ok<AGROW>
            counterexamples.(names{k})=struct('record',row,'FFT',fftinfo);
            fprintf('SINGLE-FREQ DIAGNOSTIC: actual=%gHz D=%g d=%g tracks=%d loc=%d energy=%g\n', ...
                fftinfo.actual_frequency_Hz,row.D_um2_s,row.diameter_nm,row.n_tracks,row.n_loc,fftinfo.energy_fraction);
            fft_done=true;
        end
        clear clean stress video normalized
    end
    group=struct();
    for m=1:numel(methods)
        res=inta_population_analysis(pools{m},p);
        row=make_row(names{k},methods{m},res,pools{m},p,sqrt(gain_squared(m)/gain_count(m)),applied(m));
        rows{end+1}=row; %#ok<AGROW>
        group.(methods{m})=struct('tracks',pools{m},'res',res,'record',row);
    end
    group.clean_tracks_identical=clean_same;
    d0=group.clean_median.record.D_um2_s; d1=group.clean_spectral_common.record.D_um2_s;
    group.clean_D_relative_change=d1/d0-1;
    group.single_frequency_diagnostic_available=fft_done;
    all_results.(names{k})=group;
    fprintf('\nGROUP %s clean tracks identical=%d clean D relative change=%g\n',names{k},clean_same,group.clean_D_relative_change);
end
records=[rows{:}]; details=[detail_rows{:}];
results=struct('scope','Development ablation, 3 x 1500 frames; no precision acceptance; reference truth only after inference', ...
    'seed',1,'flicker_amplitude',0.08,'flicker_frequency_Hz',423.7,'flicker_phase',0.37, ...
    'summary',records,'per_video',details,'counterexamples',counterexamples);
results.clean_tracks_identical=cellfun(@(n)all_results.(n).clean_tracks_identical,names);
results.clean_D_relative_change=cellfun(@(n)all_results.(n).clean_D_relative_change,names);
fid=fopen(fullfile(outdir,'ablation.json'),'w','n','UTF-8');
fprintf(fid,'%s',jsonencode(results,PrettyPrint=true)); fclose(fid);
writetable(struct2table(records),fullfile(outdir,'ablation_summary.csv'),'Encoding','UTF-8');
save(fullfile(outdir,'ablation_results.mat'),'results','all_results','-v7.3');
fprintf('\nCOMPLETE ablation. Clean identity flags: %d %d %d. Not an accuracy acceptance.\n',results.clean_tracks_identical);
end

function r=make_row(scenario,method,res,tr,p,gain_rms,applied)
g=res.group; use=res.sel & isfinite([tr.loc_rms_nm]); loc_rms=nan;
if any(use)
    weights=[tr(use).n_loc]; loc_rms=sqrt(sum(weights.*[tr(use).loc_rms_nm].^2)/sum(weights));
end
r=struct('scenario',scenario,'method',method,'nominal_nm',p.sample.nominal_nm, ...
    'D_um2_s',g.D_avgmsd,'diameter_nm',g.d_avgmsd,'bias_percent',100*(g.d_avgmsd/p.sample.nominal_nm-1), ...
    'n_tracks',g.n_tracks,'n_loc',g.n_loc,'active_videos',g.n_videos_with_tracks, ...
    'localization_rms_nm',loc_rms,'gain_residual_rms',gain_rms,'applied_videos',applied);
end

function err=relative_gain_error(estimated,truth)
% 只评价时间相对照明，静态尺度应由平场承担；真gain不传入估计器。
estimated=estimated(:)/mean(estimated); truth=truth(:)/mean(truth);
residual=truth./estimated; err=residual/mean(residual)-1;
end

function [tracks,info]=single_frequency_counterexample(video,p,target)
% 原test.m单频核心的诊断版本；额外保留静态平场以给它更有利条件。
% 不冒充原脚本逐行复刻：省略BMP读写/RGB/显示归一化，使用已有追踪器。
[N,~,F]=size(video); V=single(video); bg=median(V,3);
assert(all(bg(:)>0)); C=(V-bg)./bg; C=C-mean(C,3); clear V
spectrum=fft(C,[],3); freq=(0:floor(F/2))*p.sim.fps/F;
[~,index]=min(abs(freq-target)); actual=freq(index);
isDC=index==1; isNyquist=mod(F,2)==0 && index==F/2+1;
a=spectrum(:,:,index); clear spectrum
factor=2; if isDC || isNyquist, factor=1; end
phase=exp(1i*2*pi*(index-1)*(0:F-1)/F);
filtered=single(factor/F*real(a.*reshape(phase,1,1,[])));
info=struct('requested_frequency_Hz',target,'actual_frequency_Hz',actual, ...
    'frequency_resolution_Hz',p.sim.fps/F,'energy_fraction',double(sum(filtered(:).^2))/double(sum(C(:).^2)), ...
    'scope','Single stress video; static flatfield plus single-bin projection; NOT full user script or accepted estimator');
clear C
kern=rvt_kernels(N,p.ana.rvt_radii); detections=cell(F,1);
for first=1:100:F
    ix=first:min(first+99,F); R=rvt_apply(filtered(:,:,ix),kern);
    for q=1:numel(ix), detections{ix(q)}=detect_frame(filtered(:,:,ix(q)),R(:,:,q),p); end
end
tracks=link_trajectories(detections,p);
end
