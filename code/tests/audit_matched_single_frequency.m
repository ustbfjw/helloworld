function result=audit_matched_single_frequency(outdir)
%AUDIT_MATCHED_SINGLE_FREQUENCY 单频诊断公平性补充，不作正式精度验收。
% 对同一1500帧压力视频比较：原median、同平场不投影、匹配照明最近bin、10Hz。
% 匹配423.7Hz仍需映射到有限帧DFT格点；本脚本不是用户test.m完整复刻。
% 它额外提供静态平场，省略BMP/RGB/显示归一化，并使用现有RVT和反演。
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root); setup_paths(root);
if nargin<1, outdir=fullfile(root,'04_results','matched_single_frequency_diagnostic'); end
if ~exist(outdir,'dir'), mkdir(outdir); end
diary(fullfile(outdir,'matched_run.log')); cleanup=onCleanup(@() diary('off')); %#ok<NASGU>
names={'mono30','mono100','mono200'}; seeds=[1 101 101]; video_ids=[1 2 2];
base=sim_params('quick'); tab=ipsf_lookup_table(base); rows={}; all_tracks=struct();
fprintf('DEVELOPMENT DIAGNOSTIC ONLY: single selected video per group; no accuracy acceptance.\n');
fprintf('All branches share the same raw stress video; +/-8%% illumination at 423.7Hz inside exposure.\n');
for k=1:3
    p=apply_scenario(base,names{k}); p.sim.seed=seeds(k);
    p.cam.flicker_amplitude=0.08; p.cam.flicker_frequency=423.7; p.cam.flicker_phase=0.37;
    p.ana.background_mode='median'; p.out.make_plots=false;
    traj=simulate_trajectories(p,seeds(k)); [video,meta]=render_iscat_video(traj,p,tab);
    assert(meta.n_sat==0,'This diagnostic requires an unsaturated video.');
    fprintf('\n%s seed=%d, frames=%d\n',names{k},seeds(k),size(video,3));
    [tr,~]=inta_track_video(video,p);
    tr=analyze_tracks(tr,traj,p,video_ids(k)); res=inta_population_analysis(tr,p);
    row=make_row(names{k},'production_median',nan,nan,1,tr,res,p);
    rows{end+1}=row; all_tracks.(names{k}).production_median=tr; %#ok<AGROW>
    print_row(row);
    F=size(video,3); V=single(video); bg=median(V,3); assert(all(bg(:)>0));
    C=(V-bg)./bg; C=C-mean(C,3); clear V video
    original_energy=double(sum(C(:).^2));
    tr=direct_tracking(C,p); tr=analyze_tracks(tr,traj,p,video_ids(k)); res=inta_population_analysis(tr,p);
    row=make_row(names{k},'flatfield_demean_no_projection',nan,nan,1,tr,res,p);
    rows{end+1}=row; all_tracks.(names{k}).flatfield_demean_no_projection=tr; %#ok<AGROW>
    print_row(row);
    spectrum=fft(C,[],3); clear C
    target_freqs=[423.7 10]; labels={'matched_illumination_bin','magnetic_10Hz_bin'};
    for j=1:2
        fa=(0:floor(F/2))*p.sim.fps/F;
        [~,index]=min(abs(fa-target_freqs(j))); actual=fa(index);
        isDC=index==1; isNyquist=mod(F,2)==0 && index==F/2+1;
        factor=2; if isDC || isNyquist, factor=1; end
        coeff=spectrum(:,:,index); phase=exp(1i*2*pi*(index-1)*(0:F-1)/F);
        filtered=single(factor/F*real(coeff.*reshape(phase,1,1,[])));
        energy=double(sum(filtered(:).^2))/original_energy;
        tr=direct_tracking(filtered,p); tr=analyze_tracks(tr,traj,p,video_ids(k)); res=inta_population_analysis(tr,p);
        row=make_row(names{k},labels{j},target_freqs(j),actual,energy,tr,res,p);
        rows{end+1}=row; all_tracks.(names{k}).(labels{j})=tr; %#ok<AGROW>
        print_row(row); clear filtered
    end
    clear spectrum traj
end
records=[rows{:}];
result=struct('scope','Development-only single-video diagnostic; not full user test.m and not precision acceptance', ...
    'method_note','Production median baseline plus the same static flatfield/demean before no projection or one conjugate DFT-bin pair; no display normalization', ...
    'flicker_amplitude',0.08,'flicker_frequency_Hz',423.7,'flicker_phase',0.37, ...
    'frames_per_video',1500,'seeds',seeds,'records',records);
fid=fopen(fullfile(outdir,'matched_results.json'),'w','n','UTF-8');
fprintf(fid,'%s',jsonencode(result,PrettyPrint=true)); fclose(fid);
writetable(struct2table(records),fullfile(outdir,'matched_summary.csv'),'Encoding','UTF-8');
save(fullfile(outdir,'matched_results.mat'),'result','all_tracks','-v7.3');
fprintf('\nCOMPLETE matched-bin diagnostic; no particle-size accuracy assertion.\n');
end

function tracks=direct_tracking(C,p)
N=size(C,1); F=size(C,3); kern=rvt_kernels(N,p.ana.rvt_radii); det=cell(F,1);
for first=1:100:F
    ix=first:min(first+99,F); R=rvt_apply(C(:,:,ix),kern);
    for q=1:numel(ix), det{ix(q)}=detect_frame(C(:,:,ix(q)),R(:,:,q),p); end
end
tracks=link_trajectories(det,p);
end

function row=make_row(name,method,target,actual,energy,tr,res,p)
g=res.group;
row=struct('scenario',name,'seed',p.sim.seed,'method',method, ...
    'nominal_nm',p.sample.nominal_nm,'target_frequency_Hz',target,'actual_frequency_Hz',actual, ...
    'D_um2_s',g.D_avgmsd,'diameter_nm',g.d_avgmsd,'n_tracks',g.n_tracks,'n_loc',g.n_loc, ...
    'unmatched_selected_tracks',sum(res.sel & [tr.truth_id]==0),'retained_temporal_energy_fraction',energy);
end

function print_row(r)
fprintf('%s: target=%g actual=%gHz D=%g d=%g tracks=%d loc=%d unmatched=%d energy=%g\n', ...
    r.method,r.target_frequency_Hz,r.actual_frequency_Hz,r.D_um2_s,r.diameter_nm,r.n_tracks,r.n_loc, ...
    r.unmatched_selected_tracks,r.retained_temporal_energy_fraction);
end
