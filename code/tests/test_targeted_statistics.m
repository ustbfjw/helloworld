function record=test_targeted_statistics()
%TEST_TARGETED_STATISTICS 充分统计与原群体估计的数值等价、空视频回归。
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root); setup_paths(root);
p=sim_params('quick'); p.ana.fit_model='perp'; p.ana.bootstrap_n=0;
p.ana.min_len=25; p.ana.min_pairs=5; p.sim.n_videos=8;
p.cam.exposure=0; p.cam.Npx=1024; p.cam.px=1e-6;
p.part.n_per_size=[0 1 0]; p.mag.dir=[1 1 0]/sqrt(2);
p.mag.wave='sine'; p.mag.f=13; p.mag.phase0=0.37;
rs=RandStream('mt19937ar','Seed',9317); alltracks=empty_tracks();
for iv=[1 2 4 5 6 7]
    F=200+iv; times=(1:F)'/p.sim.fps;
    pos=cumsum(sqrt(4/p.sim.fps)*randn(rs,F,2),1)+ ...
        (9*field_integral(0,times,p)+2*times)*p.mag.dir(1:2);
    tr=analyze_tracks(make_raw(pos,p),[],p,iv);
    tr.truth_id=iv; tr.loc_rms_nm=2+iv/10;
    if iv==2
        % 原实现要求六个位移矩全部有限，这一滞后不得误入汇总。
        tr.wave_q2(1)=nan;
    end
    alltracks=[alltracks,tr]; %#ok<AGROW>
end
short=analyze_tracks(make_raw(zeros(12,2),p),[],p,1);
alltracks=[alltracks,short];
state=rng; result=inta_population_analysis(alltracks,p);
emptyres=inta_population_analysis(empty_tracks(),p);
zero=targeted_sufficient_stats(empty_tracks(),emptyres,p); total=zero;
active=0;
for iv=1:p.sim.n_videos
    subset=alltracks([alltracks.video]==iv);
    res=inta_population_analysis(subset,p);
    s=targeted_sufficient_stats(subset,res,p);
    assert(s.n_tracks==sum(res.sel));
    if ismember(iv,[3 8])
        assert(s.n_tracks==0 && s.n_loc==0 && ~any(s.counts));
    end
    active=active+any(s.counts>0);
    names=fieldnames(total);
    for j=1:numel(names), total.(names{j})=total.(names{j})+s.(names{j}); end
end
estimate=targeted_estimate_stats(total,p); g=result.group;
actual=[estimate.D_um2_s,estimate.diameter_nm,estimate.intercept_um2, ...
    estimate.flow_um_s,estimate.magnetic_velocity_um_s];
expected=[g.D_avgmsd,g.d_avgmsd,g.b_avgmsd,g.flow_mean,g.v_mean];
assert(all(abs(actual-expected)<=1e-11*max(1,abs(expected))));
assert(isequal(total.counts,g.n_pairs) && total.n_tracks==g.n_tracks && total.n_loc==g.n_loc);
assert(active==6 && active==g.n_videos_with_tracks);
assert(total.unmatched_selected==0 && total.matched_localizations==total.n_loc);
expected_squared=sum([alltracks(result.sel).n_loc].*[alltracks(result.sel).loc_rms_nm].^2);
assert(abs(total.localization_squared_sum_nm2-expected_squared)<1e-8);
emptyestimate=targeted_estimate_stats(zero,p);
assert(isnan(emptyestimate.D_um2_s) && isnan(emptyestimate.diameter_nm));
negative=total; tau=(1:p.ana.n_lag)'/p.sim.fps;
negative.moments(:,2)=(1-2*tau).*negative.counts;
negativeestimate=targeted_estimate_stats(negative,p);
assert(negativeestimate.D_um2_s<0 && isnan(negativeestimate.diameter_nm), ...
    '非正D必须保留，不能裁剪成正值或生成有限直径。');
assert(isequal(state,rng),'统计helper不得改变全局随机数状态。');
record=struct('test','targeted_statistics','passed',true,'n_videos',8, ...
    'active_videos',active,'n_tracks',total.n_tracks,'n_loc',total.n_loc, ...
    'D_um2_s',estimate.D_um2_s,'diameter_nm',estimate.diameter_nm, ...
    'max_absolute_equivalence_error',max(abs(actual-expected)), ...
    'empty_videos_retained',true,'negative_D_not_clipped',true);
fprintf('test_targeted_statistics PASS: active=%d/8 max_error=%.3g D=%.9g d=%.9g\n', ...
    active,record.max_absolute_equivalence_error,estimate.D_um2_s,estimate.diameter_nm);
end

function r=make_raw(pos,p)
r.frames=(1:size(pos,1))';
r.x=pos(:,1)/(p.cam.px*1e6)+(p.cam.Npx+1)/2;
r.y=pos(:,2)/(p.cam.px*1e6)+(p.cam.Npx+1)/2;
r.C=0.1*ones(size(r.x)); r.n_loc=size(pos,1);
end
