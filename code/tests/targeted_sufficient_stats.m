function stats=targeted_sufficient_stats(tracks,res,p)
%TARGETED_SUFFICIENT_STATS 提取单视频可相加统计量，不按真值筛选推断。
% res.sel和每滞后allfinite(M)口径与inta_population_analysis完全一致。
% 最后四个匹配/定位字段仅作事后诊断，不参与D、d或漂移参数估计。
nl=p.ana.n_lag;
stats=struct('counts',zeros(nl,1),'moments',zeros(nl,6), ...
    'driftA',zeros(2),'driftB',zeros(2,1),'n_tracks',0,'n_loc',0, ...
    'unmatched_selected',0,'matched_localizations',0, ...
    'localization_squared_sum_nm2',0);
assert(numel(res.sel)==numel(tracks),'轨迹与选择标记长度不一致。');
for k=find(res.sel)
    t=tracks(k);
    M=[t.msd,t.msd_perp,t.wave_q2,t.wave_qdt,t.wave_cross,t.flow_cross];
    ok=all(isfinite(M),2) & t.npairs>0;
    stats.counts(ok)=stats.counts(ok)+t.npairs(ok);
    stats.moments(ok,:)=stats.moments(ok,:)+M(ok,:).*t.npairs(ok);
    stats.driftA=stats.driftA+t.drift_A;
    stats.driftB=stats.driftB+t.drift_b;
    stats.n_tracks=stats.n_tracks+1; stats.n_loc=stats.n_loc+t.n_loc;
    stats.unmatched_selected=stats.unmatched_selected+(t.truth_id==0);
    if t.truth_id>0 && isfinite(t.loc_rms_nm)
        stats.matched_localizations=stats.matched_localizations+t.n_loc;
        stats.localization_squared_sum_nm2=stats.localization_squared_sum_nm2+t.n_loc*t.loc_rms_nm^2;
    end
end
end
