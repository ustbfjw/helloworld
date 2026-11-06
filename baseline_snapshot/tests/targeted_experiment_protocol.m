function protocol=targeted_experiment_protocol()
%TARGETED_EXPERIMENT_PROTOCOL 正式运行前固定的定向对照；不得按结果换种子或调参。
protocol.schema_version=1;
protocol.max_num_threads=4;
names={'static','global_offbin','near_inside','near_outside','same_frequency','local_left_half'};
freq=[423.7,423.7,11.4,13.4,10,423.7]; amp=[0,.08,.08,.08,.08,.08];
for c=1:numel(names)
    mode='global'; if c==6, mode='left_half'; end
    protocol.conditions(c)=struct('name',names{c},'amplitude',amp(c), ...
        'frequency',freq(c),'phase',.37,'spatial_mode',mode);
end
protocol.methods={'median','pn','full','no_spatial_gate','no_magnetic_guard','no_refinement'};
gn={'mono30','mono100','mono200'}; ds=[30,100,200]; nv=[20,60,60]; seeds=[111001,222001,333001];
for g=1:3
    protocol.groups(g)=struct('name',gn{g},'diameter_nm',ds(g),'n_videos',nv(g),'seed_base',seeds(g));
end
protocol.seed_stride=100; protocol.frames=5000;
protocol.bootstrap_replicates=2000; protocol.bootstrap_seed_base=41001; protocol.ni_margin_pp=1;
protocol.primary_comparison='global_offbin: full versus pn; three diameter groups reported separately';
protocol.inference='perpendicular MSD; displacement-pair-weighted pooled D; Stokes-Einstein; no truth-based selection';
protocol.pn_definition='sum(single(raw),1,2) per frame, convert to double, divide by temporal median; raw/gain then unchanged median pipeline';
protocol.ablation_definition=struct( ...
    'no_spatial_gate','disable only coherence and amplitude-CV rejection; retain minimum amplitude', ...
    'no_magnetic_guard','disable coarse and refined magnetic-frequency exclusion; retain frequency deduplication', ...
    'no_refinement','use DFT peak frequency without sub-bin optimization; retain all other gates');
protocol.illumination_metric='per-pixel temporal population std of (true exposure-averaged gain / estimated gain), after its own temporal-mean normalization; spatial RMS then equal-video RMS; left/right separately';
protocol.pairing='same trajectory/background and same shot/read standard-normal draws across conditions; identical raw video across methods; shot-noise scale depends on each condition';
protocol.cache_rule='reuse downstream result only when effective single(gain) vectors are exactly equal within the same raw condition';
protocol.evidence_rule='active videos >=10 and selected localizations >=5000; CI requires both methods >=4 active videos and >=95% valid replicates; retain empty videos';
protocol.comparison_rule='paired video bootstrap B=2000, same indices for every condition/method in a group; pointwise exploratory intervals, no multiplicity adjustment; NI loss=absolute full relative bias minus absolute comparator relative bias in percentage points, margin +1 pp';
protocol.acceptance='computational completeness and integrity are separate from sizing performance; report all cases, including rejection, no improvement, >5% bias and insufficient evidence';
protocol.truth_scope='generator truth is used only after inference for illumination, localization, velocity and nominal-size evaluation';
protocol.scope='synthetic scalar Debye/Mie camera model; no claim of experimental validation or patent inventiveness';
end
