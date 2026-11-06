function protocol=candidate_experiment_protocol()
%CANDIDATE_EXPERIMENT_PROTOCOL 空间留出候选的独立保留验证，运行前固定。
% 原六条件消融继续独立完成，本验证不得回改其结果或复用其种子。
protocol=targeted_experiment_protocol(); protocol.schema_version=2;
protocol.max_num_threads=2;
names={'static','global_offbin','same_frequency','gaussian_patch','tilted_smooth','oblique_step'};
modes={'global','global','global','gaussian_patch','tilted_smooth','oblique_step'};
freq=[423.7,423.7,10,423.7,423.7,423.7]; amps=[0,.08,.08,.08,.08,.08];
for c=1:6
    protocol.conditions(c)=struct('name',names{c},'amplitude',amps(c),'frequency',freq(c), ...
        'phase',.37,'spatial_mode',modes{c});
end
protocol.methods={'median','pn','full','spatial_validated'};
seeds=[911001,922001,933001]; % 与早期短帧接口fixture及原六条件对照分离。
for g=1:3, protocol.groups(g).seed_base=seeds(g); end
protocol.bootstrap_seed_base=51001;
protocol=rmfield(protocol,'ablation_definition');
protocol.primary_comparison='gaussian_patch and tilted_smooth: spatial_validated versus full; three diameter groups separately';
protocol.candidate_definition='4x4 spatial regions, disjoint separated train/validation pixel blocks, training-only spectral proposal/refinement, protected magnetic bands, quadratic spatial cos/sin coefficient fields, validation-only accept/reject, exact original full fallback';
protocol.candidate_thresholds=struct('grid',[4 4],'spatial_polynomial_degree',2, ...
    'minimum_global_validation_coefficient_rms',.002, ...
    'minimum_validation_energy_reduction',.5,'maximum_regional_amplitude_worsening',.002, ...
    'minimum_active_regions',2,'gain_lower_exclusive',.5,'gain_upper_inclusive',1.5);
protocol.validation_metric='sqrt(mean_regions(sum_frequencies(cos_error^2+sin_error^2))), coefficient amplitude RMS, not temporal RMS; validation pixels never adjust frequencies or fitted coefficient fields';
protocol.local_generator='Gaussian patch, sinusoidally varying smooth mask and unknown oblique discontinuity; none generated from the fitted quadratic polynomial. Exact masks are only in renderer and evaluation metadata.';
protocol.inclusion_rule='Require actual camera-video evidence of lower illumination error and preserved sizing; do not treat synthetic quadratic-field tests as end-to-end performance. Failed/boundary cases retained. No automatic patent novelty or inventiveness conclusion.';
protocol.comparison_rule='paired B=2000 video bootstrap, common indices per group including empty videos; candidate versus median/pn/full; exploratory pointwise intervals, no multiplicity adjustment; predeclared +1 pp absolute-bias loss margin';
end
