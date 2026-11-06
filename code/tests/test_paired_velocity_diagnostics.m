function metrics=test_paired_velocity_diagnostics(outdir,schema1_source,schema2_source)
%TEST_PAIRED_VELOCITY_DIAGNOSTICS 只复制已标注的统计fixture，不改正式实验结果。
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root); setup_paths(root);
if nargin<1, outdir=fullfile(fileparts(root),'paired_velocity_unit_tests'); end
if nargin<2, schema1_source=fullfile(fileparts(root),'targeted_statistics_fixture','summary.mat'); end
if nargin<3, schema2_source=fullfile(fileparts(root),'candidate_statistics_fixture','summary.mat'); end
sources={schema1_source,schema2_source}; expected=[90 54]; state=rng;
for schema=1:2
    source=sources{schema}; assert(exist(source,'file')==2); source_hash=targeted_sha256(source);
    folder=fullfile(outdir,sprintf('schema%d',schema)); if ~exist(folder,'dir'), mkdir(folder); end
    destination=fullfile(folder,'summary.mat');
    assert(~strcmpi(source,destination),'测试只能写入独立副本。');
    copyfile(source,destination,'f'); saved_hash=targeted_sha256(destination);
    rows=write_paired_velocity_diagnostics(folder);
    assert(numel(rows)==expected(schema));
    assert(strcmp(saved_hash,targeted_sha256(destination)) && strcmp(source_hash,targeted_sha256(source)));
    primary='full'; if schema==2, primary='spatial_validated'; end
    assert(all(strcmp({rows.method},primary)));
    static=rows(strcmp({rows.condition},'static'));
    assert(all([static.delta_flow_um_s]==0) && all([static.delta_magnetic_velocity_um_s]==0));
    assert(all([static.delta_flow_ci95_low]==0) && all([static.delta_flow_ci95_high]==0));
    assert(all([static.delta_magnetic_velocity_ci95_low]==0) && all([static.delta_magnetic_velocity_ci95_high]==0));
    tab=readtable(fullfile(folder,'paired_velocity.csv')); assert(height(tab)==expected(schema));
    if schema==2
        unavailable=rows(strcmp({rows.condition},'tilted_smooth'));
        assert(all(isnan([unavailable.delta_magnetic_velocity_um_s])) && ...
            all(isnan([unavailable.delta_magnetic_velocity_ci95_low])));
        sparse=rows(strcmp({rows.condition},'gaussian_patch'));
        assert(all(isnan([sparse.delta_magnetic_velocity_ci95_low])) && ...
            all(strcmp({sparse.magnetic_velocity_status},'descriptive_only_fewer_than_4_active_videos')));
    end
end

% 94.5%有效重抽样不能给区间；真值字段变动不得影响配对差或有效数。
folder=fullfile(outdir,'incomplete_bootstrap'); if ~exist(folder,'dir'), mkdir(folder); end
s=load(schema2_source,'report','bootstrap'); report=s.report; bootstrap=s.bootstrap;
methods=cellstr(string(report.protocol.methods)); im=find(strcmp(methods,'spatial_validated'));
B=report.protocol.bootstrap_replicates; bootstrap{1}.draws.v(1:ceil(.055*B),1,im)=NaN;
report.summary(1).truth_velocity_um_s=1e9;
save(fullfile(folder,'summary.mat'),'report','bootstrap');
rows=write_paired_velocity_diagnostics(folder);
r=rows(strcmp({rows.group},report.protocol.groups(1).name) & strcmp({rows.condition},'static'));
assert(all(isnan([r.delta_magnetic_velocity_ci95_low])) && ...
    all(strcmp({r.magnetic_velocity_status},'insufficient_valid_paired_bootstrap')));
assert(all([r.magnetic_velocity_paired_bootstrap_valid]==B-ceil(.055*B)));
assert(all([r.delta_flow_ci95_low]==0) && all([r.delta_magnetic_velocity_um_s]==0));
assert(isequal(state,rng),'事后配对诊断不得改变全局随机流。');
metrics=struct('passed',true,'scope','copied synthetic statistics fixtures only', ...
    'schema1_rows',90,'schema2_rows',54,'zero_difference_ci',true, ...
    'unidentifiable_and_low_count_ci_nan',true,'insufficient_bootstrap_ci_nan',true, ...
    'truth_not_used_for_selection',true,'source_summaries_unchanged',true,'rng_unchanged',true);
fid=fopen(fullfile(outdir,'paired_velocity_tests.json'),'w','n','UTF-8'); assert(fid>=0);
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(metrics,PrettyPrint=true));
fprintf('PASS paired velocity diagnostics: schema1 90 rows, schema2 54 rows, NaN rules, immutable source.\n');
end
