function acceptance = run_acceptance(output_dir)
%RUN_ACCEPTANCE 融合版完整验收：单元/数值、静态回归、同视频周期照明消融。
% 串行通常需要一至数小时，取决于计算机；不需要额外MATLAB工具箱。
% 失败会assert并保留日志；此脚本不会为了接近真值自动改变模型或参数。
root=fileparts(mfilename('fullpath')); setup_paths(root);
if nargin<1
    stamp=char(datetime('now','Format','yyyyMMdd_HHmmss'));
    output_dir=fullfile(root,'04_results',['acceptance_run_' stamp]);
end
if exist(output_dir,'dir')
    entries=dir(output_dir); entries=entries(~ismember({entries.name},{'.','..'}));
    assert(isempty(entries),'验收输出目录必须为空，请指定新目录，避免旧通过记录被误认为本次结果。');
end
if ~exist(output_dir,'dir'), mkdir(output_dir); end
diary(fullfile(output_dir,'acceptance_run.log'));
cleanup=onCleanup(@() diary('off')); %#ok<NASGU>
test_analysis;
test_physics_imaging;
illumination_tests=test_common_illumination(output_dir);
numerical_tests=test_numerical_convergence(output_dir);
diary('off');
names={'mono30','mono100','mono200'};
% 沿用上一版固定基准；明确属于回归/压力复测，而非全新独立保留样本。
seeds=[10001 20001 20001];
for k=1:3
    run_validation_case(names{k},fullfile(output_dir,'clean'),seeds(k));
    run_background_validation_case(names{k},fullfile(output_dir,'periodic'),seeds(k));
end
diary(fullfile(output_dir,'acceptance_run.log'));
clean=summarize_acceptance(fullfile(output_dir,'clean'));
periodic=summarize_background_validation(fullfile(output_dir,'periodic'));
acceptance=struct('passed',clean.passed && periodic.passed && ...
    illumination_tests.pass && numerical_tests.pass, ...
    'clean',clean,'periodic',periodic,'illumination_tests',illumination_tests, ...
    'numerical_tests',numerical_tests);
fid=fopen(fullfile(output_dir,'fusion_acceptance.json'),'w','n','UTF-8');
fprintf(fid,'%s',jsonencode(acceptance,PrettyPrint=true)); fclose(fid);
assert(acceptance.passed,'融合版完整验收未通过。');
fprintf('PASS fusion acceptance: static regression and paired periodic illumination stress.\n');
end
