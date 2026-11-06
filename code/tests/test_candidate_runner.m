function metrics=test_candidate_runner(outdir)
%TEST_CANDIDATE_RUNNER 小尺寸端到端记录/缓存/保存/续跑检查，不作粒径证据。
% 仅修改本测试输出目录里的协议和参数；正式协议仍为5000帧20/60/60段。
root=fileparts(fileparts(mfilename('fullpath'))); setup_paths(root);
run_candidate_validation('prepare',outdir);
protocol=jsondecode(fileread(fullfile(outdir,'protocol.json')));
s=load(fullfile(outdir,'parameters.mat'),'parameters'); parameters=s.parameters;
for g=1:3
    parameters{g}.sim.n_frames=100; parameters{g}.sim.n_videos=1;
    parameters{g}.sim.L_xy=.5e-6; parameters{g}.sim.z_init=[.9,1.1]*1e-6;
    protocol.groups(g).n_videos=1;
    % 软件fixture使用专用种子，不消耗正式保留验证的种子。
    parameters{g}.sim.seed=880001+1000*(g-1);
    protocol.groups(g).seed_base=parameters{g}.sim.seed;
end
protocol.frames=100; protocol.fixture_only=true;
save(fullfile(outdir,'parameters.mat'),'parameters');
protocol.parameters_sha256=targeted_sha256(fullfile(outdir,'parameters.mat'));
write_text(fullfile(outdir,'protocol.json'),jsonencode(protocol,PrettyPrint=true));
write_text(fullfile(outdir,'protocol.sha256'),targeted_sha256(fullfile(outdir,'protocol.json')));
for name={'mono30','mono100','mono200'}
    run_candidate_validation(name{1},outdir);
    path=fullfile(outdir,name{1},'video_001.mat'); before=targeted_sha256(path);
    run_candidate_validation(name{1},outdir);
    assert(strcmp(before,targeted_sha256(path)),'续跑改写了已完成视频。');
    a=load(path); assert(isequal(size(a.items),[6 4]));
    assert(all(arrayfun(@(x)isfield(x.stats,'counts') && isfield(x.metric,'gain_rms'),a.items),'all'));
    assert(any(arrayfun(@(x)~isempty(x.reused_from),a.items),'all'),'应覆盖精确缓存。');
    sample=load(fullfile(outdir,name{1},'sample.mat'),'samples');
    assert(numel(sample.samples)==6 && all(cellfun(@isstruct,sample.samples)));
end
metrics=struct('passed',true,'groups',3,'conditions',6,'methods',4,'frames',100, ...
    'checkpoint_and_resume',true,'scope','synthetic runner fixture, not formal performance data');
write_text(fullfile(outdir,'candidate_runner_test_record.json'),jsonencode(metrics,PrettyPrint=true));
disp(metrics);
end

function write_text(path,value)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',value);
end
