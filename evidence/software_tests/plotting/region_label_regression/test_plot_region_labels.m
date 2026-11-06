function record=test_plot_region_labels(root)
%TEST_PLOT_REGION_LABELS 实际导出图的刻度回归，不用合成图覆盖正式结果。
% 在图窗关闭前读取真实坐标轴；明确断言标签数、顺序和左右含义。
maxNumCompThreads(2);
root=char(root); addpath(fullfile(root,'work','postprocess'));
testdir=fullfile(root,'work','postprocess_label_regression');
if ~isfolder(testdir), mkdir(testdir); end
inputs={fullfile(root,'work','targeted_statistics_fixture'), ...
    fullfile(root,'work','candidate_statistics_fixture'), ...
    fullfile(root,'work','targeted_ablation_results_v2')};
outputs={fullfile(testdir,'fixture_schema1'),fullfile(testdir,'fixture_schema2'), ...
    fullfile(root,'work','research_release_assets','figures','targeted')};
run_names={'schema1_fixture','schema2_fixture','schema1_formal'};
schemas=[1 2 1];
region_expected={{'左半图423.7 Hz（左半区）','左半图423.7 Hz（右半区）'}, ...
    {'高斯局部斑（左半区）','高斯局部斑（右半区）', ...
     '平滑非均匀场（左半区）','平滑非均匀场（右半区）', ...
     '斜向阶跃（左半区）','斜向阶跃（右半区）'}};
global_expected={{'静态','全局423.7 Hz','保护带内11.4 Hz', ...
    '保护带外13.4 Hz','同频10 Hz','左半图423.7 Hz'}, ...
    {'静态','全局423.7 Hz','同频10 Hz','高斯局部斑','平滑非均匀场','斜向阶跃'}};
observed=struct('run',{},'kind',{},'tick_count',{},'labels',{});
active_run='';
old_delete=get(groot,'DefaultFigureDeleteFcn');
cleanup=onCleanup(@()set(groot,'DefaultFigureDeleteFcn',old_delete)); %#ok<NASGU>
set(groot,'DefaultFigureDeleteFcn',@capture_labels);
records=cell(1,3); input_checks=cell(1,3);
prior_path=fullfile(outputs{3},'schema1_plot_provenance.json');
prior_provenance=jsondecode(fileread(prior_path));
for k=1:3
    active_run=run_names{k}; before=input_hashes(inputs{k});
    records{k}=plot_study_summary(inputs{k},outputs{k});
    after=input_hashes(inputs{k});
    assert(isequal(before,after),'labeltest:InputChanged','绘图改变了输入CSV或协议。');
    input_checks{k}=struct('run',active_run,'unchanged',true,'sha256',after);
    assert(records{k}.fixture==(k~=3),'labeltest:Fixture','正式/fixture标记错误。');
    for kind={'region','global'}
        rows=observed(strcmp({observed.run},active_run) & strcmp({observed.kind},kind{1}));
        assert(numel(rows)==3,'labeltest:AxisCount','每组图必须捕获三个粒径分面。');
        if strcmp(kind{1},'region'), expected=region_expected{schemas(k)};
        else, expected=global_expected{schemas(k)}; end
        for j=1:numel(rows)
            assert(rows(j).tick_count==numel(expected),'labeltest:TickCount','刻度数错误。');
            assert(isequal(rows(j).labels,expected),'labeltest:Labels','实际刻度标签数、顺序或含义错误。');
            assert(~any(contains(string(rows(j).labels),newline)), ...
                'labeltest:Newline','刻度标签仍包含换行。');
        end
    end
    fprintf('LABEL_AXIS_PASS %s: region=%d ticks x 3 axes; global=6 ticks x 3 axes\n', ...
        active_run,numel(region_expected{schemas(k)}));
end
% 显式释放持有嵌套回调的根属性，防止清理闭包延至MATLAB退出才执行。
set(groot,'DefaultFigureDeleteFcn',old_delete);
record=struct('pass',true,'matlab_version',version, ...
    'generated_at',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss Z')), ...
    'test_scope','Actual exported axes labels before deletion; fixture and formal inputs read-only', ...
    'helper_sha256',file_hash(fullfile(root,'work','postprocess','plot_study_summary.m')), ...
    'test_sha256',file_hash([mfilename('fullpath') '.m']), ...
    'region_axes_checked',9,'global_axes_checked',9,'observed_axes',observed, ...
    'input_checks',{input_checks},'plot_records',{records}, ...
    'prior_formal_provenance',prior_provenance);
fid=fopen(fullfile(testdir,'region_label_regression.json'),'w','n','UTF-8');
assert(fid>=0,'labeltest:Write','无法写测试记录。');
fprintf(fid,'%s',jsonencode(record,'PrettyPrint',true)); fclose(fid);
fprintf('REGION_LABEL_REGRESSION_PASS\n');

    function capture_labels(fig,~)
        axes_list=findall(fig,'Type','axes');
        for q=1:numel(axes_list)
            ax=axes_list(q); ylabel_text=string(ax.YLabel.String);
            if isequal(ylabel_text,"区域照明增益残余 RMS（%）"), kind='region';
            elseif isequal(ylabel_text,"全图照明增益残余 RMS（%）"), kind='global';
            else, continue; end
            labels=reshape(cellstr(string(ax.XTickLabel)),1,[]);
            observed(end+1)=struct('run',active_run,'kind',kind, ...
                'tick_count',numel(ax.XTick),'labels',{labels}); %#ok<AGROW>
        end
    end
end

function hashes=input_hashes(folder)
names={'summary.csv','pairwise.csv','protocol.json','protocol.sha256'};
hashes=cell(1,numel(names));
for k=1:numel(names), hashes{k}=file_hash(fullfile(folder,names{k})); end
end

function hash=file_hash(path)
fid=fopen(path,'rb'); assert(fid>=0,'labeltest:Read','文件无法读取。');
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid), md.update(fread(fid,1024*1024,'*uint8')); end
hash=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
