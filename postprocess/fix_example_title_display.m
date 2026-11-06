function record=fix_example_title_display(source_root,example_root)
%FIX_EXAMPLE_TITLE_DISPLAY 仅修正已保存实例面板标题中的下划线显示。
% 不重算轨迹、相机、照明或粒径；从example.mat恢复三幅原像素数组。
% 对应完整的18个非空正式实例。新实例重建后可再调用本函数。
oldpath=path; oldvisible=get(groot,'DefaultFigureVisible');
cleanup=onCleanup(@() restore_environment(oldpath,oldvisible));
addpath(source_root,'-begin'); setup_paths(source_root);
assert(strcmpi(which('plot_example_frames'),fullfile(source_root,'02_iscat_imaging','plot_example_frames.m')));
set(groot,'DefaultFigureVisible','off');
manifest=jsondecode(fileread(fullfile(example_root,'example_manifest.json')));
assert(strcmp(manifest.status,'complete') && ~manifest.fixture_only);
recordfile=fullfile(example_root,'display_title_record.json');
assert(~isfile(recordfile),'请勿覆盖已有标题修正记录。');
record=struct('scope','Display titles only; inference and MAT files unchanged', ...
    'protocol_sha256',manifest.protocol_sha256,'original_replay_helper_sha256',manifest.helper_sha256, ...
    'display_helper_sha256',targeted_sha256([mfilename('fullpath') '.m']), ...
    'all_image_cdata_exact',true,'all_mat_hashes_unchanged',true,'items',[]);
items=cell(1,18); n=0;
for g=1:numel(manifest.groups)
    group=manifest.groups(g);
    for c=1:numel(group.conditions)
        condition=group.conditions(c);
        folder=fullfile(example_root,group.name,condition.name);
        matfile=fullfile(folder,'example.mat'); before=targeted_sha256(matfile);
        sample=load(matfile,'ex','pp','stats','fixture_only');
        assert(~sample.fixture_only && ~isempty(sample.ex.frame));
        tag=[group.name '_' condition.name]; sample.pp.out.dir=folder;
        pngfile=fullfile(folder,['fig2_example_frames_' tag '.png']);
        previous_png=targeted_sha256(pngfile);
        plot_example_frames(sample.ex,sample.pp,tag); fig=gcf;
        axes_list=findall(fig,'Type','axes'); seen=false(1,3);
        expected={double(sample.ex.raw),sample.ex.C,sample.ex.R};
        for k=1:numel(axes_list)
            ax=axes_list(k); caption=string(ax.Title.String);
            if contains(caption,'原始相机图像'), idx=1;
            elseif contains(caption,'背景扣除'), idx=2;
            elseif contains(caption,'径向方差'), idx=3;
            else, continue; end
            im=findall(ax,'Type','image');
            assert(isscalar(im) && isequal(im.CData,expected{idx}));
            ax.Title.Interpreter='none'; seen(idx)=true;
        end
        assert(all(seen));
        sgtitle(sprintf('%s：单视频实例，合格轨迹%d；检测圈为本帧检测位置', ...
            tag,sample.stats.n_tracks),'Interpreter','none');
        print(fig,pngfile,'-dpng','-r150'); close(fig);
        assert(strcmp(before,targeted_sha256(matfile)));
        n=n+1; items{n}=struct('group',group.name,'condition',condition.name, ...
            'mat_sha256',before,'previous_png_sha256',previous_png, ...
            'png_sha256',targeted_sha256(pngfile),'image_cdata_exact',true);
        fprintf('标题修正 %d/18：%s，三幅CData一致，MAT哈希不变\n',n,tag);
    end
end
assert(n==18); record.items=[items{:}]; record.status='complete';
fid=fopen(recordfile,'w','n','UTF-8'); assert(fid>=0);
closer=onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(record,PrettyPrint=true));
end

function restore_environment(oldpath,oldvisible)
path(oldpath); set(groot,'DefaultFigureVisible',oldvisible);
end
