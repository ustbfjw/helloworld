# 独立后处理工具

本目录包含三个原始辅助文件及一个展示标题修正工具，不属于冻结的 MATLAB 仿真源码清单。它们用于读取既有结果、制作实例或独立核对，不改变正式结果，也不自动判断算法效果、创新性或专利授权条件。

以下命令从交付根目录执行。最终目录对应关系是：候选源码 `code/`，原全局消融源码 `baseline_snapshot/`，两套正式结果 `evidence/candidate/` 与 `evidence/targeted/`。新输出使用单独目录，不覆盖交付中的图或数据。

**汇总绘图**

```matlab
release_root = pwd;
addpath(fullfile(release_root, 'postprocess'));
plot_study_summary( ...
    fullfile(release_root, 'evidence', 'targeted'), ...
    fullfile(release_root, 'local_review', 'targeted_figures'));
plot_study_summary( ...
    fullfile(release_root, 'evidence', 'candidate'), ...
    fullfile(release_root, 'local_review', 'candidate_figures'));
```

`plot_study_summary.m` 读取正式 CSV 和协议，校验行数、组合、区间及计数口径。它不重新计算置信区间，不将空值画作零，也不把软件 fixture 当作正式证据。输入包含 fixture 标记时，输出会标为 FIXTURE。

**单视频实例重建**

```matlab
release_root = pwd;
addpath(fullfile(release_root, 'postprocess'));
manifest = create_candidate_examples( ...
    fullfile(release_root, 'code'), ...
    fullfile(release_root, 'evidence', 'candidate'), ...
    fullfile(release_root, 'local_review', 'candidate_examples'));
```

`create_candidate_examples.m` 检查完整检查点、参数和源码哈希；按预设的轨迹可用性/局部模型接纳规则选实例，不按粒径误差选图。它以相同参数和种子重建六条件，只执行候选生产路径，并逐字段与检查点统计比较。数值不一致时报告错误；没有合格轨迹的组明确跳过。输出目录必须全新。原相机摘录最多 201 帧，明确标为展示摘录，不能替代完整视频的验收统计。重建会使用最多 2 个 MATLAB 计算线程，耗时和内存明显高于只读取 CSV 的绘图。

**独立 Python 核对**

```powershell
python postprocess/audit_formal_results.py evidence/targeted --source-root baseline_snapshot --output postprocess/local_targeted_audit.json
python postprocess/audit_formal_results.py evidence/candidate --source-root code --output postprocess/local_candidate_audit.json
```

`audit_formal_results.py` 依赖 NumPy、SciPy；检查冻结来源、140 个检查点和汇总 CSV，并独立计算部分点估计。输出必须位于该脚本所在的 `postprocess/` 目录，且不得位于被审计目录内。该审计不重新计算 bootstrap 或区间，也不判断某方法更优。命令中的两个新 JSON 名称用于本地复核，不替代随交付保存的正式审计记录。

前三个工具的来源和复制哈希见 `../evidence/software_tests/source_manifest.json`。其中 `asset_relative` 记录整理时的 assets 路径；打包会把其 `software_tests/` 前缀映射为最终 `evidence/software_tests/`。

**展示标题修正**

`fix_example_title_display.m` 是本次交付新增的纯绘图工具：原冻结绘图函数会把条件名中的下划线解释为下标，此工具从保存的18个 `example.mat` 读取数据，重画相机/C/RVT三联图，并把面板标题解释器设为 `none`。它核对三幅图的CData与原数组逐值一致，且MAT文件哈希不变；不重新模拟或推断。已附图完成该修正，记录为 `../examples/display_title_record.json`。

只有自行重建新的实例目录后才需要调用，不能对已含修正记录的目录重复覆盖：

```matlab
fix_example_title_display( ...
    fullfile(release_root,'code'), ...
    fullfile(release_root,'local_review','candidate_examples'));
```

该新工具的哈希在标题修正记录及交付总清单中；原 `create_candidate_examples.m` 与原重建清单保持原字节，未将新修正冒充原重建过程。
