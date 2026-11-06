# 交变梯度磁场下磁纳米颗粒 iNTA：MATLAB 工程与验证资料

> 2026-11-06科研交付：两轮完整计算与独立数值核对已完成。新增空间模式有适用边界；不宣称全条件精度通过或专利创造性已成立。

本工程分别模拟纯水中的 **30 nm、100 nm、200 nm 三组独立单分散样品**，每组只含一种标称粒径。处理顺序为三维布朗运动与磁泳、曝光内 iSCAT 光强积分、带相机噪声的二维图像、照明与背景处理、盲定位和链接、磁漂移修正、扩散系数及水动力直径估计。三组不混合，不用真实粒径标签给检测轨迹分群。

研究基础为 Kashkanova 等的 *Precision size and refractive index analysis of weakly scattering nanoparticles in polydispersions*，Nature Methods 19, 586–593 (2022)，DOI：[10.1038/s41592-022-01460-z](https://www.nature.com/articles/s41592-022-01460-z)。本工程是针对磁驱动及三组单分散样品的适配仿真，不是原论文全部装置与折射率分析的逐项复现。

**建议先阅读已附报告和结果表，无需重新运行数小时的仿真。** `quick` 只用于检查流程；本轮正式证据来自两个独立冻结协议，不来自单次演示图。

## 1. 发布包目录

以下是交付目录；两个代码快照均保留冻结的相对结构及原始字节。

```text
发布包根目录/
├── README.md
├── code/                       候选工程冻结快照：87 个 .m 与 references/
├── baseline_snapshot/          原六条件消融冻结快照：70 个 .m
├── evidence/
│   ├── targeted/               原六条件完整 PN 与三组件消融证据
│   ├── candidate/              空间候选的独立保留验证证据
│   └── software_tests/         单元、管线、渲染、统计与绘图测试记录
├── postprocess/                独立审计/绘图/示例重建辅助工具及说明
├── docs/                       算法技术说明、完整验收报告、文献与创新边界
├── figures/                    由完整汇总数据生成的正式图表
└── examples/                   有来源记录的固定视频示例
```

`evidence/targeted/` 和 `evidence/candidate/` 分别保留 `protocol.json`、`protocol.sha256`、`parameters.mat`、每组逐视频 MAT/日志、`summary.csv`、`pairwise.csv`、`video_metrics.csv` 及汇总 MAT/JSON。正式证据是否完整，以最终发布清单和审计报告为准。

`evidence/software_tests/` 中带 `fixture`、`synthetic` 的内容是合成数据自测，不属于正式性能结果。不要把这些图表混入 `figures/` 的正式结果。

## 2. 环境与路径

已验证环境为 **MATLAB R2024b，基础 MATLAB**；主算法不要求额外统计或图像处理工具箱。保留默认 JVM 以进行 SHA-256 校验，不使用 `-nojvm` 启动。其他版本尚未完成同等兼容性验收。代码为 UTF-8，主要注释和图注为中文。

先把 MATLAB 当前文件夹切换到解压后的**发布包根目录**，只需在这里设置一次：

```matlab
releaseRoot = pwd;
assert(isfolder(fullfile(releaseRoot,'code')) && ...
       isfolder(fullfile(releaseRoot,'baseline_snapshot')), ...
       '当前文件夹应为发布包根目录。');
```

**不要使用 `addpath(genpath(releaseRoot))`。** 两个快照含同名函数；递归加入整个包可能调用错误版本。推荐在两个独立 MATLAB 会话运行两轮复现。若在同一会话切换，先回到 `releaseRoot`，恢复默认路径并清除函数缓存，再只加入所需快照：

```matlab
cd(releaseRoot);
restoredefaultpath;
clear functions;                 % 保留releaseRoot变量，清除旧函数缓存
codeRoot = fullfile(releaseRoot,'code');
addpath(codeRoot);
setup_paths(codeRoot);           % 显式加入该工程的模块
assert(strcmp(which('main_run_all'),fullfile(codeRoot,'main_run_all.m')));
```

上述路径设置用于当前 MATLAB 会话，不需执行 `savepath`。运行途中不要编辑冻结快照、替换函数或向快照目录加入新的 `.m` 文件；新试验输出放在发布包根目录下独立的 `reproduction/` 目录。

## 3. 先看已附结果

在已经设置 `releaseRoot` 后，可以直接打开图表/报告，或只读结果表：

```matlab
targetedTable = readtable(fullfile(releaseRoot,'evidence','targeted','summary.csv'), ...
    'TextType','string');
candidateTable = readtable(fullfile(releaseRoot,'evidence','candidate','summary.csv'), ...
    'TextType','string');
disp(candidateTable(:,{'group','condition','method','diameter_nm', ...
    'd_ci95_low','d_ci95_high','gain_rms_pct','active_videos','n_videos'}));
```

`diameter_nm` 是组内位移矩汇总后估计 D、再经 Stokes–Einstein 换算的水动力直径；不是单条轨迹直径的简单平均。`gain_rms_pct` 是指定乘性照明残差的逐像素指标，不是相机总噪声。`active_videos` 是有有效位移对的视频数，`n_videos` 包括空视频。NaN 或缺失置信区间应保留为证据不足，不能改成 0。

**正式结果摘要：**

| 项目 | 本轮完整结果 |
|---|---|
| 两轮协议、源文件/参数 hash、逐视频完整性审计 | 各140/140检查点，70/87份源码hash一致，独立数值核对无错误；见 `evidence/audits/` |
| 原六条件主比较：全局 423.7 Hz，full 对 PN | full为30.1608/100.5743/203.7996 nm；照明RMS低于PN，三组探索性非劣，无粒径优越证据 |
| 三组件消融及保护区/同频条件 | 频率细化改善照明；空间门控有取舍；近频/同频保护保留照明盲区，见原六条件附录 |
| 候选主要比较：Gaussian / smooth，空间候选对 full | 平滑条件候选29.4873/99.4387/205.1768 nm，照明改善且探索性非劣；高斯斑全部回退无新增效果 |
| 局部边界及空间接受/回退次数 | 平滑接受20/20、49/60、42/60全部视频；active中20/20、9/16、5/17。200 nm阶跃−7.20%未达目标 |
| 最终适用范围与推荐入口 | 原默认保留；可选spatial_validated只支持所测保护带外平滑周期扰动的有限改进，不能保证任意局部照明 |
| 完整验收报告 | [验证报告](docs/验证报告.md)、[算法说明](docs/算法技术说明.md)、[技术交底摘要](docs/技术交底摘要.md) |

运行完整与性能合格是两件事；结果表行数完整不代表所有条件均满足 ±5% 误差目标，也不代表某方法优于 PN。

## 4. 小样本运行入口

按第 2 节只加入 `code/` 后运行：

```matlab
maxNumCompThreads(2);

% 全局共同周期照明方法：每组3段、每段1500帧，检查整个流程。
demoOut = fullfile(releaseRoot,'reproduction','demo_global_01');
results = main_run_all('quick',demoOut,'periodic','spectral_common');

% 可选：空间留出候选；不满足局部模型条件时回退原全局方法。
spatialDemoOut = fullfile(releaseRoot,'reproduction','demo_spatial_01');
resultsSpatial = main_run_all('quick',spatialDemoOut,'periodic','spatial_validated');
```

第四个参数可取 `median`、`spectral_common`、`spatial_validated`；不传时沿用工程配置的 `spectral_common`。第三个参数只支持 `static` 或 `periodic`；这里的 `periodic` 是全局 ±8%、423.7 Hz 照明，磁驱动仍为 10 Hz。因此这个演示不包含全部局部照明条件；空间候选回退是正常的模型选择行为，不能据此解释为程序未执行。

`quick` 不作精度验收。将第一个参数改成 `full` 会使用 20/60/60 段、每段 5000 帧的规模，但**仍不是本轮两套专用的新种子对照协议**。不同运行使用不同输出目录；`main_run_all` 不会自动阻止覆盖同名演示输出。

## 5. 重新运行原六条件正式消融

若只阅读已附正式结果，跳过本节。若需独立重算，使用 **`baseline_snapshot/`**，不要从 `code/` 调用同名 `run_targeted_ablation`。

```matlab
cd(releaseRoot);
restoredefaultpath;
clear functions;
baselineRoot = fullfile(releaseRoot,'baseline_snapshot');
addpath(baselineRoot);
setup_paths(baselineRoot);
assert(strcmp(which('run_targeted_ablation'), ...
    fullfile(baselineRoot,'run_targeted_ablation.m')));

targetedRun = fullfile(releaseRoot,'reproduction','targeted_01');
run_targeted_ablation('prepare',targetedRun);   % 新目录，只执行一次
run_targeted_ablation('mono30',targetedRun);
run_targeted_ablation('mono100',targetedRun);
run_targeted_ablation('mono200',targetedRun);
run_targeted_ablation('summarize',targetedRun); % 三组完整后才汇总
```

这一轮包括：中值背景、完整全帧 PN、完整全局法 `full`、去空间一致性门控、去磁频保护、去格点内细化，共 6 种方法；条件为静态、全局 423.7 Hz、保护带内 11.4 Hz、保护带外 13.4 Hz、同频 10 Hz、左半图 423.7 Hz。三个粒径的种子基数为 111001、222001、333001，按视频递增 100。实际已输出108行汇总和90行配对比较。

runner 内部使用 4 个计算线程。默认按上述命令串行运行三组，避免额外并行占用内存。若中断，保持同一代码和输出目录，重新调用未完成的组即可；不要再次 `prepare`，也不要编辑协议以绕过校验。只有完整且兼容的逐视频记录会被跳过。

## 6. 重新运行空间候选的独立保留验证

这一轮使用 **`code/`**。建议打开新的 MATLAB 会话；若继续使用当前会话，应重新执行以下路径切换，并确认 `releaseRoot` 仍指向发布包根目录：

```matlab
cd(releaseRoot);
restoredefaultpath;
clear functions;
codeRoot = fullfile(releaseRoot,'code');
addpath(codeRoot);
setup_paths(codeRoot);
assert(strcmp(which('run_candidate_validation'), ...
    fullfile(codeRoot,'run_candidate_validation.m')));

candidateRun = fullfile(releaseRoot,'reproduction','candidate_01');
run_candidate_validation('prepare',candidateRun); % 新目录，只执行一次
run_candidate_validation('mono30',candidateRun);
run_candidate_validation('mono100',candidateRun);
run_candidate_validation('mono200',candidateRun);
run_candidate_validation('summarize',candidateRun);
```

这一轮包括中值背景、完整 PN、原 `full`、`spatial_validated`，共 4 种方法；条件为静态、全局 423.7 Hz、同频 10 Hz、高斯局部斑、平滑非均匀场、斜向阶跃。后三种局部形状不是按候选的二次空间模型生成。三个粒径的新种子基数为 911001、922001、933001。实际已输出72行汇总和54行配对比较。runner 内部使用 2 个计算线程，断点续跑规则与上一轮相同。

两轮均为每组 20/60/60 段、每段 5000 帧，采用同视频方法配对及 2000 次独立视频 bootstrap。两轮是不同的实验，不得将结果或 checkpoint 混入同一个目录。重新 `prepare` 会记录当前时间与运行环境，参数中的工程目录也可能不同，因此协议 hash 不要求与发布包的原始 hash 相同；包内原始证据及其 hash 应原样保留，不覆盖。

## 7. 汇总图与示例

正式汇总经完整性审计后，可直接用 `postprocess/` 重新绘图，不需要重新模拟相机视频：

```matlab
addpath(fullfile(releaseRoot,'postprocess'));
plot_study_summary( ...
    fullfile(releaseRoot,'evidence','targeted'), ...
    fullfile(releaseRoot,'reproduction','figures_targeted_01'));
plot_study_summary( ...
    fullfile(releaseRoot,'evidence','candidate'), ...
    fullfile(releaseRoot,'reproduction','figures_candidate_01'));
```

该工具只读完整 CSV/协议，输出中文 300 dpi PNG、矢量 PDF 和来源 hash；缺 CI、缺尺寸、对数轴真零、左右区域及空间接受次数均明确标记，不自动把统计 status 写成优越或验收结论。固定示例的生成方式、所用视频索引及一致性检查见 `postprocess/` 与 `examples/` 中实际保存的说明和来源记录；示例图不能替代全部视频汇总。

## 8. `run_acceptance` 的定位

**`run_acceptance.m` 是历史回归入口，不是本次新种子六条件消融或空间候选的独立保留验收。** 它执行早期单元/数值测试、静态回归及既有周期压力条件，沿用历史基准种子。文件名中的 acceptance 或历史报告中的 PASS 不能替代第 5、6 节的证据。

如确需运行历史回归，也应指定新目录，不覆盖包内证据；通常阅读 `evidence/software_tests/` 已经足够。两个代码快照只复制冻结的 `.m` 文件，`code/` 另附 `code/references/` 中的 5 个来源文件，不包含旧 README、旧方法说明或旧验收报告；**本发布包根目录 README、最终技术说明与本轮完整验收报告共同界定此次结论**。

## 9. 适用范围与解释边界

- 主尺寸估计使用垂直于固定磁梯度方向的位移 MSD，斜率为 2D，避免将沿梯度磁漂移直接当作布朗扩散。它依赖已知固定方向、各向同性体相扩散及当前稀悬浮模型，不能自动推广到磁成链、旋转梯度、近壁各向异性或未知横向流。
- 全局法估计空间共同的周期乘性照明；空间候选只在固定训练/留出像素检查及安全增益门槛满足时采用低阶空间场，否则精确回退原全局法。留出像素互不重叠，但不保证无颗粒污染或统计独立。
- **同频盲区仍存在。** 两种方法均保护第 1 至第 5 倍磁驱动频率（含基频）附近；空间候选不会自动解除保护，因此没有实现照明与真实磁运动的同频盲分离。
- 原 `test.m` 的逐像素单频重建未被强行用作 MSD 输入。保留的是频谱辅助估计照明的思想，校正作用于原始逐帧视频。
- 当前相机生成器是标量 Debye/Mie、有限曝光、高计数高斯散粒近似及 ADC 模型；尚无真实装置、真实相机标定或真实样品的完整验收。默认包覆层厚度为 0，实际水动力直径不必等于厂家核直径。
- 完整 PN 基线使用全视场逐帧总功率，时间中值只提供常数尺度；PiSCAT 公开实现采用时间均值。该常数在理想静态对比度中抵消，但不能称为 PiSCAT 源码逐位复刻。
- 软件测试通过、照明误差下降、尺寸精度保持和专利创造性是不同问题。仿真结果不能保证中国发明专利授权，也不能证明已优于全部现有算法；已有组件与本项目具体组合的区别及检索边界见 `docs/` 和 `code/references/`。

解释性能时请同时保留无改善、拒绝校正、边界失效和证据不足的条件，不以最接近 30/100/200 nm 的输出反向选择算法、种子或阈值。
