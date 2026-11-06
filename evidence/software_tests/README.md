# 软件与机制自测证据

本目录归档已经产生的测试记录，全部按原始字节复制。**它们证明的范围由各记录明确限定：接口、数值一致性、错误处理或合成机制行为；不能证明真实实验有效，也不能代替正式视频统计、科学性能评价或专利判断。** `passed`、`pass` 及其他结论字段保留原值，没有根据本次整理重新赋值。

本目录在准备区名为 `software_tests/`，最终打包位置为 `evidence/software_tests/`。正式的两套结果位于同级 `../targeted/` 和 `../candidate/`，不在此处混排。

| 子目录 | 现存记录覆盖的范围 | 解释边界 |
|---|---|---|
| `targeted_renderer/` | 六条件渲染与原生成器的逐值一致性、局部左右区域、噪声/RNG、输入检查；照明误差指标参考值核对 | 渲染或指标正确，不表示照明估计器、粒径反演已获性能验收 |
| `targeted_ablation/` | 原全局模块默认一致性、三个开关和合成条件、代码检查；小型 MAT 单元结果 | 开关行为与合成数据测试；不能替代新种子正式消融 |
| `targeted_pipeline_tests/` | 外置归一化与生产管线、完整 PN 定义的一致性 | 居中初始粒子的接口 fixture，仅作数值集成核对 |
| `candidate_renderer_tests/` | 三种保留形状生成、曝光/加噪顺序、元数据及 RNG | 生成器参考核对，不是候选估计器的保留条件效果 |
| `candidate_application_tests/` | 空间/全局增益应用、分块与 single 算序、增益误差指标、多频接口 | 只测试应用及推断后指标 |
| `spatial_candidate/` | 静态/全局精确回退、已知平滑场、多频、留出区恶化拒绝、磁频保护、输入与 RNG | 已知合成模型的机制测试；其照明残差数值不是正式相机视频的粒径效果 |
| `candidate_pipeline_tests/` | 内外空间增益路径与回退路径的一致性 | 使用后 ADC 构造增益的软件 fixture，不当作真实曝光物理验证 |
| `targeted_runner/`、`candidate_runner/` | 三组各 100 帧的小型端到端保存、缓存及续跑；附各自 fixture 协议和小型参数 MAT | 明确为 fixture；不附逐视频检查点或 sample MAT，也不是正式 5000 帧统计 |
| `statistics/targeted/`、`statistics/candidate/` | 汇总的合成统计 fixture，含空视频、负 D、NaN、有效样本与配对口径等测试 | 只保留 `test_record.json`；未复制合成 summary、bootstrap 输出或 140 个小检查点 |
| `plotting/` | 候选空间/回退诊断绘图日志；两套汇总绘图自测、缺值/区间/log 与零值处理、输入拒绝和来源记录 | 图形软件测试；本归档不重复附带 fixture 的图片/PDF，正式图表另见交付的 `figures/` |
| `paired_velocity/` | 两 schema 的配对差、零差区间、不可辨识/样本不足时 NaN、来源不改写 | 仅现存测试 JSON；没有预设速度非劣界限，也不作运动保真结论 |
| `independent_python/` | Python 审计在合成统计数据上的自测，含缺失、协议/参数/源码哈希错误和缓存核对 | 记录明确 `formal_results_audited=false`；正式审计 JSON 另置于 `../audits/` |
| `example_reproduction/` | 100 帧实例重建清单，以及缺检查点/统计篡改拒绝、全组跳过和 RNG 恢复记录 | 清单标为 fixture；两有轨迹组采用全局回退，无轨迹组跳过。未复制这些测试实例的 MAT、相机摘录或图片 |
| `failed_initial_attempt/` | 首次定向运行的报错协议、日志、参数和 runner 文本 | 失败历史，0 个正式检查点；不能计入任一成功结果 |

**记录完整性与版本边界**

`source_manifest.json` 记录每份复制文件的工作区相对来源、准备区相对目标、字节数和 SHA-256；不把它当成测试运行清单或性能判定。README 是本次整理说明，不冒充原始运行输出。后续发布总清单会另行覆盖整个交付的文件哈希。

候选绘图测试的现存结果是控制台日志，没有单独的通过状态 JSON。配对速度和 Python 自测有结果 JSON，但没有在本归档中找到对应的专门控制台日志；没有补写日志。实例重建的第一轮清单形成后，仅移除了辅助脚本中多余的代码分析抑制注释，随后静态检查为零问题；此事在 `example_reproduction/guard_record.json` 的 `comment_only_suppression_cleanup` 中记录，不能把早期清单里的 helper 哈希冒充当前文件哈希。

如来源记录中提及未归档的 fixture 输出文件名，这是当时测试过程的来源描述；本目录有意仅保留必要的记录和少量参数/单元 MAT。重新执行测试需要使用相应冻结源码并指定新的独立输出目录，不能写入正式证据目录。
