# Embleya 13 蛋白组：以冻结 OrthoFinder3 OG 为参考

日期：2026-10-04。按用户要求，本轮不使用 OpenBench/Orthobench 成绩评价修改。
OrthoFinder3 分区是方法一致性参考，不作为独立生物学真值。

## 输入与版本边界

- 参考：远端 S0 `real_embleya/orthofinder/Results_Sep20/Orthogroups/`，
  `Orthogroups.tsv` + `Orthogroups_UnassignedGenes.tsv`。本轮不混用 N0 HOG。
- 13 物种、128,483 蛋白；按 `(species_name, original_id)` 对齐，校验 species_id/
  protein_id 映射与唯一性。OF assigned 和 unassigned 的不交并集与 V2 输入完全相等。
- OF 已分组 122,954 蛋白、13,124 OG；另有 5,529 未分组蛋白。
- 已有完整 V2 是 job1410721/P5，基础源码 f69c14cdfe9b 的冻结快照，**不是当前 HEAD**。
  SSN 来源为 job1410705 的已验收 mean 投影。P5 最终16,178 OG，无未分配蛋白。
- 当前 HEAD 为405ce3e。当前默认 kway/depth20 与实验 soft/depth42 分别评估；
  另一数据集上近期完成的 H5 结果不代替本数据集结果。

## 计分口径

主结果限定在 OF assigned 蛋白全集，预测组投影到此全集；因此不使用 OF 未分组
蛋白构造主指标的负例。另报每个 OF unassigned 蛋白视为独立 singleton 的敏感性分析，
这仅是一项显式假设。缺失预测独立记为 singleton，并单独计数；重复/未知身份报错。
计分采用 contingency table 组合计数，不枚举全部蛋白对。

pair precision/recall 是“同 OG 蛋白对”的一致性，不是物种对正交关系准确率。
同时报告 B-cubed、macro best-group F1、exact/split/merged 数量，避免大组支配结论。
指标工具：`benchmarks/og_extraction/embleya_reference.py`。

## P5 历史锚点：实测结果

| 指标 | P5 final OG | P5 terminal families |
|---|---:|---:|
| Pair precision | 5.9327% | 5.8953% |
| Pair recall | 92.1115% | 91.4903% |
| Pair F1 | 11.1475% | 11.0769% |
| B-cubed precision | 73.0181% | 73.0594% |
| B-cubed recall | 94.5483% | 93.7207% |
| Macro best-group F1 | 75.0227% | 74.5290% |
| 完整匹配参考 OG | 8,229 | 7,992 |
| 被拆分参考 OG | 1,102 | 1,358 |

主口径：773,634参考组内对，712,606保留，61,028丢失；
11,298,812预测同组对落在不同参考OG。825个预测组混合多个参考OG。
这不是“所有家族都很差”：8,229/13,124参考OG完整匹配，pair指标主要受大组支配。
跨参考组的预测同组对：前5组占62.4302%，前10占79.4091%，前30占94.3930%。
OF未分组作为singleton时，pair P/R/F1 为5.8518%/92.1115%/11.0046%。

## 第一偏离定位：五个大 OG 是搜索失败终止叶

读取 P5 component0 的真实 nodes/members 表，并按 protein_id 与最终OG交叉核对：

| 最终 OG | 全部成员数 | 组件:叶节点 | depth | 原停止原因 |
|---|---:|---|---:|---|
| OG000000001 | 2,481 | 0:1 | 1 | UNSTABLE |
| OG000000035 | 1,568 | 0:12 | 1 | UNSTABLE |
| OG000000078 | 1,411 | 0:17 | 1 | GAMMA_LIMIT |
| OG000000048 | 1,378 | 0:13 | 1 | GAMMA_LIMIT |
| OG000000026 | 1,366 | 0:9 | 1 | UNSTABLE |

每个最终OG恰好来自一个完整terminal leaf，非最终OG跨叶合并所造。
最大组在主计分全集中剩2,470蛋白，混合227个OF OG，贡献3,031,868个不一致合并对。
旧层级将搜索失败作为普通TERMINAL，后续发布为大OG，是这五个案例的直接机制。
此证据不证明所有其余差异都来自相同机制。

参考组内对只有9对跨V2连通组件，组件边界允许的乐观召回上界99.99884%。
这排除了“组件断开”是主瓶颈，但不排除组件内部边权/mean投影影响划分。
P5 terminal已有11,298,324个不一致合并对；final仅增加488个，
因此首先修层级失败处理/分裂选择，而不是先重写通过parity的OG提取。

## 当前版本评估：提交记录（结果见下节）

- JOB_ID：1411104。
- RUN_ID：20261004T081659Z_405ce3ea83d6_0366c702_12341。
- HEAD：405ce3ea83d617d57c0a496197f6755c529f9a81，dirty=1；新增评估代码、测试及脚本
  均保存在提交快照，源码SHA256为
  0366c702c2154ba1d25dabe5257dbdc8ad578d80ec9415c351fb1a61cc9073a0。
- 脚本：`slurm/embleya_of_reference.sh`；使用项目 `dev/slurm-submit`。
- 56 CPU、128 GB、24小时上限，与P5资源一致；无array，两个配置顺序运行。
- 配置A：当前 DEFAULT_CONFIG（kway/depth20）；配置B：同配置显式soft/depth42。
  B为组合实验，不将差异单独归因于binary或depth；不迁移生产默认。
- 只读复用P5 input/edges/components，不重搜、不重建SSN、不复用旧hierarchy缓存。
  先验证SSN manifest、身份全集，复制OF参考；分别构建新hierarchy、注释、提取、导出、计分。
  搜索失败遵守当前UNRESOLVED门禁，失败配置不产生正式OG分数，也不以缺失当好成绩。
- 记录参考/输入/边的hash、完整配置、环境、分阶段time/stdout/stderr和结束后输入完整性。
- 本地5个计数fixture测试通过、Ruff及shell语法检查通过；作业开头另跑小型测试。
  本文在提交后撰写，不属于该次计算源码快照。

## 修正优先级与验收

1. **先验证已实现的失败语义修正**：原UNSTABLE/GAMMA_LIMIT大叶应真正解析，
   或保留UNRESOLVED阻断发布。阻断不是科学准确度改善，要单独报告完成率。
2. **对剩余大合并做定点层级诊断**：跟踪上述蛋白集合在新树的去向，检查候选覆盖、
   稳定性、gamma边界、子群比例、停止原因和事件资格；保留所有参数与失败证据。
   不按OF标签直接强拆/合并，不仅因组大而删除，不先降低ARI或盲目提高gamma。
3. **再处理错误拆分**：按损失对数列出参考OG，区分跨SSN组件、首次hierarchy分离、
   事件判定与最终OG选择；例如P5 OG0000005 的61蛋白被分到5组、丢失1,386组内对。
4. **条件性隔离图/聚类差异**：若同组件内大合并仍持续，进一步固定hits比较
   OF方向图/MCL与V2 mean图/Leiden，而不是直接重跑搜索或宣布SSN完全无关。
5. **避免单指标/单数据过拟合**：同时看pair与macro/B-cubed、exact/split/merge、
   大组贡献、唯一覆盖、unresolved、资源消耗。按组件划分开发/保留验证集合后再调参；
   本轮不以OpenBench分数作修改的选择标准。

本地紧凑证据：`.provenance/embleya-of-reference/`，包括P5两个comparison JSON、
top-merge-terminals JSON、输入映射及OG TSV。大型搜索和图保留远端。

## 2026-10-05 验收：job1411104 完成

远端 accounting 实测：COMPLETED、ExitCode 0:0、耗时06:33:17，
调度时间2026-10-04 16:15:33–22:48:50（Asia/Shanghai）。
两配置的 hierarchy-all/annotate-network/orthogroups/export 均正常退出；
各自层级日志7611组件completed、0 failed、5005单点，OG阶段12616组件完成。
当前OG发布门禁通过；本轮未另外重读所有节点表统计终止原因。
远端评估器fixture为5 passed。input-integrity.json的8项输入/边/参考检查全部true。

取回两个members.tsv和run.yaml，逐项匹配计分报告记录的SHA256；
以本地既有身份表/参考表重新计算，两个计分口径的全部标量指标均复现。
两个实际配置的差异仅topology_policy及max_depth；其余配置相同。
本地接收检查记录在该RUN_ID的results/local-acceptance.json。

### 主口径结果（OF assigned 122,954蛋白）

| 指标 | P5历史 | 当前default | 当前soft42 |
|---|---:|---:|---:|
| Pair precision | 5.93% | 90.23% | **95.71%** |
| Pair recall | 92.11% | 24.09% | **69.08%** |
| Pair F1 | 11.15% | 38.03% | **80.24%** |
| B-cubed precision | 73.02% | 97.50% | 97.01% |
| B-cubed recall | 94.55% | 37.25% | 81.47% |
| Macro best-group F1 | 75.02% | 59.68% | **91.14%** |
| 完整匹配参考OG | 8,229 | 4,329 | **8,600** |
| 被拆分参考OG | 1,102 | 8,442 | 3,989 |
| 不一致合并对 | 11,298,812 | 20,175 | 23,934 |
| 丢失参考组内对 | 61,028 | 587,265 | 239,231 |

所有配置最终完整分配128,483蛋白；组数等完整输出统计如下（区别于主口径投影组数）：

| 指标 | P5历史 | 当前default | 当前soft42 |
|---|---:|---:|---:|
| 最终OG总数 | 16,178 | 82,819 | 25,126 |
| 当前实测singleton OG | — | 70,983 | 9,464 |
| 最大OG | 2,481 | 37 | 40 |

OF未分组视作独立singleton的敏感性口径，default P/R/F1为
89.3291%/24.0901%/37.9467%，soft42为95.1311%/69.0770%/80.0371%；结论方向一致。

### 差距已从超大合并转为剩余过度拆分

1. 两个当前版本都消除了P5千蛋白大组；并非仅靠阻止发布，实际完成了全量输出。
2. default过度拆分严重：70,983单蛋白OG、召回24.09%；其macro best-group F1
   还低于P5，因此不应仅看pair F1就宣称所有方面改善。
3. soft42相对default净增加348,034保留参考对：新恢复350,147对，同时失去2,113对。
   不一致合并对净增3,759：新增6,732、消除2,973；pair precision提高，
   但B-cubed precision略降，残余合并仍需追踪。
4. soft42仍丢失239,231参考对，其中仅9对跨SSN组件，剩余239,222对位于组件内部。
   当前证据定位到“组件内部层级/事件/OG选择”，尚未逐对区分三者。
   不能沿用另一数据集的RefOG事件归因作为本数据集结论。
5. soft42是拓扑政策与深度预算的组合实验，不将收益单独归因于binary或depth42。
   同时P5与当前版本跨越多项算法变更，历史比较也不是单变量因果实验。

### 下一轮修正应以首次分离诊断为先

优先检查最大剩余拆分案例，而非继续整体提高gamma/降低ARI：

| OF参考OG | 蛋白数 | default碎片数 | soft42碎片数 | soft42丢失组内对 |
|---|---:|---:|---:|---:|
| OG0000000 | 153 | 22 | 18 | 10,001 |
| OG0000001 | 117 | 46 | 32 | 6,414 |
| OG0000002 | 112 | 42 | 28 | 5,906 |

对这些成员定位两棵树的首次分离节点、已选candidate、binary/kway/fallback、
事件标签及OG选择trace，区分：
- 过早层级分离：诊断有限候选覆盖与选择目标，保持准入/预算可审计。
- 有共同适格祖先却未合并：检查OG选择/消耗顺序，并用冻结V1策略区分实现偏差与规则局限。
- 缺少适格祖先：检查层级与事件规则的衔接，以及OF OG和V1式输出家族的定义差异。
  大型多拷贝家族尤其应先验证语义，避免把方法定义差异直接解释为生物学错误。
同时保留soft42新增不一致合并对作为反例约束，不直接按OF标签拼接预测组。

资源观察：default层级1:07:11，soft42层级1:32:25；两者OG阶段约77/80分钟，
导出约11分钟。最大阶段RSS分别931,152/947,364 KiB；仅为本次阶段计量，非独立性能基准。

结论：本数据上soft42明显优于当前default的一致性表现，仍有显著召回缺口。
保留soft42作为后续诊断基线，不自动迁移生产默认。此次仅取回、核验、分析和更新文档，
未修改科学算法、未重提作业、未使用OpenBench结果。
