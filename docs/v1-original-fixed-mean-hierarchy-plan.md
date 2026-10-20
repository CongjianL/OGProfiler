# 固定mean SSN的V1原始hierarchy对照

2026-10-04；针对job1410868后剩余低召回与新增污染，用户要求优先RefOG012/015/019/036。

## 实验边界

- 从冻结job1410868 `new-hierarchy/components` 读取原字节mean SSN。
- 完整递归组件470（46蛋白）、943（31）、3976（12）、500（45）；
  另保留28085/26818的两个SSN isolated蛋白。不是从已分散子节点开始搜索。
- 用冻结 `legacy/OGProfiler_v1.py`，SHA256
  `35607fb7cdbd5f8db6c6a3bc3f22df953884b24af2df023f79307bdfcdc4ac3b`。
  AST原定义直接执行HHN.__init__/SpecificBoard/BipartiteGraphs/RunCommunityDetection，
  配套GetAttribution/MappingGammaForCC；只观测代理la.find_partition，无方法重写。
- rber、同mean边权、V1 gamma coefficient1、原始每次10轮、seed参数缺省。
  保留小节点exactly2目标、重复首点评估、原1000轮循环与失败terminal、
  2蛋白直接拆分及ONE_SPECIES停止。没有V2稳定性/比例门槛、24点预算或深度42截断。
- 外层串行BFS替代V1多进程队列调度，仅编排原方法。组件vertex顺序使用当前固定SSN顺序，
  不是历史V1重建SSN后的全局vertex顺序；记录图ID映射、版本和所有实际membership。
- 不给V1偷偷加seed。单组件1次seed-free realization；原算法随机且调度/软件环境影响结果，
  本对照不声称历史bitwise重放，也不凭一次成功推断算法普遍恢复。
- V1原始事件/OG选择在新V1树上用既有冻结oracle执行，species_overlap=0；
  当前V2树和OG只读取，不重新执行、不修改OG/scorer、停止尺寸或稳定性阈值。

## 输出及验收

1. 每个组件全部节点/成员/终止理由、原始事件、V1实际OG组与unassigned。
2. 每次Leiden调用的节点、gamma、seed缺省、10轮、count、quality、membership和耗时，
   包括V1重复评估与失败点；不只保存最终binary。
3. 按protein ID/clade匹配新V1树与V2树，定位根到叶的首次完整成员分区偏离。
4. 四个RefOG全部真值对逐对比较SSN阶段、各树LCA事件、同叶/跨叶与最终保留状态，
   及best-group P/R/F1、额外成员；明确这是局部诊断，不是全数据官方指标。
5. hierarchy严格子群、子群完整唯一覆盖及最终OG无重复校验；输入checksum前后不变。
6. 原始V1算法失败terminal单独标记，不当作V2政策完成，不移植到生产。

工具 `benchmarks/og_extraction/v1_original_hierarchy.py`；脚本
`slurm/v1_original_mean_ssn_targets.sh`。
资源为新小诊断1 CPU/8G/1小时，无array，与H5资源无关。
通过dev/slurm-submit不可变快照执行。现有dirty工作树精确捕获diff、untracked与hash，
不将run冒充HEAD的clean版本。来源H5固定且输出在新run，原产物不改动。

本地及remote fixture测试7 passed，覆盖V1二蛋白零调用、首点重复/缺省seed、
1001调用后原始失败terminal与现有召回/污染审计。后续等待job完成再分析真实结果。

## 首轮已提交

- Slurm job **1411083**。
- RUN_ID `20261004T025041Z_27a2533e1334_28f01e46_79637`。
- HEAD `27a2533e13340b9192e39251dc18db3af028fc7c`，**dirty=1**，不是clean HEAD实验。
- Snapshot SHA256 `28f01e46e17314702b15f67369649e4d75775ee965c24af8ce5b8be103cb36c3`。
- 精确源码/patch/untracked/hash与参数记录在本地`.provenance/slurm/RUN_ID`及远端run的provenance。
- 原始mean SSN来自冻结job1410868；本轮没有提交全数据V1 hierarchy或官方benchmark。

2026-10-04首轮完成：job1411083 COMPLETED/0:0，11秒。四个完整组件验收通过，
具体结果见 [同mean SSN的原始V1对照结果](v1-original-fixed-mean-hierarchy-results.md)。
019比V2多保留54真值对；012相同，015同污染，036反而少保留44对。
