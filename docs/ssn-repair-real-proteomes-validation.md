# 真实蛋白组的 SSN 修复验收

输入由用户指定：
`/home/mselab/licj/projects/running/ogprofiler-runs/s0_v1_v2_of_regression/proteomes`。
已检查为 13 个 FAA，约 57 MB。原始目录只读，不复用或覆盖旧 S0 结果。

## 执行与来源

通过项目 wrapper 提交 `slurm/verify_ssn_repair.sh PROTEOMES OF_SOURCE_ROOT`。
资源沿用 S0：56 CPU、128 GB、24 小时，无 array。
生产算法为提交 `3d75da0`；新增加的验收工具随 dirty source snapshot 保存，
其 diff、source SHA256、Slurm 参数和输入来源均进入 provenance。

1. 在 Slurm 下复制蛋白组到独立 RUN_ID 目录，比较原文件与复制品校验和。
2. 使用新源码、显式 mean 配置重新执行 prepare → search → edges。
3. 同一份新 hits 和蛋白表进入 V2 与实际 OF3.1.5 原始函数。
4. 逐层比较 max bitscore、拟合样本与参数、B、BH/RBH、cutoff、connect、方向 W。
5. 直接读取 CLI 生成的 retained_edges.parquet，验证其 score_uv/score_vu
   等于 OF 完整方向 W，weight 等于 `(Wuv+Wvu)/2`，不止验证内存重建结果。
6. 比较后复核 hits、蛋白表、生产边文件和原/复制 FASTA 的校验和。

## 验收条件

`--require-equivalence` 在任一必需层支持不同或浮点差异超容差时返回失败。
拟合样本多重集必须一致，参数比较使用既定数值容差。
缺失必需阶段同样失败；非 canonical、重复、非有限/非正生产权重直接报错。
容差沿用审计：`1e-8 + 1e-6 * max(abs(V2),abs(OF))`。

`OF_directional_W_vs_V2_symmetric_lift` 与 max 投影仍仅作诊断，
不将方向矩阵与无向图的预期差异误判为失败。

成功要求：Slurm/进程正常退出，report completed/input_integrity_verified/validation.passed
均为 true，最终 verification.json 的 passed 为 true，并且原蛋白组保持未变。
组件统计本身接近、作业正常退出，均不足以单独证明修复通过。

## 范围

同一份 hits、共享蛋白 ID 顺序下验证 OF 默认完整运行 SSN，以及 Leiden mean 投影。
不把原始 FASTA 排序和独立搜索的差异归为 SSN 实现差异。
不运行 MCL/Leiden/OG 提取，因此不据此判定最终 OG precision 或聚类等价性。

结果位于独立 RUN_ID 下的 `verification.json`、`ssn-audit/report.json`、
`ssn-audit/stage_totals.tsv` 和 `v2-run/edges/edge-manifest.json`。
异常时先读 Slurm、v2.stderr 和逐层失败列表，不盲目重复提交。

## 本次提交

- JOB_ID：`1410705`。
- RUN_ID：`20260930T160827Z_3d75da046ac4_82fa8f69_64915`。
- Source SHA256：`82fa8f698fec09e8ddef0828609c1dc2e3e1ab6297eeccfd063e63df7f215ab8`。
- 远端运行目录：`/home/mselab/licj/projects/running/ogprofiler-runs/20260930T160827Z_3d75da046ac4_82fa8f69_64915`。
- 提交前本地/远端小型测试各 24 passed，包含刻意损坏的生产权重失败检查。
- 首次查询观察到 `RUNNING`，日志已进入 fresh prepare/search/edges。

## 已完成结果（2026-10-01）

**这份真实蛋白组的 SSN 修复验收通过。**
Slurm `COMPLETED`，exit `0:0`，耗时 `00:16:01`；
report 的 completed/input_integrity_verified/validation.passed 均为 true；
verification.json 的 passed、proteomes_snapshot_verified、original_proteomes_unchanged 均为 true。

输入规模：13 物种、128,483 蛋白，重新搜索得到 10,941,213 原始 hits；
有效非自身正分 hits 为 10,812,752，重复有效 query-target 为 0。

| 验收层 | 结果 |
| --- | --- |
| 最大 bitscore 去重 | 10,812,752 项完全一致 |
| NBS 样本/参数 | 169/169 方向物种对一致；两侧样本均 540,243，参数完全相同，无拟合 warning |
| 完整方向 B | 支持完全一致，超容差项 0，最大绝对误差 8.88e-16 |
| BH | 1,381,274 个方向，集合完全一致 |
| 跨物种 RBH | 1,077,246 个方向，集合完全一致 |
| cutoff | 128,483 蛋白，超容差项 0，最大绝对误差 6.66e-16 |
| 固定 B 的 cutoff/control | 完全一致，包括近似并列 RBH 边界 |
| connect | 1,280,455 个方向，集合完全一致 |
| 方向 W | 1,394,818 非零方向，支持一致，超容差项 0，最大绝对误差 1.78e-15 |
| CLI 落盘 score_uv/score_vu | 与 OF W 支持一致，超容差项 0，最大绝对误差 1.78e-15 |
| CLI 落盘 mean 权重 | 与 OF `(W+W.T)/2` 支持一致，超容差项 0，最大绝对误差 1.78e-15 |

浮点结果并非所有项 bitwise 一致：B 有 896,578 项 exact difference，cutoff 有 9,575，
方向 W 有 103,241；误差仅浮点舍入量级，均远小于验收容差，不造成 BH/RBH/connect 差异。

最终生产无向图有 **701,947** 条边、12,616 个组件、5,005 个单点，最大组件 34,982。
这些统计也与 OF W 的无向支持一致。实际 OF 方向 W 有 1,394,818 非零项，
mean 图的诊断对称展开有 1,403,894 项；二者直接相减及 max 投影仍出现差异是预期行为，
不代表验收失败。所有必需验收层失败列表为空。

结论限定于此次数据、共享 hits/ID 顺序、OF 默认完整运行方向 SSN 及完整精度 mean 投影。
它验证了审计指出的 NBS、cutoff、W 组装及 CLI 落盘路径修复，
不是独立搜索结果、MCL 与 Leiden 聚类结果或最终 OG 准确度的等价证明。

紧凑原始结果已取回本地：
`.provenance/slurm/20260930T160827Z_3d75da046ac4_82fa8f69_64915/results/`。
其中 report.json 保存 hits/蛋白表/生产边文件校验和、参考源码 hash、依赖版本和完整阶段统计；
verification.json 保存实际来源、源码快照身份和总验收结论。大型输出保留远端。
