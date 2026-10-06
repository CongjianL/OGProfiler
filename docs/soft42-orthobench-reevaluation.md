# Open Orthobench：已选 soft42 的独立复评

本轮不调参、不更新生产默认。使用 `presets/embleya-soft42.yaml`，仅显式将
runtime.workers 从 1 改为 56；subtree_workers=1、seed=42 及全部科学参数保持不变。
单作业沿用此前完整 Orthobench 的 56 CPU / 250 GB / 72 h，无 job array。

数据来源是已完成 job1410751 的 repaired-mean：12 物种、251378 蛋白。
复用冻结搜索和 mean-SSN，复制 input/edges/components 到独立运行目录，逐文件
校验哈希，全部层级组件从头计算。此次是冻结上游的层级与 OG 复评，不是重新
运行 DIAMOND 的端到端实验。prepare 检查上游科学参数一致，搜索线程数不参与匹配。

执行入口 `slurm/soft42_orthobench.sh`，参数为原 P6 运行根目录。由 dev/slurm-submit
创建源码快照并记录 Git/源码哈希/资源/参数。保留原 preset 和实际 run.yaml，记录
线程覆盖、数据来源、环境与树 manifest。历史树和结果均不覆盖。

沿用已验证的 hierarchy_regression score：
- soft42 的 v1_compatible 为主要结果；terminal 仅为切分诊断对照。
- 与同一冻结输入的历史 repaired-mean V2、独立完整 V1 结果比较。
- 使用冻结官方评分器；输出 precision/recall/F1、全输入覆盖率、RefOG 覆盖率、
  逐 RefOG 分裂/误合并诊断。完整导出参与评分，不先筛参考基因。
- 先通过完整树 gate、OG parity、ID唯一性/官方解析器 roundtrip，再报告结果。
- 结果保存在 h5-metrics/report.json；通用报告标签 new_bounded_hierarchy 指本轮 soft42。

Slurm COMPLETED 与科学结论分开：完成后读取实际分数及诊断再评估，不预先声称提升。
Open Orthobench 是本轮跨数据集复评，不反向选择 soft42 参数。

## job1411546 结果（2026-10-06）

- RUN_ID：20261005T232747Z_8bcccc9d3302_512d821e_28231。
- 源码 8bcccc9d3302c86ba899471fe17c5bb12b592b3f，dirty=false。
- sacct：COMPLETED / 0:0，耗时 07:08:44。
- evaluation_completed、strategy_parity_passed、fixed_ssn_unchanged、benchmark_unchanged
  均为true；完整树 resolved，UNRESOLVED=0，候选预算检查通过。
- parity检查30165组件，其中10712组件有层级树，其余单例单独持久化。
- 所有评分输出覆盖251378蛋白与全部raw RefOG成员（两项覆盖率均100%），
  original ID映射与官方reader roundtrip均通过。

| 方法 | 官方Precision | 官方Recall | 官方F1 | OG数 |
|---|---:|---:|---:|---:|
| soft42 / v1_compatible（主结果） | 92.6961% | 33.8220% | 49.5608% | 67916 |
| 历史repaired-mean V2 / v1_compatible | 0.0588% | 80.8705% | 0.1176% | 33283 |
| 历史完整V1 | 78.1215% | 45.5011% | 57.5075% | 58218 |
| soft42 / terminal（诊断对照） | 62.5849% | 1.0901% | 2.1429% | 216820 |

soft42相比历史完整V1：Precision +14.5745个百分点、Recall -11.6791个百分点、
F1 -7.9467个百分点。V1是独立完整流水线对照，不是仅层级参数变化的因果对照。
历史V2极低精确率由本轮相同官方评分器重新得到，不是百分比/比例转换错误。

70个RefOG诊断：soft42为37 SPLIT、22 MERGE_AND_SPLIT、8 EXACT、3 OVERMERGE；
V1为29 SPLIT、25 MERGE_AND_SPLIT、12 EXACT、4 OVERMERGE。
soft42的RefOG碎片关联总数488，中位4.5；V1为409、中位3。
macro best-group F1：soft42 71.5764%、V1 76.5947%（辅助诊断，非官方pairwise F1）。
soft42相对V1的best-group F1：10个RefOG提升、27持平、33下降。
下降较大者包括RefOG019（12 vs 2碎片，best F1 .2667 vs .9600）、
RefOG007（16 vs 2，.2727 vs .9244）、RefOG060（17 vs 2，.4545 vs .9926）。
这支持当前不足主要表现为过度切分/召回不足，但尚未把原因归结到特定评分或树参数。

结论：保留已选soft42作为当前Embleya方案；Orthobench证明其高精确率优势，
但F1尚低于完整V1，不据此宣称通用最佳。不同参考体系的分数不直接横比。
本轮不扩大参数搜索、不新增过滤门槛、不修改生产默认、不自动提交后续实验。
如未来重启研发，应先明确召回提升目标，并定位损失来自树表示还是OG选择层。

紧凑原始指标与诊断汇总：`docs/soft42-orthobench-1411546-results.json`。
完整小报告缓存于本地 `.benchmark_cache/soft42-orthobench-1411546/`；
大结果仍留在上述远端RUN_ID的独立目录，历史实验均保留。
