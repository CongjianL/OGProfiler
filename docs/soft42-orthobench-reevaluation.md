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
