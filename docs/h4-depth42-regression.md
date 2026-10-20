# H4-only depth42 回归结果

2026-10-03；job1410845；源码801b68b，干净不可变快照。
对照job1410810，仅max_depth20→42；其他准入、停止尺寸、SSN/seed及OG/scorer不变。

## 验收

- Slurm COMPLETED/0:0，11:51:40→13:30:19（Asia/Shanghai），1:38:39。
- component0的69,642蛋白唯一完整覆盖，unresolved=0，实际最大深度33。
- 旧depth<20的98,429节点：成员clade、节点字段、子分区及全部候选trace保持不变。
- 新增777节点（上界854）、15,189次调用（上界30,744）；完整总调用1,797,261。
- 99,664节点，46,043分裂、53,621政策叶：49,659 SINGLETON、3,962 ONE_SPECIES。
- 34,481 BINARY、9,189 KWAY、2,373 FALLBACK_KWAY；74次refinement截断。
- 完整双进程回放nodes/members/candidates相同；manifest、成员覆盖、叶尺寸、每节点候选
  预算、robust3调用计数、选中候选原准入及fallback证据全部通过。
- 串行命令墙时27:54.96、最大RSS994,788KB；双进程17:18.56、最大RSS1,261,412KB。

依据是冻结作业验收报告与日志，已取回小型JSON核查；本轮未在本地重新加载大型SSN。
报告保存于`.provenance/slurm/20261003T035305Z_801b68b38892_b93b0699_26447/results/`：
`h4-acceptance.json`、`depth-prefix-acceptance.json`、`h4-h5-completion.json`及time/log。

## 结论边界与下一步

本轮消除了component0的全部未解析，并验证深度20之前拓扑/候选前缀不变。
H5状态为NOT_REQUESTED_H4_ONLY，scoring_started=false；没有新的P/R/F1或准确度结论。
生产默认max_depth仍20，深度42只在显式实验配置使用。

下一步应显式决定将soft+depth42用于H5全组件实验（不是默认迁移）。复用已验证component0，
固定其manifest与输入，处理其余组件后验收完整hierarchy；有任何未解析先审计，不静默加深。
完整hierarchy resolved后再运行实际V1 OG parity与官方Orthobench，比较P/R/F1、coverage及大组污染。
