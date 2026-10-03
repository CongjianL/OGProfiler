# ADR0005：soft hierarchy 深度资源预算的显式对照

- Status：**Accepted for H4-only experiment（2026-10-03用户“进入下一轮”）；全局默认仍20。**
- Date：2026-10-03。
- Would amend：ADR0002的深度资源预算取值，不改变DEPTH_LIMIT的UNRESOLVED语义。
- Evidence：[job1410810路径审计](../depth-limit-binary-policy-audit.md)。

## 背景

soft policy在component0留下143个depth20未搜索节点，570蛋白、最大27。
binary及改变候选日程的fallback共同增加路径长度；不以降低停止尺寸/ARI解决资源截断。
全部已选分裂仍通过比例.95及ARI.9。旧kway树已政策完成，但不是新拓扑后代的搜索保证。

## 决议

1. **仅下一次component0 H4诊断**将max_depth从20改为42，以job1410810冻结soft树为对照。
   其他配置、SSN、seed、预算、优化轮数及OG/scorer保持不变。
2. 42来自整数门槛 m<=floor(.95*n)：27到1最多22边，原depth20+22。
   无隐藏“深度失败转terminal”、无V1二蛋白人造拆分、无放松child比例/稳定性。
3. 运行使用独立不可变快照。保持depth<20已获准分裂的成员集合与子分区前缀一致；
   numeric ID因DFS插入而改变不构成拓扑不一致。验收整个新后代，不只143原节点。
4. 待处理子树的额外展开上界为854新节点、30,744次Leiden调用；不含重跑和并行回放。
   搜索本身仍24 unique gamma、robust3、每调用10轮；参数版本/cache身份必须记录42。
5. H4验证resolved、唯一完整覆盖、候选预算、准入与显式fallback、调用数/耗时、完整
   serial/parallel一致性；任一真实搜索/无边失败仍UNRESOLVED并阻断OG与评分。
6. **不自动串接H5**：42仅针对component0截断证据，不是全部组件的充分上界。
   H4之后再确定其他组件的预算适用范围及H5执行条件。

## 尚未决定/明确不做

不把42设置为全局默认，不声称准确度改善，不改停止尺寸1或ARI.9，不改比例.95，
不运行官方评分来掩盖未解析节点。生产算法不变；回归wrapper通过显式mode生成depth42配置，仅运行H4。
全局默认仍20，H5需后续单独决议。


## 执行接线

`slurm/h4_h5_hierarchy_regression.sh` 第4参数 `depth42-h4-only`，第3参数必须soft_binary_24_v2。
基线明确使用job1410810，配置从冻结基线复制，仅修改max_depth。默认full旧回归路径保持。
工具depth_budget_regression按成员clade验证depth<20的节点字段、子分区、全部候选trace，
排除numeric ID/parent ID差异；同时检查额外节点854/调用30744上界与无DEPTH_LIMIT。
若仍有其他UNRESOLVED，保留报告并阻止双进程/OG/评分；resolved后做H4完整双进程回放。
H4验收完成即退出，不启动H5或官方评分。Slurm资源保持原56CPU/250G/72小时。
