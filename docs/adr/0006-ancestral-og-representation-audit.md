# ADR0006：先审计固定层级表示能力，再设计独立祖先OG策略

- 日期：2026-10-05。
- 状态：固定层级诊断和实验树切分基础已实施；生物学评分与生产接线为Proposed。
- 来源：用户提供方案并要求建立PR继续优化。
- 范围：Embleya 13基因组，冻结job1411104/current-soft42与同次OF参考文件。

## 决议与隔离

保留v1_compatible及其I/II/III语义，生产src、默认配置、缓存与输出schema保持原状。
当前网络层级不是协调基因树；本次不把多叉等同复制，不把同物种多拷贝自动排除在候选家族外。
失败层级继续拒绝正常选组，即使有可覆盖它的大祖先节点也不掩盖UNRESOLVED。

新增benchmarks实验基础：
1. 固定全层级/合格层级/实际OG三层每参考组最佳成员F1。
2. 完整、不相交树切分的全局pair F1 oracle；仅为参考标签驱动诊断。
3. 接受外部keep/split分数的无标签动态规划接口；不自创生物学分数，不开放生产CLI策略。

这不修订ADR0002的OG语义保持约束；未来ancestral_og应以新策略、独立算法身份进入，
不得静默改变v1_compatible。祖先层级、物种树输入、定根/协调不确定性也须进入该身份。

## 审计口径

- primary为OF assigned蛋白，OF unassigned不当主指标负例；输出另保存节点完整成员数。
- 每组最佳节点采用投影成员F1=2*intersection/(reference_size+node_assigned_size)，
  并列按更小节点、稳定字符串ID确定。它不是可部署策略；独立最优节点可互相嵌套。
- 合格节点依照实际V1事件函数：I且n_species>1，或None；显式SSN singleton组件单列允许。
- 输出continuous eligibility_gap/selection_gap及exact匹配计数，不凭任意“高质量阈值”分型。
  两种gap不宣称可相加为丢失蛋白对，也不强制非负；实际组可受active-view消耗影响。
- 完整切分oracle最大化全森林的2TP/(PP+R)，R是冻结参考总同组对。
  用整数分式迭代求全局最优比率；每轮树DP解可分的2TP-q*PP目标，收敛以整数相等判定。
  设置迭代上限，未收敛则失败，不输出伪上界。
- oracle只界定固定树可行完整切分的pair F1，不是macro F1/B-cubed各自最大值，
  也不一定覆盖现有active-view产生的所有非完整子树组。
- 每个蛋白唯一覆盖，组间不嵌套；同/跨物种同组pair单列，不能称为ortholog准确率。
- 节点成员、拓扑、species counts、输入/产物manifest/checksum与结束后完整性逐层核对。

## 性能与资源

计数按组件postorder聚合参考ID频次，不存所有内部节点的完整后代蛋白列表。
全局DP保留紧凑拓扑/TP/PP和terminal映射；不枚举128k蛋白全体两两组合。
全组件文件检查和审计放Slurm，2CPU/8GB/2h，无array，无Leiden/搜索/SSN构建。
资源时间限额高于六子树实验的1h是因为文件/校验范围扩至全部组件，不改科学参数。

## 后续里程碑

- [x] 完成全量表示审计，确定优先独立选择策略；因果贡献仍保持限定。见[审计结果](../embleya-representation-audit-results.md)。
- [ ] 按相关家族块制定开发/保留集，保留30637/8863反例；历史panel标准不改。
- [ ] 核对实际OF3.1.5 Orthogroups.tsv的写出阶段，冻结基准不替换。
- [ ] 设计标签无关keep/split分数，显式目标祖先与证据缺失状态。
- [ ] 用祖先前后复制、丢失、不平衡、多叉、结构域桥接的演化fixture验收生物学模型。
- [ ] 单独接入物种树/协调证据，处理细菌HGT与定根不确定性。
- [ ] 通过保留验证后才考虑ancestral_og生产配置及serial/parallel/cache/schema验收。

本次fixture只证明DP数学、覆盖、门禁与计分，不声称已实现这些生物学判定。
本PR不采用job1411157失败的统一回退细化政策，也不使用OpenBench成绩选择方法。

## 2026-10-05 增量：资格受限完整切分对照

`pair_f1_oracle(..., restrict_eligible=True)` 使用每个组件显式资格集合；
默认仍为全节点 oracle。`run_audit` 同时输出两者及计数、覆盖验证。
这只收紧诊断可行域，不放开生产事件规则，也不引入参考标签到推理路径。
全量数值待下一轮审计；见结果文档中的新增字段与解释限制。
