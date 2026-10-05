# 固定树加权图评分基线（实验）

算法身份：`fixed-tree-weighted-pair-penalty-v1`。
入口：`benchmarks.og_extraction.fixed_tree_graph`；不注册生产 CLI，不修改
`v1_compatible`、生产配置、事件规则、缓存或输出 schema。

## 定义

对当前组件的每个完整子树 v 计算：

```
keep(v) = internal_weight(v) - pair_penalty * n(v)*(n(v)-1)/2
split(v) = sum(best(child))
best(v) = max(keep(v), split(v))
```

`internal_weight` 为两端均位于子树的无向边权之和，每条边只计一次。
缺失边贡献零；节点大小包括所有蛋白，不依赖 OF assigned/unassigned 身份。
`pair_penalty` 必须显式给出、有限且非负，量纲与边权相同，没有隐式归一化。
这是一项预先明确的图目标，不是从审计最优标签拟合的生物学评分。

所有已解决节点均可作为候选，忽略 I/II/III 标签；动态规划选择完整、互不重叠的
树切分，同分保留父节点。结构叶不可再分，即使其 keep 分数为负也保持覆盖。
任意 UNRESOLVED 节点会阻断执行，不因保留大祖先而隐藏失败。

不修改固定树，也不调用 Leiden。边按两端叶的最低公共祖先归桶，后序累计边权；
不保存内部节点完整后代成员表，时间 O(V + E*H)，空间 O(V + P + E)。
没有参考 OG、物种树、基因树或结构域证据输入。

## 使用

组件 JSON 必须只有三个顶层字段：

```json
{
  "nodes": [
    {"cluster_id": 0, "parent_id": null, "n_genes": 2},
    {"cluster_id": 1, "parent_id": 0, "n_genes": 1},
    {"cluster_id": 2, "parent_id": 0, "n_genes": 1}
  ],
  "membership": [[10, 1], [20, 2]],
  "edges": [[10, 20, 2.0]]
}
```

```
python -m benchmarks.og_extraction.fixed_tree_graph \
  --input component.json --pair-penalty 1 --out cut.json
```

这里的 1 仅为两蛋白示例，不是 Embleya 推荐值。输入须来自同一组件：
所有蛋白恰好映射到结构叶，包括孤立蛋白；边端点须存在；无自环/重复边，
使用 u<v 的规范方向；边权有限且非负。存在 n_genes 时核对实际成员计数。
输出含算法身份、显式参数、逐节点分数、所选节点、唯一蛋白归属、根/全叶退化标记，
以及输入、评分代码和 DP 代码 SHA256；已有输出文件拒绝覆盖。

## 实验边界

- 这是网络候选家族基线，**不是 ancestral OG 判定**，也不证明复制发生时间。
- 很低惩罚可能全保留根，高惩罚可能全取结构叶；报告显式标记两种退化，不静默
  修改参数规避。参数依赖边权尺度，不与旧 Leiden gamma 混用。
- 单条很重的桥接边仍可能导致合并；本模型没有连接覆盖广度、结构域或协调证据。
  测试强桥保留父节点是数学目标测试，不是生物学正确性验收。
- 未以 13 基因组 OF 标签挑选参数；尚无全量性能成绩。
- 当前提供单组件 API/JSON CLI。冻结 Parquet/manifest 到该输入的全量适配、
  外部评估与正式 Slurm 批量入口仍待实现；禁止在登录节点运行全量分析。
- 未来评估应先固定参数与相关家族块划分，保留 30637/8863 反例，同时报告
  precision、recall、同/跨物种指标、过合并/拆分及退化情况。
