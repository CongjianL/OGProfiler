# V2 SSN 对齐 OF3.1.5 与 Leiden mean 投影

依据固定 hits 审计作业 1410698，修改默认完整运行 SSN；OG 和层级算法保持原样。

## 修复契约

1. NBS 前按方向 query-target 取最大 bitscore。拟合输入按 query/target
   protein_id 排序，复现共享 ID 顺序下 OF 的稀疏矩阵遍历，包括长度积相同时的顺序。
2. NBS 使用完整、不重叠 bin，舍弃不足一个 bin 的尾部；95th percentile
   使用 NumPy。小于 100 hits 的分支保持全部样本。
3. 参数通过 SciPy curve_fit 求解，包括等长度积的退化问题；不再用最大 score
   替代该退化拟合。默认少于两个拟合点仍输出零矩阵；显式 v2_max 选项继续保留。
4. BH/RBH 规则保持一致；LRB cutoff 精确复现 OF 的重复行索引赋值，
   每个方向物种对取最后一个 RBH target 的值，再跨目标物种取 min，
   保留 OF 的 sentinel 与无 RBH fallback。
5. 方向 connect 判定后，从完整 B 组装 `Wuv=(Cuv+Cvu)*Buv`，
   而非丢弃未过阈值的反向 B。双向 connect 乘数为 2。
6. 默认 Leiden 权重：`Uuv=(Wuv+Wvu)/2`，不是 forward+reverse 的和。
   真正缺失的方向为 0；该方向存在但未过自己的阈值时仍使用完整 B。

例：双向 connect、B=(8,8) → W=(16,16) → U=16。
仅正向 connect、完整 B=(9,7) → W=(9,7) → U=8。
仅正向 connect、反向 hit 缺失、B=(9,0) → W=(9,0) → U=4.5。

Leiden 保持 `directed=False`。W 使用完整浮点精度后做 mean，不先模拟 OF writer
的三位小数字符串量化；这是明确的 Leiden 投影契约，区别于实际 MCL 文件。
`score_uv/score_vu` 在 LRB retained_edges 中存储组装后的 W；normalized_hits 存储 B。
显式 forward/max/min 等配置仍可选择 W 的其他投影。
RBH/AR/ARB 等非默认实验方法保持其原有 retained-B 组装行为。

## 缓存与迁移

- `EDGE_ALGORITHM_VERSION` 升为 `of315-nbs-lrb-directional-mean-v4`，旧 edge manifest
  即使输入和显式配置相同也需重建；修复前的 SSN 缓存将失效。
- 默认配置与 EdgeBuildConfig 均改为 mean。已有配置显式设 forward 时仍以配置为准，
  新运行请检查 `edges.symmetrization=mean`。
- 下游 component/hierarchy 必须对应新 retained_edges checksum，避免混用旧图结果。
- `legacy_nbs` 配置名为兼容 CLI 保留，但现在指 OF 对齐模型，不再代表冻结 V1 分箱。
- 新增运行依赖 `scipy>=1.11`；复现需记录 NumPy/SciPy 版本及参考源码 hash。

## 验证

- 边界测试覆盖完整/尾 bin、最大 HSP 去重、等长度积、并列 RBH cutoff、
  双向连接乘数、未通过阈值的反方向分数保留。
- 参考差分使用实际 OF3.1.5 原始函数，不以生产实现自行生成 expected。
  合成数据包括长度积 ties、乱序输入、重复 HSP、同物种 hits、方向不对称。
- 预期 max bitscore、拟合样本/参数、B、BH/RBH、cutoff、connect、方向 W
  在容差内一致；V2 图等于 OF W 的 mean 投影。
- 审计报告中的方向 W 与无向 U 直接比较仍可能不同，这是预期的图对象区别。
- 修复后全量固定 hits 重跑和 OG 准确度验证尚待执行；本次不自动提交大规模作业。

已观察的验证结果：本地完整测试集 96 passed；远端小型差分与 CLI 回归 30 passed；
Ruff、similarity 模块 Mypy、pip check 均通过。少量两点/退化拟合产生 OF 同类
OptimizeWarning（协方差估计警告），差分仍在容差内通过。
本地使用 NumPy 1.26.4 / SciPy 1.17.1；远端使用现有环境的 NumPy 2.5.2 / SciPy 1.18.0。
