# SSN 构建与 OrthoFinder 3.1.5 的一致性审计

日期：2026-09-30。OGProfiler HEAD：`f49f2a5`。

## 结论

**当前默认 SSN 与 OrthoFinder 3.1.5 默认全量推断流程并不严格一致。**

搜索架构、覆盖率开关及主要 BH/LRB 判定思想接近；但 NBS 拟合数据选择、
重复 hit 处理顺序、退化拟合、图权重语义及一个并列 RBH 的实现边界存在差异。
之前 S0 的边数接近不等于边集或加权图等价，也不足以把残余差异归于 DIAMOND 版本。

本次不调整 OG、层级或 SSN 生产代码，不部署、不提交 Slurm 作业。

## 审计对象和方法

- OGProfiler 当前默认：`legacy_nbs`、`lrb`、`forward`、coverage disabled。
- 参考源码：远端实际安装的 `orthofinder-3.1.5` 环境。
- OrthoFinder 常规全量运行默认：`gathering_version=(1,0)`、`v2_scores=False`。
  `--assign` 的新物种 clade-specific 阶段会使用另一种 scoring，本报告不声称覆盖该分支。
- 参考函数：`tools/waterfall.py` 的 `scnorm` 和 `WaterfallMethod`；
  `utils/blast_file_processor.py`；`orthogroups/gathering.py::WriteGraph_perSpecies`。
- 从实际源码 AST 提取原始类，执行本地小型差分；调用当前 OGProfiler 函数。
  图权重测试使用与参考 writer 相同的稀疏矩阵加法和逐元素乘法。
- 本地差分数值环境：NumPy 2.1.1、SciPy 1.14.1；退化 curve-fit 的具体小数可能随依赖变化。
- 既有 similarity 单测：11 passed。单测通过并未覆盖以下严格等价要求。

参考源码 SHA256：

| 文件 | SHA256 |
| --- | --- |
| waterfall.py | `304324b1967f917dd523c9d98a7d0c7cffba7cd730ac8b31176ef4f0f80918c8` |
| gathering.py | `d942df722f3d9c94271e357683e080603ee3f958e6f805869d7d6d6cb5a40e3f` |
| blast_file_processor.py | `c7a0aed3a773b2fb377a953c17dfde581a74b299b678f3b174659a2615f73e0c` |
| matrices.py | `0d7452e0c935e6364a3833b1e31393a53344d46f4fa4fdd225c3dcf05ac5520b` |

## 1. NBS 分箱：V1 的切片与 OrthoFinder 不同

OGProfiler `similarity/normalization.py:23-24`：每个 bin 取 `scale+1`，相邻 bin 重叠；
还处理不足一个完整 bin 的尾部。

OrthoFinder `scnorm.GetTopPercentileOfScores`：bin 数使用整除；
每个完整 bin 取 `scale`，不重叠；余数尾部不参与拟合样本选择。

两者使用相同的大体 bin 宽度与 95th percentile，但拟合样本不同。

实测，长度积和 bit score 同为 `1..100`：

```text
OrthoFinder top lengths: 20,40,60,80,100
OGProfiler top lengths: 20,21,40,41,60,61,80,81,100
```

101-hit fixture 中把尾部 score 设为 10000：OrthoFinder 不用尾部拟合；
OGProfiler 把该样本纳入且重复选中，导致拟合及分数明显变化。
对第 1、50、100、101 个 hit，测试结果为：

```text
OrthoFinder: 1, 1, 1, 99.00990
OGProfiler: 209.24492, 0.51473, 0.17754, 17.31144
```

这是受控边界样本，不代表真实数据平均误差；它证明两种算法不等价。
现有 `test_v1_top_bin_logic_keeps_overlapping_95th_percentile_members` 明确锁定的是 V1 行为。

## 2. 最大 HSP 去重发生在不同阶段

OrthoFinder 先通过 `GetBLAST6Scores` 为每个 query-target 取最大 bitscore，再拟合 NBS。

OGProfiler 的 `run_edge_stage` 先将原始 hits 全部送入 NBS，
之后在 `build_retained_edges` 才去重；重复 HSP 已经影响拟合。

fixture：三个独立 hit 的 score 为 40、150、100，另有第一个 hit 的重复 score 10。

```text
第一个最大 hit 的 NBS：OrthoFinder 0.806267；OGProfiler 1.711186
```

影响是否出现在默认 DIAMOND 实验，需要实际确认重复 query-target 数；
它对 multi-HSP 输入和 BLAST 后端是明确的流程差异。

## 3. 相同长度积的退化拟合不等价

OGProfiler 检测自变量方差为零后，使用 `a=0,b=log10(max score)`。
OrthoFinder 对两个以上样本仍调用 `curve_fit`，不是使用最大值归一化。

同长度积、scores=40,80：

```text
OrthoFinder: 约 0.707109,1.414217
OGProfiler: 0.5,1.0
```

非退化、相同拟合样本的 log-linear 模型中，闭式最小二乘与 curve_fit
数学目标相同；那部分主要是数值容差问题。退化分支则是算法选择差异。

## 4. 最大差异：加权图不是同一对象

令 `Cuv` 表示方向 u→v 是否通过阈值，`Buv` 表示完整归一化分数。

OrthoFinder 默认 writer：

```text
connect2 = Cuv + Cvu
Wuv = connect2 * Buv
Wvu = connect2 * Bvu
```

`connect` 是数值型稀疏矩阵，不是布尔逻辑 OR；两方向通过时 `connect2=2`。
writer 为每个方向分别输出，并保留该方向完整 B 分数。
文件写出保留三位小数，交给 MCL。

OGProfiler：先丢弃未通过阈值的方向，再折叠为单一无向边；
`forward` 取低 protein_id→高 protein_id 的通过分数，缺失则取反向。

差分 fixture：

| 场景 | OrthoFinder 图 | OGProfiler 无向边 |
| --- | --- | --- |
| 两方向都通过，Buv=8,Bvu=7 | Wuv=16,Wvu=14 | weight=8 |
| 仅反向通过，Buv=4,Bvu=6 | Wuv=4,Wvu=6 | score_uv=0,score_vu=6,weight=6 |

因此有三种独立偏离：

1. 丢失双向通过的连接倍数；不是所有边统一缩放两倍。
2. 丢失未过阈值但实际存在的前向归一化分数。
3. 用单个无向标量替代两个方向的权重；结果还依赖 canonical ID 顺序。

文档里的“forward 匹配 OrthoFinder3”应理解为不完整的近似说明，不是严格等价。
V1 无向图语义与 OrthoFinder 的 MCL 输入矩阵语义也不应混为一谈。

## 5. BH/LRB 主规则一致，但并列 RBH 有实现边界差异

相同完整 B 下，主要规则相符：

- 按目标物种取 best hit，严格 `score > max-1e-3`。
- 跨物种 reciprocal best hits。
- 查询蛋白的最低 RBH score 用作方向阈值。
- 无 RBH 时使用最佳跨物种 score 加 `1e-6`。
- 同物种非 self hit 使用同一 LRB 阈值。
- retained hit 使用 `>=`。

但实际参考 `GetMostDistant_s` 使用带重复 I 索引的 NumPy 赋值，
同一目标物种存在多个近似并列 RBH 时，该语句并非严格的逐 query minimum reduction。
OGProfiler 使用 Python `min(scores)`，实现了严格的数学最小值。

fixture：query 0 对同一外物种两个目标分别为 7.9995 和 8，都属于 BH；
同物种 hit 0→1 score=7.99975。

```text
参考实现 cutoff=8，排除 0→1。
OGProfiler cutoff=7.9995，保留 0→1。
```

这是参考实现的索引边界问题，不建议为追求一致盲目复制；
应记录为明确的兼容性差异并测量真实数据影响。

## 6. 已对齐的部分及适用范围

| 环节 | 判断 |
| --- | --- |
| 每物种建库、逐有向物种对搜索 | 主架构一致，当前执行 n² searches |
| DIAMOND more-sensitive,evalue=0.001,每对单线程 | 核心参数一致 |
| max-target-seqs/max-hsps=0 时不传参数 | 行为对齐，依赖具体工具版本默认 |
| matrix/gap penalties | 参考显式 BLOSUM62/11/1，V2 使用工具默认，需锁定版本核验 |
| self hits、非正 bitscore | 主要处理相符 |
| 覆盖率与 identity 过滤 | 默认都不作为保边条件 |
| 少于两个拟合样本 | 默认都忽略该组；分箱差异可改变是否落入该分支 |
| BH/LRB 主数学规则 | 大体一致，有上节并列索引例外 |
| NBS、最终边权、图方向性 | 存在实质差异 |

本结论针对默认全量 OrthoFinder 流程，不覆盖所有 search backend、
`v2_scores=True`、homology graph 或增量 assign 分支。

## 7. 原 S0 回归报告的证据与缺口

读取既有 real_embleya 报告：

```text
v2_edges=701245
orthofinder_edges=701880
ssn_edges.v2_vs_orthofinder.equal=false
```

报告已有边集差异、边权差异。例如：

```text
NZ_BIFH01000001.1_1 / NZ_BIFH01000016.1_193
V2=0.7000397860270643；OF graph parser=1.427
```

该 OF 数值来自旧 parser 的无向折叠，不是严格方向对方向比较。

现有回归脚本问题：

1. `parse_orthofinder_graph` 将两个方向折叠到同一个 key，后读方向覆盖前者。
   实测输入 `0→1:16,1→0:14`，结果只剩 `{(0,1):14}`。
2. `weight_diff_count` 在 20 条时提前 break，实际是样例数量而非总数。
3. `only_left/only_right` 也只保留前 20 个，未记录完整差异计数。
4. hit 比对的 `diff_count` 是 key 并集大小，不是真实不一致数量；
   `diffs=[]` 只表示被检查的前 20 个 key 没有差异。
5. 总边数净差 635 不是边集对称差，不应据此计算“边集差异率”。

因此旧报告没有建立“相同 SSN”前提。

## 8. 在修改 OG 前建议完成的验收

1. 固定实际 OrthoFinder 版本、依赖及默认分支。
2. 相同输入 hits 分别走两套 edge 流程，排除搜索版本差异。
3. 按 `(species,original_id)` 对齐身份，比较去重 bitscore、拟合样本和参数、
   完整方向 B、BH、RBH、cutoff、方向 connect。
4. 独立验收无向 support 与加权方向图；明确 V2 Leiden 所需无向投影规则。
5. 记录完整边集对称差、权重误差分布、方向缺失、双向通过率及组件规模，
   样例限制与总体计数分开。
6. 将 SSN 差异影响与层级验收问题分别测试；目前未量化这些差异对
   Open_Orthobench precision 崩溃的独立贡献。

## 本地证据位置

临时源码快照与差分 harness：
`/private/tmp/ogp-ssn-audit-of315/`。

- `audit.py`：八项差异断言，调用当前 OGProfiler 与参考源码类。
- `results.json`：完整小型差分输出。
- `s0-report.json`：既有 real_embleya 报告副本。
- `src/orthofinder/`：实际参考源码快照。

在仓库目录运行：

```bash
PYTHONPATH=src python /private/tmp/ogp-ssn-audit-of315/audit.py
```

该命令使用本地已装 SciPy 的 Python，不需要远端运行或大规模计算。
