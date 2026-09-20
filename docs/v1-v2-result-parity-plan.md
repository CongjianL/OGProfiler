# V1 ↔ V2 结果一致性改动方案

> 状态：草案（待评审）
> 目标读者：重构与回归测试负责同学
> 前置结论：`src/ogprofiler/`（OGProfiler 2）不是 `OGProfiler.py`（V1）的等价移植，而是算法级重写。
> 本文把导致结果不一致的 8 类差异整理成可执行的改动方案，并明确每一处是「配置对齐」还是「算法回迁」。

---

## 0. 背景与目标

经过 `OGProfiler.py`（V1，1832 行）与 `src/ogprofiler/`（V2）逐函数对比，结果不一致由 8 类差异造成，按影响从大到小：

1. 边构建阶段新增了覆盖率过滤（V1 从不使用覆盖率）。
2. 同源搜索（search）阶段的参数不一致，导致 hit 集合本身不同。
3. 输入准备（prepare）阶段的文件/基因排序与序列规范化差异。
4. 社区检测 / 层级构建算法整体重写。
5. 演化事件标注算法整体重写。
6. OG 提取与导出算法整体重写。
7. 边权重语义从「单向 NBS」变成「双向 max 对称化」。
8. 若干小差异（`ar` 方向 bug、单 hit 退化组 NBS、Leiden 迭代次数、blastp 多 HSP 去重）。

按差异性质分两类处理（已决议）：

- **对齐项（第 1/2/7/8 项）**：关闭/对齐 V2 中可配置的行为差异（覆盖率过滤、搜索参数、边权重、NBS 细节），使上游产物（hits、SSN 边）与 V1 一致。
- **保持项（第 3/4/5/6 项）**：保留 V2 的输入准备（确定性排序、序列规范化）与「层级 → 事件 → OG」算法链，作为**有意差异/算法替换**在文档与 `KNOWN_ISSUES.md` 中公示，不做 V1 还原。

即：真正动代码的只有第 1/2/7/8 项；输入准备与层级/事件/OG 差异只公示、不还原。

### 改动原则（贯穿全文）

1. **不破坏 V2 生态**：不引入 V1 的过程式副本（如 `legacy_v1.py`），不做「按 mode 分流到 legacy 过程式代码」的改法；所有改动落在 V2 现有的 stage / parquet / manifest 架构内。
2. **还原语义、优化实现**：对齐目标是 V1 的**算法语义**（判定规则、阈值、顺序），不是 V1 的代码；实现允许用更高级、更成熟的包或更高效的写法（如集合位图、向量化、更优的 gamma 搜索）。
3. **影响大的优先**：先处理覆盖率过滤、搜索参数、边权重等对结果影响大的项，层级/事件/OG 等算法项在 S5 决策后按序推进。

---

## 1. 边构建：覆盖率过滤

### 现状
`src/ogprofiler/similarity/engine.py` 的 `filter_coverage()` 在 NBS 归一化后执行，默认丢弃 query/target 覆盖率 < 50% 的 hit。V1 在 `read_blast_results()` 中只计算覆盖率、从不使用。

### 依据（OrthoFinder3）
对照 OrthoFinder3（`scripts_of/gathering.py` 的 WaterfallMethod、`scripts_of/blast_file_processor.py`）：从 BLAST 结果到最终图**全程不按覆盖率/identity 过滤**，只取 bit score（并排除同种 self-hit）。V1 当年即按此实现，V2 的 50% 覆盖率过滤是额外引入的行为，应移除。

### 最终设定
**不按覆盖率过滤**（与 OrthoFinder3 / V1 一致）：

1. 默认值：`min_query_coverage: 0`、`min_target_coverage: 0`、`min_bidirectional_coverage: 0`（`src/ogprofiler/config.py` 的 `DEFAULT_CONFIG["edges"]`）。
2. 保留过滤器能力，但增加显式开关 `edges.apply_coverage_filter: false`，`false` 时直接跳过 `filter_coverage()`。
3. `edge-manifest.json` 的 `parameters` 记录 `apply_coverage_filter`。
4. 更新 `docs/edge-engine.md`，删除「默认 50% 是 V2 有意行为」的表述。

### 涉及文件
- `src/ogprofiler/config.py`
- `src/ogprofiler/similarity/engine.py`（`EdgeBuildConfig`、`build_retained_edges`）
- `src/ogprofiler/similarity/stage.py`（`run_edge_stage` 的参数与 manifest）
- `docs/edge-engine.md`

### 验收
- 用同一份 `hits.parquet`，`apply_coverage_filter=false` 与 V1 的 retained hit 数一致（在 lrb 且无其他差异的前提下）。
- 与 OrthoFinder3 比对：同一数据集下，V2 的 retained edge 集合（去重后无向边）与 OrthoFinder3 的 graph 边集（`connect2` 连通性）一致。
- manifest 中 `coverage_filtered_hits == directional_hits`。

---

## 2. 同源搜索：参数对齐

### 现状
V1、V2 与 OrthoFinder3 的搜索命令对比：

| 参数 | V1 (`OGProfiler.py`) | V2 (`src/ogprofiler/`) | OrthoFinder3 (`scripts_of/config.json`) |
|---|---|---|---|
| DIAMOND 灵敏度 | `--more-sensitive` | 默认 `sensitive` | `--more-sensitive` |
| `--max-target-seqs` | 未设置 | `0`（不传） | 未设置 |
| `--max-hsps` | 未设置 | `0`（不传） | 未设置 |
| evalue | `0.001` | `0.001` | `-e 0.001` |
| 输出字段 | `-f 6`（12 列） | 显式 8 列 | `-f 6`（12 列，读第 12 列 bit score） |
| MMseqs 灵敏度 | `-s 7.5` | `sensitive` → `-s 5.7` | 未设置（mmseqs 默认 5.7） |

### 依据（OrthoFinder3）
- DIAMOND：`diamond blastp ... --more-sensitive -p 1 --quiet -e 0.001 --compress 1`（无 `--max-target-seqs`、无 `--max-hsps`）。
- blastp：`blastp -outfmt 6 -evalue 0.001 -query ... -db ...`（无 `-max_target_seqs`、无 `-max_hsps`）。
- mmseqs：`mmseqs search ... --threads 1`（未设 `-s`，即 mmseqs 默认 5.7）。

### 最终设定
1. `DEFAULT_CONFIG["search"]`：
   - `sensitivity: more-sensitive`（对齐 OrthoFinder3 DIAMOND）
   - `max_target_seqs: 0`（不传该参数，用各后端默认值，对齐 OrthoFinder3/V1 不传）
   - `evalue: 0.001`（不变）
2. 新增 `search.max_hsps: 0`（`0` = 不传，对齐 OrthoFinder3/V1；非 0 才传）。MMseqs/BLAST+ 同理。
3. MMseqs：保持 `sensitivity: sensitive`（→ `-s 5.7`），与 OrthoFinder3 的 mmseqs 默认一致；**不采纳 V1 的 `-s 7.5`**。
4. **搜索架构改为逐物种对**（`search/stage.py` 重写）：每个物种建一个 DB，逐对（n_species² 次）diamond，线程池并行、每对单线程，合并输出。这是对齐 V1/OrthoFinder 的关键——单库全对全会因 DIAMOND evalue 依赖库大小而得到不同 hit 集（实测差 2~8 倍）。
5. 输出字段：V2 保持 8 列内部格式；第 1 项关闭覆盖率过滤后此差异对结果无影响。在 `parse_tabular_hits()` 中注明 8 列与 12 列的对应关系。
6. `search-manifest.json` 记录 `sensitivity`、`max_target_seqs`、`max_hsps`、`evalue`、`search_mode=per_species_pair`、`species_ids`、`n_pairs`。

### 涉及文件
- `src/ogprofiler/config.py`
- `src/ogprofiler/search/base.py`（`SearchParameters` 增加 `max_hsps` 字段）
- `src/ogprofiler/search/diamond.py`
- `src/ogprofiler/search/mmseqs.py`
- `src/ogprofiler/search/blast.py`
- `src/ogprofiler/search/hits.py`
- `src/ogprofiler/search/stage.py`
- `docs/search-backend.md`

### 验收
- 同一数据集，V2 `search` 产出的 raw `hits.tsv` 与 V1 对应 `BlastResults/*.out` 行数/内容一致（按 genome-pair 合并比对）。
- 与 OrthoFinder3 比对：V2 `hits.parquet`（按 genome-pair 还原）与 OrthoFinder3 同参数跑出的 `Blast*.txt` 行数、bit score 一致。
- DIAMOND：验证 `--more-sensitive`（不传 `--max-target-seqs`）与 V1/OrthoFinder 未传参数时的输出一致。
- MMseqs：验证 `-s 5.7` 与 OrthoFinder3 默认一致。
- S0 实测（13 物种真实蛋白组，128K 蛋白）：V2 10,941,213 vs V1 10,941,049 vs OrthoFinder 10,941,104（误差 < 0.002%）。

---

## 3. 输入准备：保持 V2（不退回到 V1）

### 现状
| 项目 | V1 | V2 |
|---|---|---|
| 基因组排序 | `os.listdir()`（文件系统顺序） | `sorted_input_paths()`（UTF-8 字节序） |
| 基因排序 | FASTA 文件内出现顺序 | `sorted_original_ids()` |
| 序列处理 | 原样保留 | 大写化 + 非法字符校验（默认 `error`） |
| header 切分 | `split(' ')[0]` | `split(maxsplit=1)[0]` |
| 基因组数下限 | ≥ 4，否则报错退出 | 无限制 |

### 决策
**保持 V2 现状，不提供 V1 兼容开关，也不回退。** 理由：

1. V2 的 `sorted_input_paths()` / `sorted_original_ids()` 是确定性排序，优于 V1 的 `os.listdir()` 文件系统顺序（非确定性、跨平台不可复现）。
2. V2 的大写化 + 非法字符校验更严格、更可诊断；V1 原样保留序列存在隐患。
3. `split(maxsplit=1)` 对含 tab 的 header 更稳健；V1 的 `split(' ')` 是历史遗留。

因此该项**不作为 parity 对齐项**，列为「有意差异」，在文档与 `KNOWN_ISSUES.md` 公示。

### 改动方案（仅文档/回归脚本，无源码改动）
1. 源码保持 V2 现状，不新增 `input.order` / `gene_order` / `preserve` / `min_species` 等开关。
2. 在 `KNOWN_ISSUES.md` 与 `docs/` 相关文档中，把「输入顺序/规范化差异」列为有意差异，并写明理由。
3. 回归比对脚本（S0）按 **原始基因 ID** 做 canonical 重映射后再比较 V1/V2 产物，而不是按 `gN` / `protein_id` 编号直接比较，从而消除编号差异的干扰。

### 影响与注意
- 由于保持 V2 排序，V1 的 `G0/G1/...` 与 V2 的 `species_id` 编号不再一一对应。
- S4（边权重）的「与 V1 一致」验收必须通过原始 ID 重映射后进行；层级/事件/OG 已决议保持 V2，不受编号差异影响。

### 涉及文件
- 无源码改动
- `KNOWN_ISSUES.md`
- `docs/` 相关文档（如 `docs/OGProfiler2_Design_Specification.md`）
- S0 回归比对脚本（`dev/` 或 `tests/`）

### 验收
- 文档明确列出该差异为「有意差异」及理由。
- S0 回归脚本能按原始基因 ID 对齐比较 V1/V2 产物，不因编号顺序差异误报。

---

## 4. 社区检测 / 层级构建：保持 V2 算法（不用 V1 算法）

### 现状
V2（`src/ogprofiler/hierarchy/engine.py` + `resolution.py`）是 DFS + 稳定性分辨率搜索（固定 seed、多 seed ARI/NMI、`min_child_size` / `max_child_fraction` / `max_depth` 验收）；V1（`HHN.RunCommunityDetection` / `BipartiteGraphs` / `SpecificBoard`）是 `la.find_partition` + 二分/递增 gamma 搜目标社区数，`n_iterations=10`，无固定 seed。

### 决策
**保持 V2 算法，不用 V1 算法。** 理由：

1. 层级算法是 V2 重写的核心科学价值（可复现、有基准设施），还原 V1 会破坏 V2 生态。
2. V1 的层级算法无固定 seed、靠启发式目标社区数，可复现性与科学性弱于 V2。
3. 这是**有意算法替换**，不是 bug 或遗漏，不应回退。

### 改动方案（仅公示，无源码改动）
1. 不改 `src/ogprofiler/hierarchy/` 的算法。
2. 在 `KNOWN_ISSUES.md` 与 `docs/` 相关文档中，把「层级算法为有意替换」作为已知差异公示，写明 V1/V2 的算法差异与理由。

### 涉及文件
- 无源码改动
- `KNOWN_ISSUES.md`
- `docs/OGProfiler2_Design_Specification.md`（V1 regression 章节）

### 验收
- 文档明确记录「层级算法有意替换、不还原 V1」及理由。

---

## 5. 演化事件标注：保持 V2 算法（不用 V1 算法）

### 现状
V2（`src/ogprofiler/evolution/network.py`）在层级树上用 species bitmap + overlap score 标注；V1（`HHN.GetEvolutionEvents`）在 hhn 图上按顶点 degree + 邻居 `genomeIDs` 集合标注 `I/II/III-1/III-2/III-3`。

### 决策
**保持 V2 算法。** V1 的事件标注依赖 V1 的 hhn 图结构（degree + 邻居），与 V2 的层级树形态不同，无法在不破坏 V2 的前提下单独还原；且第 4 项已决议保持 V2 层级，故本项随之保持 V2。

### 改动方案（仅公示，无源码改动）
1. 不改 `src/ogprofiler/evolution/` 的算法。
2. 在 `KNOWN_ISSUES.md` 与 `docs/` 中公示事件标注语义差异（`I/II/III-*` ↔ `SPECIATION_LIKE/DUPLICATION_LIKE/MIXED/...`）。

### 涉及文件
- 无源码改动
- `KNOWN_ISSUES.md`
- `docs/network-evolution-annotation.md`

### 验收
- 文档明确记录事件标注差异及理由。

---

## 6. OG 提取与导出：保持 V2 算法（不用 V1 算法）

### 现状
V2（`src/ogprofiler/output/results.py`）从层级 terminal family 直接导出；V1（`ExtractOGSorted` / `ExtractOG` / `WriteOGFiles` / `GetGenesIDs`）从基因组数高到低按 `Event + genomesNum` 选节点、递归吞并子节点，输出 `OGFile_coalescence_SameGenome.txt`。

### 决策
**保持 V2 算法。** OG 提取依赖层级与事件语义（第 4/5 项已决议保持 V2），故本项随之保持 V2。

### 改动方案（仅公示，无源码改动）
1. 不改 `src/ogprofiler/output/` / `src/ogprofiler/orthology/` 的算法。
2. 在 `KNOWN_ISSUES.md` 与 `docs/` 中公示 OG 导出差异（terminal family 导出 ↔ V1 coalescence 导出）。
3. 若下游仍需 V1 的 `OGFile_coalescence_SameGenome.txt` 输出**格式**，作为独立格式适配项另议，不属于算法还原。

### 涉及文件
- 无源码改动
- `KNOWN_ISSUES.md`
- `docs/final-result-export.md`

### 验收
- 文档明确记录 OG 导出差异及理由。

---

## 7. 边权重语义：对齐 OrthoFinder3 的「对称连通性 × 前向分数」

### 现状
- V1：SSN 边权重 = `get_connections` 写入的 `NBS`（低编号基因组→高编号方向，缺失才回退反向）。
- V2：`_symmetrize` 默认 `max(score_uv, score_vu)`。
- OrthoFinder3（`scripts_of/gathering.py` 的 `WriteGraph_perSpecies`）：边存在 = 对称连通性 `connect2 = connect(i,j) + connect(j,i).T`；边权重 = 前向分数 `B_connect = connect2 × B[i→j]`。

### 影响
即使 edges 相同，Leiden 权重也不同（max vs 前向分数），分区结果不同。

### 最终设定（对齐 OrthoFinder3）
1. 边存在：任一方向通过 LRB/RBH 阈值即保留（对称连通性），V2 现有 union 逻辑不变。
2. `EdgeBuildConfig.symmetrization` 新增 `forward` 并设为默认：
   - 双向保留：`weight = score_uv`（前向分数，对齐 OrthoFinder 的 `B[i→j]`）。
   - 仅反向保留：`weight = score_vu`（V1 的 NBS 回退语义）；OrthoFinder 原实现此时取前向分数 `B[u→v]`（可能低于阈值），如需严格一致，edge 阶段应额外保留「反向-only 边的前向归一化分数」，差异先记入 `KNOWN_ISSUES.md`。
3. `max/min/mean/geometric_mean` 保留为可选；`RetainedEdge` 始终保留 `score_uv/score_vu` 便于追溯。
4. 注意：`forward` 依赖「u=低编号、v=高编号」的 canonical 顺序，验收时按原始基因 ID 重映射后比较（第 3 项已决议保持 V2 编号）。

### 涉及文件
- `src/ogprofiler/similarity/engine.py`
- `src/ogprofiler/config.py`
- `docs/edge-engine.md`

### 验收
- 同一 `hits.parquet`，`forward` 权重与 V1 `ssn.gml` 的 `NBS` 边属性逐一一致（按原始基因 ID 重映射后比较）。
- 与 OrthoFinder3 比对：V2 的 `score_uv` / `score_vu` 与 OrthoFinder3 的 `B[i→j]` / `B[j→i]` 逐值一致（同 hits、同 legacy_nbs 归一化）。

---

## 8. 小差异清单（一并修复或记录）

| # | 差异 | V1 | V2 | 方案 |
|---|---|---|---|---|
| 8.1 | `-d ar` 方向 bug | `sorted_connected_rh` 误加载同一 `(I,J)` 矩阵，条件恒真 → 保留全部 hit | 正确反向 join | 保留 V2 修复，但写入 `KNOWN_ISSUES.md` 说明这是**有意修复**，会导致 `ar` 结果不同 |
| 8.2 | 单 hit 退化组 NBS | 全 0 矩阵 | `a=0, b=log10(max)`，组最大映射为 1 | 新增 `similarity.nbs_fallback: v1_zero | v2_max`，默认 `v1_zero`（对齐 V1 与 OrthoFinder3：`NormaliseScores` 对「Too few hits」返回全 0） |
| 8.3 | Leiden 迭代次数 | `n_iterations=10` | leidenalg 默认（2） | 随第 4 项「保持 V2 算法」作为有意差异公示，不再单独处理 |
| 8.4 | blastp 多 HSP 去重 | 矩阵覆盖，保留最后一行（V1 的历史行为） | `_deduplicate` 保留 max | 第 2 项关闭 `--max-hsps` 后，**保持 V2 的 max 去重**（对齐 OrthoFinder3 的 `if score > B[i,j]`），不复现 V1 的 last-wins |

### 涉及文件
- `src/ogprofiler/similarity/normalization.py`（8.2）
- `src/ogprofiler/similarity/engine.py`（8.4）
- `src/ogprofiler/config.py`
- `KNOWN_ISSUES.md`

---

## 9. 建议实施顺序

按依赖关系与影响面排序：

1. **第 1 项**（覆盖率过滤）——独立、影响最大、纯配置。
2. **第 2 项**（搜索参数）——决定 hit 集合，先对齐才能做后续回归。
3. **第 3 项**（输入顺序/规范化）——保持 V2，仅文档公示 + 回归脚本按原始 ID 重映射。
4. **第 7 项 + 第 8 项**（边权重与 NBS 细节）——决定 SSN 边属性。
5. **第 4/5/6 项**（层级→事件→OG 链）——已决议保持 V2 算法，作为有意差异写入 `KNOWN_ISSUES.md` 与相关文档公示（无源码改动、无后续依赖）。

---

## 10. 回归验证与可追溯性

1. 新增「V1 ↔ V2 ↔ OrthoFinder3」三路回归测试集，用固定小数据集分别跑三个实现，逐产物比对：
   - `BlastResults/*.out` ↔ `search/hits.tsv` ↔ OrthoFinder3 `Blast*.txt`
   - `ssn.gml` ↔ `edges/retained_edges.parquet` ↔ OrthoFinder3 graph（MCL matrix / `connect2`）
   - （层级/事件/OG 已决议保持 V2，不与 V1/OrthoFinder3 逐节点比对，只公示差异）
2. 所有可配置差异（coverage、search、symmetrization、nbs_fallback）必须写入各自 stage 的 manifest `parameters`，保证缓存命中与差异可归因。
3. 在 `CHANGELOG.md` 中记录每一项默认值的变更及理由。
4. 文档同步：更新 `docs/edge-engine.md`、`docs/search-backend.md`、`docs/OGProfiler2_Design_Specification.md` 第 29.3 节，明确「V1/OrthoFinder3 兼容默认值」与「有意算法替换」的边界。

---

## 11. 需评审的开放问题

1. ~~目标是「严格复现 V1」还是「尽量对齐 + 可归因差异」？~~ **已决议**：第 1/2/7/8 项对齐 V1；第 4/5/6 项保持 V2 算法、只公示差异。
2. 第 2 项 `max_target_seqs` 默认是否真的改回 25？V1 的 25 是 DIAMOND 隐式默认，若 V1 实际运行环境 DIAMOND 版本不同，该值可能不同——需要现场用 `diamond --help` 确认并记录版本。
3. ~~第 3 项是否值得为「文件系统顺序」引入非确定性？~~ **已决议**：保持 V2 的 `sorted` 确定性顺序，不引入 `v1_filesystem`，输入顺序/规范化差异作为有意差异公示。
4. 第 8.1 项（`ar` bug）是修复还是复现？建议修复并作为已知差异公示，而非复现 V1 的 bug。

---

## 12. PR / Issue 拆分索引

已按「tracer-bullet 垂直切片」拆成可独立合入的 PR，并发布到 issue tracker（`CongjianL/OGProfiler`）。依赖关系与第 9 节实施顺序一致。原 9 个切片中的 S7/S8 已随「保持 V2」决议并入 S6（见下表）。

| Issue | 标题 | 类型 | 阻塞于 |
|---|---|---|---|
| #2 | S0: 建立 V1↔V2 回归比对脚手架 | AFK | 无 |
| #3 | S1: 覆盖率过滤对齐 V1（默认关闭 + 显式开关） | AFK | #2 |
| #4 | S2: 搜索参数对齐 V1 | AFK | #2 |
| #5 | S3: 输入准备保持 V2（有意差异公示 + 回归脚本按原始 ID 重映射） | AFK | #2 |
| #6 | S4: 边权重 forward（前向分数）+ NBS 细节对齐 | AFK | #3, #5 |
| #7 | S5: 层级/事件/OG 链保持 V2（已决议，公示差异） | HITL | #3, #4, #5, #6 |
| #8 | S6: 层级/事件/OG 差异公示（合并原 S6/S7/S8） | AFK | #7 |

> 原 #9（S7）、#10（S8）随本决议并入 #8，后续在 tracker 上关闭。

标签：AFK 切片使用 `ready-for-agent`；#7（决策）保留 `ready-for-human`。
