# OGProfiler 2 分阶段详细开发任务

> 文档定位：把 OGProfiler 2 总体设计转化为可执行开发计划。每个阶段包含目标、任务、交付物、测试要求、验收标准与进入下一阶段的条件。

---

## 0. 总体开发策略

OGProfiler 2 不采用“从旧脚本第一行开始逐段改写”的方式，而采用：

> **冻结 V1 + 新建 V2 + 垂直切片 + Regression 驱动开发。**

总优先级：

```text
Hierarchy Engine
    >
Edge Engine
    >
Scientific Validation
    >
Phylogenetic Refinement
    >
UI / Export
```

推荐版本节奏：

```text
0.x  原型和性能验证
1.x  可重复运行的内部开发版
2.0  Network hierarchy MVP
2.1  Orthology + phylogenetic refinement
2.2+ Extreme-scale / distributed extensions
```

---

# Phase 0：冻结 OGProfiler 1 与建立基线

## 目标

把当前 `OGProfiler.py` 固定为 reference implementation，建立后续所有重构和算法修改的比较基础。

## 任务 0.1：Legacy 归档

- [ ] 将当前版本复制到 `legacy/OGProfiler_v1.py`
- [ ] 记录 Python 版本
- [ ] 记录主要依赖版本
- [ ] 记录 DIAMOND/MMseqs/BLAST/FastTree/MAFFT 等外部工具版本
- [ ] 保存一份完整运行命令模板
- [ ] 添加 `legacy/README.md`
- [ ] 明确 V1 后续只允许 critical bug annotation，不继续功能开发

## 任务 0.2：整理已知问题清单

至少记录：

- [ ] 单脚本高耦合
- [ ] hierarchy 期间重复 `subgraph()`
- [ ] resolution search 中反复 `partition.subgraphs()`
- [ ] multiprocessing 反复创建 Pool/Manager.Queue
- [ ] graph 在 worker 之间重复 pickle
- [ ] internal hierarchy node 重复保存完整 geneIDs
- [ ] pairwise ortholog combinations + shortest paths 高复杂度
- [ ] `sorted_connected_rbh/rh` 中 CM dump 缩进问题
- [ ] reciprocal-hit 方向加载/索引问题
- [ ] undirected edge weight 依赖输入方向的问题
- [ ] `os.listdir()` 非确定性顺序
- [ ] 随机 seed 未统一
- [ ] hard-coded FastTree path
- [ ] `eval()` + `methodcaller()` pipeline
- [ ] `shell=True` 外部命令调用
- [ ] coverage 已计算但未正式进入过滤逻辑

## 任务 0.3：建立 baseline datasets

### Dataset A：Small sanity dataset

- [ ] 4 个小 proteomes
- [ ] 每个 proteome 10–100 个蛋白
- [ ] 存在简单 one-to-one homologs
- [ ] 用于完整 pipeline smoke test

### Dataset B：Paralog dataset

- [ ] 明确存在 recent duplication
- [ ] 存在 RBH 与非 RBH homologs
- [ ] 用于测试 BH/RBH/LRB 行为

### Dataset C：Gene family expansion

- [ ] 一个大型扩增 family
- [ ] 多个层级 subdivision
- [ ] 用于测试 hierarchy depth 和 resolution search

### Dataset D：Fusion / multidomain

- [ ] 人工包含 fusion protein
- [ ] 包含局部 domain-only hits
- [ ] 用于测试 coverage 和 bridge sensitivity

### Dataset E：Large connected component

- [ ] 构造一个明显大于常规 family 的 CC
- [ ] 用于 runtime / peak RSS benchmark

## 任务 0.4：记录 V1 baseline

对 A–E 数据集保存：

- [ ] SSN node count
- [ ] SSN edge count
- [ ] connected component count
- [ ] hierarchy node count
- [ ] hierarchy depth
- [ ] terminal family count
- [ ] family membership
- [ ] ortholog pair count
- [ ] runtime
- [ ] peak RSS
- [ ] disk usage

## 交付物

```text
legacy/
  OGProfiler_v1.py
  README.md
benchmarks/
  datasets/
  baseline_v1/
KNOWN_ISSUES.md
```

## 验收标准

- V1 可以在至少 Dataset A/B/C 上稳定重现结果。
- 所有 baseline 结果可版本控制或通过 manifest 校验。
- 后续任何 V2 结果变化都可以和 baseline 对照。

---

# Phase 1：项目骨架、配置与核心数据模型

## 目标

建立干净、可测试、可复现的 OGProfiler 2 工程基础。

## 任务 1.1：项目初始化

- [ ] 创建 `pyproject.toml`
- [ ] 创建 `src/ogprofiler/`
- [ ] 创建 `tests/`
- [ ] 建立基础 lint/type/test 配置
- [ ] 建立版本号体系
- [ ] 建立 package CLI entry point

## 任务 1.2：核心模型

实现：

```text
Protein
Species
DatasetManifest
RunManifest
ComponentSummary
HierarchyNode
SplitResult
```

要求：

- [ ] 内部 protein/species 使用 integer ID
- [ ] `original_id` 仅作为 metadata
- [ ] 所有 dataclass/model 都可序列化
- [ ] ID 不依赖 Python object 地址

## 任务 1.3：Deterministic ID system

- [ ] 输入文件 deterministic sort
- [ ] species ID 生成固定
- [ ] protein ID 生成固定
- [ ] 同一输入重复 prepare 结果完全一致
- [ ] 保存 input checksum

## 任务 1.4：Config system

实现 YAML + CLI override：

```text
search
similarity
edges
hierarchy
evolution
output
runtime
```

- [ ] 参数验证
- [ ] 默认值集中管理
- [ ] resolved config 写入 run directory

## 任务 1.5：Workspace

创建：

```text
run/
  run.yaml
  manifest.json
  run.db
  input/
  search/
  edges/
  components/
  hierarchy/
  evolution/
  results/
```

## 任务 1.6：Logging 与 Exceptions

- [ ] `OGProfilerError`
- [ ] `InputError`
- [ ] `SearchError`
- [ ] `EdgeConstructionError`
- [ ] `HierarchyError`
- [ ] `PhylogenyError`
- [ ] `CheckpointError`
- [ ] 文件日志 + console 日志
- [ ] 记录 stage / task / component context

## 任务 1.7：FASTA prepare

- [ ] FASTA parser
- [ ] duplicate ID 检查
- [ ] empty sequence 检查
- [ ] 非法字符处理策略
- [ ] 生成 `proteins.parquet`
- [ ] 生成 `species.parquet`
- [ ] 生成标准化 `proteins.faa`

## 测试

- [ ] 同一输入两次 prepare 的 ID 完全一致
- [ ] 输入顺序变化后，按既定排序规则结果一致
- [ ] duplicate protein ID 正确报错
- [ ] 空 proteome 正确报错

## 交付物

```text
src/ogprofiler/core/
src/ogprofiler/input/
src/ogprofiler/config.py
src/ogprofiler/logging.py
src/ogprofiler/exceptions.py
```

## 验收标准

```bash
ogprofiler prepare --proteomes testdata/ --out run/
```

可重复生成稳定 metadata。

---

# Phase 2：Hierarchy Engine 原型优先验证

> 这是最关键的技术验证阶段。建议在完整 Search/Edge Engine 重写前完成。

## 目标

直接读取 V1 `ssn.gml` 或简化 edge table，验证新的 component-centric hierarchical Leiden 是否显著改善性能。

## 任务 2.1：Legacy SSN importer

- [ ] 读取 V1 `ssn.gml`
- [ ] 提取 integer vertex IDs
- [ ] 提取指定 edge weight
- [ ] 转成简单 edge table

## 任务 2.2：Connected component extraction

- [ ] 从 edge table 得到 CC
- [ ] 记录 component size
- [ ] 可以按 component 导出 edge subset

## 任务 2.3：Minimal hierarchy model

先只实现：

```text
cluster_id
parent_id
component_id
depth
n_genes
resolution
quality
terminal_reason
```

## 任务 2.4：Leiden wrapper

实现：

```python
run_leiden(graph, gamma, method, weights, seed)
```

要求：

- [ ] 返回 membership，不返回 subgraphs
- [ ] 支持 RBER/RBConfiguration/CPM/Modularity
- [ ] 统一 random seed
- [ ] 记录 Leiden call count

## 任务 2.5：Resolution search prototype

实现第一版：

1. exponential search
2. local log-grid
3. lowest valid split selection

记录每次 candidate：

```text
gamma
k
quality
min_child_size
max_child_fraction
```

## 任务 2.6：DFS hierarchy engine

- [ ] 每个 component 单独加载
- [ ] stack/DFS
- [ ] parent-child 建立
- [ ] k-way split
- [ ] terminal reason
- [ ] 不在 hierarchy 结构中保存完整 descendant geneIDs

## 任务 2.7：Prototype benchmark instrumentation

记录：

- [ ] runtime
- [ ] peak RSS
- [ ] Leiden calls
- [ ] subgraph construction count
- [ ] hierarchy node count
- [ ] terminal family count

## 测试

- [ ] children disjoint
- [ ] union(children) = parent
- [ ] every protein belongs to exactly one terminal family
- [ ] no hierarchy cycles
- [ ] every non-root node has one parent

## 关键对照实验

同一 `ssn.gml`：

```text
V1 hierarchy engine
V2 prototype hierarchy engine
```

比较：

```text
runtime
peak RSS
Leiden calls
subgraph constructions
family count
family membership
hierarchy depth
```

## 验收标准

必须证明至少以下一项明显改善，且没有产生结构错误：

- runtime 显著下降；
- peak RSS 显著下降；
- subgraph constructions 大幅下降；
- Leiden calls 受到严格控制。

如果该阶段无法证明新引擎价值，应先调整 hierarchy strategy，再继续后续重写。

---

# Phase 3：Similarity Search Backend

## 目标

建立独立、可替换、可记录版本的 homolog search layer。

## 任务 3.1：SearchBackend protocol

定义：

```python
class SearchBackend:
    def build_database(...): ...
    def search(...): ...
    def version(...): ...
```

## 任务 3.2：DiamondBackend

- [x] database build
- [x] all-vs-all search
- [x] configurable sensitivity
- [x] configurable e-value
- [x] 明确 max-target-seqs 行为
- [x] 不依赖 shell redirection
- [x] `subprocess.run(..., check=True)`

## 任务 3.3：MMseqsBackend

- [ ] easy-search / search 模式设计
- [ ] temp directory management
- [ ] output schema 一致化

## 任务 3.4：BlastBackend

- [ ] 作为兼容/验证 backend
- [ ] 明确性能不作为主路径

## 任务 3.5：Search manifest

记录：

```text
backend
version
command
parameters
input checksum
output checksum
```

## 测试

- [ ] Dataset A 三种 backend 均能产生可解析 hits
- [x] failed subprocess 正确抛异常
- [x] resume 不重复执行已完成 search

## 验收标准

```bash
ogprofiler search --backend diamond
```

可以生成标准方向性 hit table。

---

# Phase 4：Hit Parser、NBS 与 Edge Engine

## 目标

替代 V1 大量 species-pair pickle/sparse matrix，形成统一 numeric edge pipeline。

## 任务 4.1：Hit table schema

字段：

```text
query_id
target_id
query_species
target_species
bitscore
identity
query_coverage
target_coverage
evalue
```

- [x] Arrow schema 固定
- [x] numeric dtype 明确
- [x] 分块读取支持

## 任务 4.2：Legacy NBS 实现

- [x] 完整复现 V1 top-bin 逻辑
- [x] 复现拟合公式
- [x] 明确小样本 fallback
- [x] 输出 `normalized_score`
- [x] 单元测试与 V1 数值对照

## 任务 4.3：Coverage filtering

配置：

```text
min_query_coverage
min_target_coverage
min_bidirectional_coverage
```

- [x] 默认值先设为 legacy-compatible 或明确记录改变
- [x] Dataset D 做 fusion sensitivity test

## 任务 4.4：Best hit / RBH

- [x] query group best hit
- [x] accepted tolerance
- [x] reciprocal join
- [x] same-species paralog handling
- [x] 无矩阵实现

## 任务 4.5：LRB threshold

- [x] 计算 per-query most-distant RBH threshold
- [x] 无 RBH fallback
- [x] 保留所有满足 threshold 的 homologs

## 任务 4.6：AR / ARB compatibility

- [x] 明确原始定义
- [x] 修正 V1 direction bug 后重新定义 expected behavior
- [x] 不追求错误行为 100% 兼容

## 任务 4.7：Edge canonicalization

将：

```text
A → B
B → A
```

join 为：

```text
u=min(A,B)
v=max(A,B)
score_uv
score_vu
```

## 任务 4.8：Symmetrization

实现：

```text
max
min
mean
geometric_mean
```

- [x] 默认策略暂不硬编码为“科学最优”
- [ ] benchmark 后选择默认
- [ ] edge weight 不依赖文件顺序

## 任务 4.9：retained_edges.parquet

最终字段建议：

```text
u
v
u_species
v_species
score_uv
score_vu
weight
coverage
edge_type
```

## 测试

- [x] V1/V2 legacy NBS 数值回归
- [x] RBH 边界 case
- [x] reverse-only hit
- [x] same-species paralog
- [x] symmetrization order invariance
- [x] Dataset D coverage sensitivity

## 验收标准

V2 可以不构造任何 BH/RBH/CM sparse matrix，生成与设计一致的 retained edge table。

---

# Phase 5：Connected Component Engine 与磁盘分区

## 目标

彻底移除“必须先构造全局 igraph SSN”的依赖。

## 任务 5.1：Union-Find

- [x] path compression
- [x] union by rank/size
- [x] integer protein IDs
- [x] 大规模 edge stream 支持

## 任务 5.2：Component index

输出：

```text
protein_id
component_id
```

## 任务 5.3：Component statistics

输出：

```text
component_id
n_vertices
n_edges
n_species
```

## 任务 5.4：Edge partitioning

- [x] 按 component 写 parquet partitions
- [x] 支持只读取单个 component
- [x] largest component 优先排序

## 任务 5.5：Singleton handling

无边 protein：

- [x] 作为 singleton terminal family
- [x] 不进入 Leiden
- [ ] 正确进入最终 members/families（Phase 6/7 汇总时验收）

## 测试

- [x] 与 igraph connected components 对照
- [x] random graph property test
- [x] singleton
- [x] disconnected edge chunks

## 验收标准

在没有 global igraph 的情况下，能够得到所有 component 并按需局部加载。

---

# Phase 6：正式 Hierarchical Leiden Core

## 目标

把 Phase 2 prototype 升级为正式生产内核。

## 任务 6.1：ComponentGraphLoader

- [x] 读取 component edge partition
- [x] remap global protein_id → local vertex index
- [x] 保存 local↔global mapping
- [x] 只在 worker 内构造 igraph

## 任务 6.2：Hierarchy store

实现 component-local writer：

```text
component_x.nodes.parquet
component_x.members.parquet
```

## 任务 6.3：Split acceptance framework

基础 metrics：

- [x] child_count
- [x] min_child_size
- [x] largest_child_fraction
- [x] tiny_fragment_fraction
- [x] quality
- [x] optional edge separation metrics

## 任务 6.4：Terminal reasons

实现：

```text
NO_SPLIT
ONE_SPECIES
MIN_SIZE
UNSTABLE
LOW_QUALITY
MAX_DEPTH
NO_EDGES
GAMMA_LIMIT
```

## 任务 6.5：Leiden stochasticity modes

### fast

- [x] 1 seed

### robust

- [x] 3 seeds
- [x] membership consistency metric

### publication

- [x] 5+ seeds
- [x] ARI/NMI
- [x] consensus or best-supported split

## 任务 6.6：Species bitset integration

每个 cluster 动态维护：

- [x] species bitmap
- [x] n_species

不保存 genomeIDs string list。

## 任务 6.7：DFS memory lifecycle

- [x] 完成 child 后释放不再需要对象
- [x] 避免 membership 的不必要 copy
- [x] 监测 worker RSS

## 测试

- [x] hierarchy invariants 全部通过
- [x] seed reproducibility
- [x] k-way split
- [x] one-species stop
- [x] max-depth stop
- [x] gamma-limit stop

## 验收标准

Dataset C/E 上 hierarchy engine 无异常内存增长，并比 V1 显著减少 subgraph/Leiden 重复调用。

---

# Phase 7：Scheduler、并行与 Resume

## 目标

实现 component-level parallel execution 和可靠断点续跑。

## 任务 7.1：Component scheduler

- [x] 按 component size descending 排序
- [x] worker 数由 config 控制
- [x] 主进程只分发 component_id

## 任务 7.2：Worker isolation

- [x] worker 自行加载 component
- [x] worker 自行写结果
- [x] 主进程只收轻量 summary
- [x] 禁止 Manager.Queue 传大对象

## 任务 7.3：run.db

SQLite 表至少包括：

```text
stage
task_id
status
started_at
completed_at
input_hash
output_path
error
algorithm_version
```

## 任务 7.4：Resume semantics

- [x] DONE + hash match → skip
- [x] RUNNING from crashed process → recover/reset
- [x] FAILED → configurable retry
- [x] algorithm/config change → INVALID

## 任务 7.5：Fault isolation

一个 component 失败时：

- [x] 其他 component 继续
- [x] 最终报告 failed task
- [x] 可以只重跑失败 component

## 测试

- [x] kill process 后 resume
- [x] worker exception
- [x] 部分 output 损坏
- [x] config 改变后 invalidation

## 验收标准

大型 run 可以安全中断和继续，不重复完成的 component。

---

# Phase 8：Network Evolution Annotation

## 目标

在不宣称 network hierarchy 等价于 gene tree 的前提下，对 hierarchy 节点提供可解释的进化倾向标注。

## 任务 8.1：Event schema

字段：

```text
network_event
overlap_score
pairwise_overlap_summary
confidence
```

事件：

```text
SPECIATION_LIKE
DUPLICATION_LIKE
MIXED
POLYTOMY
AMBIGUOUS
SPECIES_SPECIFIC
```

## 任务 8.2：Binary overlap

- [x] `S1 & S2`
- [x] overlap count
- [x] normalized overlap
- [x] threshold config

## 任务 8.3：Polytomy overlap

- [x] child pair overlap matrix
- [x] all-zero overlap detection
- [x] mixed overlap detection
- [x] duplication-like pattern detection

## 任务 8.4：Legacy event mapping

为了兼容旧结果，可增加 export mapping：

```text
I
II
III-1
III-2
III-3
```

但内部逻辑不再依赖这些标签。

## 测试

- [x] disjoint child species
- [x] identical child species
- [x] partial overlap
- [x] k-way polytomy
- [x] one species family

## 验收标准

每个 internal hierarchy node 都可得到明确 `network_event` 或 `AMBIGUOUS`，并保留连续 overlap metrics。

---

# Phase 9：Terminal Families 与结果导出

## 目标

形成 OGProfiler 2.0 MVP 的稳定最终输出。

## 任务 9.1：Family ID assignment

- [x] deterministic OG IDs
- [x] 与 component traversal 顺序解耦
- [x] 可重复运行一致

## 任务 9.2：families.tsv

建议字段：

```text
family_id
component_id
cluster_id
n_genes
n_species
terminal_reason
network_event
```

## 任务 9.3：members.tsv

```text
family_id
protein_id
species_id
original_id
```

## 任务 9.4：hierarchy.tsv

```text
cluster_id
parent_id
component_id
depth
n_genes
n_species
resolution
quality
child_count
terminal_reason
```

## 任务 9.5：events.tsv

```text
cluster_id
network_event
overlap_score
confidence
```

## 任务 9.6：FASTA export

- [x] 按 family 导出 fasta
- [x] 支持只导出选定 family
- [x] 避免默认生成百万个小文件时失控

## 任务 9.7：Graph export

按需：

```bash
ogprofiler export graph --component 123 --format graphml
```

而不是默认输出整个 global GML。

## 验收标准

OGProfiler 2.0 MVP 可以从 FASTA 输入稳定生成：

```text
families.tsv
members.tsv
hierarchy.tsv
events.tsv
```

---

# Phase 10：Orthology Engine

## 目标

利用已经建立的 rooted hierarchy 推导 ortholog candidate，完全移除 V1 中 leaf-combination + shortest-path 模式。

## 任务 10.1：Hierarchy traversal

对于 speciation-like node：

- [x] 获取 child descendant terminal families
- [x] 获取 child gene sets
- [x] cross-child pairing

## 任务 10.2：Same-species filtering

- [x] 不输出同 species pair
- [x] 明确 in-paralog / co-ortholog 策略

## 任务 10.3：Streaming writer

- [x] generator 输出
- [x] `.tsv.zst`
- [x] 不缓存全部 pair
- [x] progress 以 pair chunk 为单位

## 任务 10.4：可选输出

默认允许关闭：

```text
--emit-pairwise-orthologs false
```

避免用户只需要 family/hierarchy 时生成巨大二次方输出。

## 测试

- [x] 小型已知 hierarchy 手工验证 pairs
- [x] 与 V1 在简单数据上的 pair overlap
- [x] 大 family 内存测试

## 验收标准

Pairwise ortholog generation 的内存与输出规模解耦，峰值内存不随 pair 数线性增长。

---

# Phase 11：Phylogenetic Refinement

## 目标

为 network hierarchy 提供可选的第二证据层。

## 任务 11.1：Family selection policy

默认优先：

- [x] AMBIGUOUS
- [x] MIXED
- [x] large family
- [x] duplication-rich
- [x] user-selected families

## 任务 11.2：AlignmentBackend

- [x] MAFFT adapter
- [x] external binary discovery
- [x] version capture
- [x] subprocess standardization

## 任务 11.3：TreeBackend

至少实现一种：

- [x] FastTree

后续可扩展：

- [ ] IQ-TREE
- [ ] RAxML-NG

## 任务 11.4：Rooting strategy

支持：

```text
midpoint
species-tree-aware
user-provided/outgroup
```

Midpoint 不再是唯一默认科学结论来源。

## 任务 11.5：Species-tree input

- [x] Newick parser
- [x] species ID mapping
- [x] missing species validation
- [x] tree pruning policy

## 任务 11.6：ReconciliationBackend

接口：

```python
class ReconciliationBackend:
    def annotate(gene_tree, species_tree, mapping): ...
```

输出：

```text
phylo_event
confidence
supporting_node
```

## 任务 11.7：Network/Phylo event 并列

禁止覆盖：

```text
network_event
phylo_event
```

如果冲突，应保留冲突状态。

## 验收标准

用户可对指定 family 运行：

```bash
ogprofiler annotate --phylogenetic-refinement
```

得到独立于 network_event 的 phylo_event。

---

# Phase 12：科学 Benchmark 与算法选择

## 目标

从“软件能跑”进入“方法可信”。

## 任务 12.1：Parameter benchmark matrix

测试：

- [x] different normalization
- [x] different coverage thresholds
- [x] symmetrization max/min/mean/geometric_mean
- [x] different Leiden methods
- [x] different gamma search strategy
- [x] different split acceptance thresholds
- [x] different random seeds

## 任务 12.2：与其他方法比较

至少规划：

```text
OGProfiler 1
OGProfiler 2
OrthoFinder
其他适用方法
```

## 任务 12.3：Family metrics

- [x] family number
- [x] family size distribution
- [x] species coverage
- [x] single-copy family recovery
- [x] large-family fragmentation

## 任务 12.4：Evolution metrics

- [x] species-overlap consistency
- [x] gene-tree concordance
- [x] reconciliation consistency

## 任务 12.5：Orthology metrics

有 reference/simulation 时：

- [x] precision
- [x] recall
- [x] F1

## 验收标准

为 OGProfiler 2 默认参数提供数据驱动的理由，而不是凭经验选择。

---

# Phase 13：Synthetic Evolution Benchmark

## 目标

建立拥有 ground truth genealogy 的测试体系。

## 任务 13.1：模拟场景设计

系统改变：

- [ ] sequence divergence
- [ ] duplication rate
- [ ] gene loss
- [ ] family expansion
- [ ] horizontal transfer（如果未来纳入）
- [ ] fusion/domain architecture variation

## 任务 13.2：Ground truth

每个模拟数据记录：

```text
true family
true genealogy
true duplication nodes
true speciation nodes
true ortholog pairs
```

## 任务 13.3：Recovery metrics

- [ ] terminal family recovery
- [ ] hierarchy similarity
- [ ] event classification accuracy
- [ ] ortholog precision/recall

## 任务 13.4：Leiden applicability map

最终希望回答：

> 在什么 divergence、duplication rate、loss rate 和 fusion rate 下，hierarchical Leiden 仍能可靠恢复进化结构？

## 验收标准

形成可用于论文 Methods/Results 的系统 benchmark 数据。

---

# Phase 14：Extreme-scale Optimization

> 只有在 2.0 科学逻辑稳定后才开始。

## 目标

支持超大型 protein datasets 与异常巨大 CC。

## 任务 14.1：Subtree scheduling

- [ ] 大 component 首次 split 后释放 child tasks
- [ ] scheduler 可以递归提交大 subtree
- [ ] parent-child ID 仍保持 deterministic

## 任务 14.2：Memory-mapped data

评估：

- [ ] Arrow memory map
- [ ] NumPy memmap
- [ ] local index arrays

## 任务 14.3：Disk I/O profiling

- [ ] parquet row-group tuning
- [ ] compression benchmark
- [ ] sequential vs random access

## 任务 14.4：HPC/SLURM

- [ ] component task manifest
- [ ] array jobs
- [ ] results merge
- [ ] failure resume

## 任务 14.5：Distributed option

仅在实际规模证明必要时考虑：

- [ ] Ray/Dask/自定义 scheduler 评估
- [ ] 不提前引入分布式复杂度

## 验收标准

超大型数据的主要限制从 Python object overhead 转变为真正的算法/数据规模本身。

---

# Phase 15：用户体验、文档与发布

## 目标

把研究代码变成可维护工具。

## 任务 15.1：CLI polish

- [ ] `ogprofiler run`
- [ ] `ogprofiler status`
- [ ] `ogprofiler inspect`
- [x] `ogprofiler export`

## 任务 15.2：README

至少包括：

- [ ] project concept
- [ ] installation
- [ ] quick start
- [ ] output explanation
- [ ] scientific caveats

## 任务 15.3：Method documentation

必须明确：

- [x] Network hierarchy ≠ gene tree
- [x] `network_event` 是 heuristic
- [x] `phylo_event` 来自系统发育分析
- [x] terminal family 的正式定义

## 任务 15.4：Reproducibility guide

- [ ] seed
- [ ] version pinning
- [ ] manifest
- [ ] run.yaml
- [ ] external tool versions

## 任务 15.5：Release checklist

- [ ] unit tests pass
- [ ] integration tests pass
- [ ] regression tests pass
- [ ] benchmark completed
- [ ] changelog updated
- [ ] version tagged

---

# OGProfiler 2.0 MVP 最小任务集合

如果希望尽快做出第一个真正可用版本，可把 2.0 边界压缩为以下任务：

## 必须完成

- [ ] Phase 0：V1 freeze + baseline
- [ ] Phase 1：core/input/config/workspace
- [ ] Phase 2：hierarchy prototype
- [ ] Phase 3：至少 DiamondBackend
- [ ] Phase 4：legacy NBS + LRB + edge canonicalization
- [ ] Phase 5：Union-Find components
- [ ] Phase 6：正式 hierarchy engine
- [ ] Phase 7：component scheduler + resume
- [x] Phase 8：network events
- [x] Phase 9：families/members/hierarchy/events output

## 可以推迟到 2.1

- [x] pairwise ortholog generation
- [x] full phylogenetic refinement
- [x] species-tree reconciliation
- [ ] subtree-level parallelism
- [ ] distributed/HPC extension

---

# 建议的开发顺序（实际执行版）

如果从下一次 commit 开始，我建议严格按以下顺序推进：

```text
1. Freeze V1
2. Build regression datasets
3. New project skeleton
4. Legacy SSN importer
5. New hierarchy prototype
6. Benchmark V1 vs V2 hierarchy
7. Finalize hierarchy architecture
8. Build search backend
9. Build edge engine
10. Build Union-Find components
11. Connect full FASTA → hierarchy pipeline
12. Add network events
13. Stabilize 2.0 outputs
14. Scientific benchmarks
15. Orthology engine
16. Phylogenetic refinement
17. Extreme-scale optimization
```

这里最关键的是：

> **第 5–6 步必须早于完整 homology-search 重写。**

因为 OGProfiler 当前最重要的技术风险不是“DIAMOND wrapper 是否漂亮”，而是新的 hierarchy engine 能否真正解决性能瓶颈。

---

# 每个 Pull Request 的基本验收规则

所有影响核心算法的 PR 至少回答：

1. 这个改动改变了什么行为？
2. 是否影响 family membership？
3. 是否影响 hierarchy topology？
4. 是否影响 reproducibility？
5. 是否影响 runtime / peak RSS？
6. 哪个 regression dataset 覆盖了该行为？
7. 结果变化属于 bug fix、performance-only 还是 scientific-model change？

如果无法回答第 6 和第 7 项，则不建议合并核心算法改动。

---

# Definition of Done：OGProfiler 2.0

OGProfiler 2.0 可以定义为完成，当满足：

## 工程层面

- [ ] 不再依赖单文件主程序
- [ ] 不再通过 `eval()` 调度 pipeline
- [ ] 不再把全局 igraph 传给 multiprocessing workers
- [ ] hierarchy worker 直接写 component-local 结果
- [ ] 支持 checkpoint/resume
- [ ] 同一 input/config/seed 可重复

## 性能层面

- [ ] hierarchy 峰值内存主要取决于当前 component
- [ ] gamma search 不调用 `partition.subgraphs()`
- [ ] Leiden call count 有明确上界和日志
- [ ] large output 支持 streaming

## 科学层面

- [ ] network hierarchy 与 phylogeny 概念分离
- [ ] terminal family 有正式定义
- [ ] network_event 明确为 heuristic
- [ ] k-way split 被原生支持
- [ ] edge weight 与 input ordering 无关

## 输出层面

至少稳定生成：

```text
families.tsv
members.tsv
hierarchy.tsv
events.tsv
run.yaml
manifest.json
```

## 测试层面

- [ ] unit tests
- [ ] integration tests
- [ ] V1 regression tests
- [ ] hierarchy invariant tests
- [ ] performance benchmark

全部通过。

---

# 最终开发原则

OGProfiler 2 的开发过程中始终坚持三条原则：

### 原则 1：先证明算法执行模型，再扩大功能范围

先证明新的 hierarchy engine 可扩展，再增加更多系统发育和导出功能。

### 原则 2：结构优化与科学模型修改分开进行

先做 legacy-compatible implementation，再逐项引入 coverage、symmetrization、new resolution criterion 等科学变化。

### 原则 3：任何大对象都必须有明确生命周期

对于 graph、membership、pair list、sequence、edge table，必须明确：

```text
什么时候创建
在哪里保存
谁拥有
什么时候释放
是否需要跨进程传输
```

如果一个对象没有明确生命周期，就不应该进入核心 pipeline。

这将是避免 OGProfiler 2 再次演化成 V1 式高耦合、高内存程序的最重要工程纪律。
