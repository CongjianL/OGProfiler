# OGProfiler 2 设计方案

> 文档定位：OGProfiler 2 的总体设计说明书，定义项目目标、科学边界、架构原则、数据模型、核心算法、并行与存储模型、输出体系、测试原则与版本边界。
>
> 基础版本：以当前 `OGProfiler.py` 作为 legacy/reference implementation，不在原脚本上持续增量修改，而是重新实现 OGProfiler 2 内核。

---

## 1. 项目目标

OGProfiler 2 的目标不是简单重构原有脚本，而是重新设计一个能够支持大规模蛋白序列数据的层级蛋白家族推断框架。

建议将 OGProfiler 2 定义为：

> **基于序列相似性网络、多尺度社区检测和可选系统发育验证的层级蛋白家族推断框架。**

英文可表述为：

> **Hierarchical protein family inference from sequence similarity networks with phylogenetic refinement.**

OGProfiler 2 应同时满足两个目标：

1. 在计算层面，能够扩展到远大于当前版本的数据规模，并显著降低图复制、进程间通信和层级成员重复存储造成的内存与时间开销。
2. 在科学层面，将“网络层级结构”和“真实系统发育关系”明确分层，避免把 Leiden 社区层级直接等同于基因树。

---

## 2. 核心科学模型

OGProfiler 2 将整个推断过程拆分为两个证据层。

### 2.1 Network hierarchy

负责回答：

> 哪些蛋白在 sequence space 中形成稳定的 family/subfamily 结构？

证据来源：

- homologous sequence search
- normalized similarity score
- sequence similarity network
- Leiden community detection
- hierarchical partition

输出是：

> **Network Family Hierarchy**

它表示蛋白序列相似性空间中的多尺度家族结构，而不是严格意义上的系统发育树。

### 2.2 Evolutionary interpretation

负责回答：

> 哪些层级节点更可能对应 speciation、duplication 或复杂/不确定的进化事件？

证据来源可以分为两级：

1. 基于 child species composition 的 network-level heuristic；
2. 对重要或模糊家族进行 gene tree + species tree reconciliation。

因此核心关系是：

```text
Sequence similarity
        │
        ▼
Homology graph
        │
        ▼
Hierarchical Leiden
        │
        ▼
Network Family Hierarchy
        │
        ├───────────────┐
        ▼               ▼
Terminal families   ambiguous/internal nodes
                        │
                        ▼
                   Gene phylogeny
                        │
                        ▼
             Evolutionary annotation
```

---

## 3. OGProfiler 2 的六项核心设计目标

### 3.1 峰值内存由“全局数据规模”降低到“当前最大 connected component”

旧执行模型大致为：

```text
全部序列
  ↓
全部 pairwise matrices
  ↓
全局 SSN
  ↓
不断复制 subgraph
  ↓
层级网络
```

OGProfiler 2 改为：

```text
全部序列
  ↓
edge table
  ↓
connected component index
  ↓
CC1 → hierarchy
CC2 → hierarchy
CC3 → hierarchy
...
```

Hierarchy 阶段的目标峰值内存应近似取决于：

\[
O(V_{largestCC}+E_{largestCC})
\]

而不是所有网络及其重复子图之和。

### 3.2 Hierarchy internal node 不重复保存完整 gene list

内部节点仅保存结构和统计信息，例如：

```text
cluster_id
parent_id
component_id
depth
n_genes
n_species
resolution
quality
split_status
network_event
phylo_event
```

蛋白 membership 只在 terminal family 层保存一次：

```text
protein_id → terminal_cluster_id
```

### 3.3 Leiden 不再被强制产生二叉树

Leiden 的职责是识别 community，而不是制造 phylogenetic bifurcation。

如果稳定结果是 3、4 或更多 communities，则保留真正的 k-way split：

```text
parent
 ├─ child A
 ├─ child B
 ├─ child C
 └─ child D
```

二叉拓扑需求应交由 gene tree 模块处理。

### 3.4 并行单位改为 connected component

每个 worker 获取一个 `component_id`，自行加载该 component 的 edge data，并从 root 一直处理到 terminal families。

原则：

- 不在每个 hierarchy level 重建进程池；
- 不通过 multiprocessing 参数传递完整 igraph；
- 不用 Manager.Queue 回传大规模图或 membership 列表；
- worker 直接写 component-level 结果文件。

### 3.5 Global SSN 不再是必须长期存在的数据对象

完整 SSN 只是一种概念表示，而不是强制数据结构。

OGProfiler 2 的核心流程应为：

```text
homology hits
     ↓
filtered/symmetrized edge stream
     ↓
Union-Find
     ↓
component ID
     ↓
partitioned edge storage
```

Leiden worker 只在需要时加载当前 component 并构造 igraph。

### 3.6 全流程可复现

每个 run 必须记录：

- OGProfiler version
- algorithm version
- input checksums
- search backend version
- command line
- resolved configuration
- random seed
- Leiden method
- resolution strategy
- normalization method
- edge symmetrization method
- coverage threshold
- event inference mode

同样的输入、参数、seed 和软件版本应得到可复现结果。

---

## 4. 推荐的软件架构

```text
OGProfiler2/
│
├── pyproject.toml
├── README.md
├── LICENSE
├── CHANGELOG.md
│
├── src/
│   └── ogprofiler/
│       ├── cli.py
│       ├── config.py
│       ├── logging.py
│       ├── exceptions.py
│       │
│       ├── core/
│       │   ├── models.py
│       │   ├── ids.py
│       │   ├── manifest.py
│       │   └── workspace.py
│       │
│       ├── input/
│       │   ├── fasta.py
│       │   └── proteomes.py
│       │
│       ├── search/
│       │   ├── base.py
│       │   ├── diamond.py
│       │   ├── mmseqs.py
│       │   └── blast.py
│       │
│       ├── similarity/
│       │   ├── parser.py
│       │   ├── normalization.py
│       │   ├── best_hits.py
│       │   ├── edge_filter.py
│       │   └── symmetrize.py
│       │
│       ├── graph/
│       │   ├── components.py
│       │   ├── union_find.py
│       │   └── loader.py
│       │
│       ├── hierarchy/
│       │   ├── leiden.py
│       │   ├── resolution.py
│       │   ├── engine.py
│       │   ├── scheduler.py
│       │   └── validation.py
│       │
│       ├── evolution/
│       │   ├── events.py
│       │   ├── overlap.py
│       │   ├── alignment.py
│       │   ├── genetree.py
│       │   └── reconciliation.py
│       │
│       ├── orthology/
│       │   ├── families.py
│       │   └── pairs.py
│       │
│       ├── storage/
│       │   ├── edges.py
│       │   ├── hierarchy.py
│       │   ├── checkpoint.py
│       │   └── metadata.py
│       │
│       └── export/
│           ├── tables.py
│           ├── fasta.py
│           ├── graph.py
│           └── legacy.py
│
├── tests/
│   ├── unit/
│   ├── integration/
│   ├── regression/
│   └── performance/
│
└── legacy/
    └── OGProfiler_v1.py
```

### 4.1 模块依赖原则

依赖方向固定为：

```text
Input
  ↓
Search
  ↓
Similarity edges
  ↓
Connected components
  ↓
Leiden hierarchy
  ↓
Evolution annotation
  ↓
Orthology / Export
```

禁止跨层直接访问内部实现。

---

## 5. 核心数据模型

### 5.1 Protein

内部主键使用整数，不再把 `G123|g4567` 作为核心主键。

```python
Protein:
    protein_id: int
    species_id: int
    original_id: str
    length: int
```

示例：

```text
protein_id   species_id   original_id   length
0            0            geneA         287
1            0            geneB         412
2            1            XP_001        295
```

整数 ID 有利于：

- igraph
- NumPy
- Arrow/Parquet
- join/sort
- memory efficiency

### 5.2 Species

```text
species_id
species_name
source_file
```

输入文件必须 deterministic sort 后再生成 species ID。

### 5.3 Sequence storage

序列 metadata 与序列本身分离：

```text
proteins.parquet
species.parquet
proteins.faa
```

避免把所有 sequence 保存到 Python dictionary/pickle。

---

## 6. Homology Search 设计

### 6.1 Backend abstraction

定义统一接口：

```python
class SearchBackend:
    def build_database(...): ...
    def search(...): ...
    def version(...): ...
```

实现：

- `DiamondBackend`
- `MMseqsBackend`
- `BlastBackend`

### 6.2 默认搜索模型

优先评估：

```text
all_proteins.faa
       ↓
global DIAMOND/MMseqs database
       ↓
all-vs-all search
```

而不是强制 species × species 建库和搜索。

如果 legacy normalization 确实要求 species-pair context，可以在统一 hit table 上按 species pair 分组处理，而不是创建大量独立文件和矩阵。

### 6.3 搜索输出规范

方向性 hit table：

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
normalized_score
```

---

## 7. Similarity / Edge Engine

### 7.1 内部存储从 matrix/pickle 改为 columnar table

取消大规模：

```text
similarity*.pic
bitScores*.pic
NB*.pic
BH*.pic
RBHs*.pic
CM*.pic
```

推荐使用 Arrow/Parquet 保存数值型边数据。

### 7.2 NBS 初期保留 legacy 模型

OGProfiler 2 第一阶段不要同时改变工程架构和生物学评分模型。

提供：

```text
normalization = legacy_nbs
```

先尽量复现 V1 的 NBS：

\[
BS_{norm} = \frac{BS}{10^b L^a}
\]

后续版本再增加：

```text
normalization:
    legacy_nbs
    none
    new_model
```

### 7.3 Coverage 成为正式过滤条件

建议支持：

```text
--min-query-cover
--min-target-cover
--min-bidirectional-cover
```

并正式定义：

\[
coverage = \min(C_q, C_t)
\]

减少 domain-only hit、fusion/chimera 对 SSN 的错误桥接。

### 7.4 RBH / LRB 无矩阵化实现

核心思想可以保留：

\[
T_i = \min_{j \in RBH(i)} NBS_{ij}
\]

保留：

\[
NBS_{ik} \ge T_i
\]

但实现方式改为：

```text
hits
 ↓
groupby(query)
 ↓
best hit
 ↓
RBH join
 ↓
threshold
 ↓
filter
```

不再生成 BH/RBH/CM sparse matrices。

### 7.5 Undirected edge 对称化

最终 canonical edge 必须满足：

```text
u < v
w(u,v) = w(v,u)
```

定义：

```python
symmetrize(forward, reverse, method)
```

候选方法：

- max
- min
- mean
- geometric_mean

默认方法通过 benchmark 决定，而不是依赖 genome 输入顺序。

---

## 8. Connected Component Engine

采用 Union-Find / Disjoint Set Union。

扫描每条 retained edge：

```python
union(u, v)
```

输出：

```text
components.parquet

protein_id
component_id
```

边文件根据 `component_id` 分区：

```text
components/edges/
    component=000000/
    component=000001/
    component=000002/
```

这样 worker 可以按 component 局部加载。

---

## 9. Hierarchy 数据结构

Hierarchy 不再用 igraph 长期存储。

推荐结构：

```python
HierarchyNode:
    cluster_id
    parent_id
    component_id
    depth
    n_genes
    n_species
    resolution
    quality
    child_count
    split_status
    terminal_reason
    network_event
    phylo_event
```

Hierarchy forest 保存为：

```text
hierarchy_nodes.parquet
```

Terminal membership 保存为：

```text
protein_terminal_family.parquet
```

### 9.1 避免 ancestor 重复 membership

建议 DFS 后为 terminal family 定义 leaf order，internal node 记录：

```text
leaf_start
leaf_end
```

需要恢复 ancestor membership 时，通过 terminal family 区间再 join protein membership。

---

## 10. Hierarchical Leiden Engine

这是 OGProfiler 2 的核心。

### 10.1 Split API

```python
result = split_cluster(
    graph,
    cluster_context,
    config,
)
```

返回：

```python
SplitResult(
    accepted=True,
    membership=...,
    resolution=...,
    quality=...,
    stability=...,
    terminal_reason=None,
)
```

不返回 `partition.subgraphs()`。

### 10.2 每个 cluster graph 只物化一次

Gamma 搜索期间：

```python
g = build_induced_graph(vertex_ids)

for gamma in candidates:
    partition = run_leiden(g, gamma)
    membership = partition.membership
    k = len(partition)
```

禁止在 resolution search 中反复创建 subgraph。

### 10.3 不再搜索“恰好两个 community”

新的目标是：

> 找到最早出现的、稳定的、符合最小社区约束的有效 split。

允许：

```text
k = 2, 3, 4, ...
```

### 10.4 Resolution strategy

建议第一版使用三步策略。

#### Step A：Exponential search

从 `gamma0` 开始：

\[
\gamma, 2\gamma, 4\gamma, 8\gamma, ...
\]

直到第一次出现 `k > 1` 或达到 `gamma_max`。

#### Step B：局部 log-grid refinement

在 `[gamma_low, gamma_high]` 测试有限数量候选值。

记录：

```text
gamma
k
quality
min_child_size
max_child_fraction
stability
```

#### Step C：选择最低有效 resolution

优先接受最粗粒度、但已出现稳定合理 subdivision 的 resolution。

---

## 11. Split Acceptance Rules

Resolution search 与 split acceptance 必须分开。

推荐基础规则：

1. `child_count >= 2`
2. 每个 child ≥ `min_family_size`
3. `largest_child / parent_size < max_dominance`
4. tiny fragments 比例不超过阈值
5. robust 模式下 stability ≥ threshold
6. split 后确实产生新的结构信息

可附加统计：

- inter-community edge fraction
- intra-community edge retention
- modularity / partition quality
- size entropy

---

## 12. Terminal Reason

终止切割必须显式记录原因：

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

这比 V1 中简单的 `+` 更具有可解释性和可调试性。

---

## 13. Leiden stochasticity

推荐三种运行模式。

### fast

- 每个 gamma 运行 1 个 seed
- 用于大规模探索

### robust

- 候选 gamma 使用 3 个 seed
- 比较 partition consistency

### publication

- 关键 family 使用 5 个或更多 seed
- 计算 ARI/NMI/consensus

如果不同 seed 结果不稳定：

```text
terminal_reason = UNSTABLE
```

这本身也是对 sequence-space 边界模糊性的描述。

---

## 14. Hierarchy Engine 执行模型

推荐 DFS：

```python
def process_component(component_id):
    graph = load_component_graph(component_id)

    root = create_root_cluster()
    stack = [root]

    while stack:
        cluster = stack.pop()

        subgraph = build_subgraph(
            graph,
            cluster.vertices,
        )

        result = split_cluster(subgraph, config)

        if not result.accepted:
            mark_terminal(cluster, result.terminal_reason)
            continue

        children = create_children(
            cluster,
            result.membership,
        )

        save_split(cluster, children)
        stack.extend(children)
```

优先 DFS 的原因：

- 更早释放中间对象；
- 峰值内存更容易控制；
- 适合 component-local processing。

---

## 15. Scheduler 与并行模型

### 15.1 第一阶段：component-level parallelism

按照 component size 从大到小调度：

```text
Worker 1 → CC1
Worker 2 → CC2
Worker 3 → CC3
Worker 4 → CC4
```

减少长尾等待。

### 15.2 Worker 输出原则

worker 直接写：

```text
component_000123.hierarchy.parquet
component_000123.members.parquet
```

只向主进程返回轻量 summary：

```python
ComponentSummary(
    component_id=123,
    n_clusters=537,
    n_terminal=269,
    status="complete",
)
```

### 15.3 第二阶段可选优化：subtree scheduling

对于超大型 component，在第一层/第二层分裂后，将大 child subtree 重新提交 scheduler。

此功能属于后续性能优化，不作为 2.0 MVP 前置条件。

---

## 16. Species Set 数据结构

内部不使用大量 `set(str)` 或空格拼接的 genomeIDs。

优先使用：

- integer bitset
- Roaring bitmap（如未来有需要）

如果 species 数量适中，Python integer bitset 已非常有效：

```python
overlap = child1_species & child2_species
n_overlap = overlap.bit_count()
```

---

## 17. Network Evolution Annotation

V1 的 I / II / III-* 事件保留 legacy mapping，但 V2 内部使用更加明确的术语：

```text
SPECIATION_LIKE
DUPLICATION_LIKE
MIXED
POLYTOMY
AMBIGUOUS
SPECIES_SPECIFIC
```

字段必须命名为：

```text
network_event
```

而不是直接 `event`。

### 17.1 Binary node

令两个 child species sets 为 `S1`, `S2`。

如果：

\[
S_1 \cap S_2 = \emptyset
\]

可记为：

```text
SPECIATION_LIKE
```

并记录连续 overlap 指标：

\[
Overlap = \frac{|S_1 \cap S_2|}{\min(|S_1|, |S_2|)}
\]

如果 overlap 较高则偏向 `DUPLICATION_LIKE`，中间状态为 `MIXED`。

### 17.2 Polytomy

对于多 child node，计算 pairwise species-overlap matrix：

\[
O_{ij} = |S_i \cap S_j|
\]

根据整体 overlap pattern 判断：

- `SPECIATION_LIKE_POLYTOMY`
- `DUPLICATION_LIKE`
- `AMBIGUOUS`

不为了适配二叉事件模型而先强制 Leiden 二分。

---

## 18. Phylogenetic Refinement

Phylogenetic refinement 是第二证据层，而不是所有 family 的强制步骤。

优先处理：

- `AMBIGUOUS`
- `MIXED`
- 大型 family
- duplication-rich family
- 用户指定 family

Pipeline：

```text
family sequences
       ↓
MSA
       ↓
gene tree
       ↓
rooting
       ↓
species mapping
       ↓
event inference / reconciliation
```

模块接口：

```python
class AlignmentBackend: ...
class TreeBackend: ...
class ReconciliationBackend: ...
```

外部程序路径不得硬编码。

---

## 19. Species Tree 支持

CLI 可接受：

```bash
ogprofiler run \
    --proteomes genomes/ \
    --species-tree species.nwk
```

无 species tree：

```text
network hierarchy
+
species-overlap heuristic
```

有 species tree：

```text
network hierarchy
+
gene tree
+
species-tree reconciliation
```

输出并列保存：

```text
network_event
phylo_event
event_confidence
```

两者不能互相覆盖。

---

## 20. Terminal Family 的正式定义

OGProfiler 2 中的 terminal family 应定义为：

> 在当前 edge definition、resolution strategy、split acceptance 和 stability criterion 下，没有发现可信进一步 network subdivision 的蛋白集合。

推荐术语：

> **terminal network family**

避免直接称为 evolutionary family。

---

## 21. Orthology Engine

V1 的 leaf combinations + shortest path 模型应废弃。

如果 hierarchy 中节点 S 为 `SPECIATION_LIKE`：

```text
        S
       / \
      A   B
```

则 ortholog candidate 可直接来自：

\[
Genes(A) \times Genes(B)
\]

并过滤 same-species pairs。

### 21.1 Streaming 输出

禁止构建巨大的：

```python
ortholog_pairs = []
```

改为 generator + writer：

```python
for pair in generate_pairs(...):
    writer.write(pair)
```

输出：

```text
ortholog_pairs.tsv.zst
```

Pairwise ortholog 输出应为可选功能，因为其理论输出规模本身可能接近 `O(n²)`。

---

## 22. Workspace 与存储

推荐目录：

```text
run/
│
├── run.yaml
├── manifest.json
├── run.db
│
├── input/
│   ├── proteins.parquet
│   ├── species.parquet
│   └── proteins.faa
│
├── search/
│   └── hits/
│
├── edges/
│   └── retained_edges.parquet
│
├── components/
│   ├── index.parquet
│   └── edges/
│
├── hierarchy/
│   ├── components/
│   ├── nodes.parquet
│   └── members.parquet
│
├── evolution/
│   ├── alignments/
│   ├── trees/
│   └── events.parquet
│
└── results/
    ├── families.tsv
    ├── hierarchy.tsv
    ├── members.tsv
    ├── events.tsv
    └── ortholog_pairs.tsv.zst
```

### 22.1 数据格式职责

- Parquet：内部数据存储
- TSV：跨软件交换与最终表格输出
- GML/GraphML：可视化导出
- FASTA：序列输出
- Newick：系统发育树

GML 不再作为核心中间存储格式。

---

## 23. Checkpoint / Resume

推荐在 workspace 中维护 SQLite：

```text
run.db
```

记录：

```text
stage
task_id
status
started
completed
input_hash
output_path
algorithm_version
```

状态：

```text
PENDING
RUNNING
DONE
FAILED
INVALID
```

支持：

```bash
ogprofiler run --resume
```

Resume 逻辑基于任务状态和 input hash，而不是简单判断文件是否存在。

---

## 24. CLI 设计

完整流程：

```bash
ogprofiler run
```

单阶段命令：

```bash
ogprofiler prepare
ogprofiler search
ogprofiler build-edges
ogprofiler components
ogprofiler cluster
ogprofiler annotate
ogprofiler orthologs
ogprofiler export
```

状态：

```bash
ogprofiler status
```

检查对象：

```bash
ogprofiler inspect family OG000123
ogprofiler inspect component 517
```

---

## 25. 配置体系

示例：

```yaml
search:
  backend: diamond
  evalue: 1e-3
  threads: 32

similarity:
  normalization: legacy_nbs

edges:
  method: lrb
  min_query_coverage: 50
  min_target_coverage: 50
  symmetrization: max

hierarchy:
  method: rber
  seed: 42
  min_family_size: 2
  resolution_strategy: adaptive
  stability_mode: robust

evolution:
  network_overlap_threshold: 0.0
  phylogenetic_refinement: false
```

命令行参数覆盖 YAML。

---

## 26. Pipeline orchestration

禁止使用 `eval()` 和动态 methodcaller 驱动主流程。

推荐显式 orchestration：

```python
def run_pipeline(config):
    dataset = prepare(config)
    hits = run_search(dataset, config)
    edges = build_edges(hits, config)
    components = build_components(edges)
    hierarchy = infer_hierarchy(components, config)
    events = annotate_hierarchy(hierarchy, config)
    export_results(...)
```

保持简单、显式、可测试。

---

## 27. 外部命令调用规范

禁止：

```python
subprocess.Popen("diamond ... > output", shell=True)
```

使用：

```python
subprocess.run(
    ["diamond", "blastp", "--query", query, ...],
    check=True,
)
```

并记录：

- command
- stdout/stderr
- return code
- tool version

---

## 28. 异常体系

定义：

```text
OGProfilerError
InputError
SearchError
EdgeConstructionError
HierarchyError
PhylogenyError
CheckpointError
```

CLI 最外层统一处理。

库代码不得随意 `sys.exit(0)`。

---

## 29. 测试体系

### 29.1 Unit tests

至少覆盖：

```text
test_union_find
test_rbh
test_lrb_threshold
test_nbs
test_symmetrization
test_species_overlap
test_resolution_search
test_terminal_reason
```

### 29.2 Hierarchy invariants

每次 hierarchy 完成后必须自动验证：

1. children 互不重叠；
2. `union(children) = parent`；
3. 每个 protein 最终只属于一个 terminal family；
4. `child.n_genes < parent.n_genes`；
5. hierarchy 无环；
6. 每个非 root cluster 只有一个 parent；
7. 所有 component proteins 都被覆盖。

### 29.3 V1 regression dataset

至少建立：

- Dataset A：4 个小 proteomes，验证完整流程；
- Dataset B：明确 paralogs，验证 RBH/LRB；
- Dataset C：大型 gene expansion，验证 hierarchy；
- Dataset D：fusion/multidomain，验证 bridge sensitivity；
- Dataset E：大型 connected component，验证性能。

### 29.4 Regression 类型

区分：

- compatibility regression
- scientific regression

当 V2 结果与 V1 不同时，必须可以归因到明确改变，例如：

- bug fix
- coverage filtering
- edge symmetrization
- RH implementation correction
- new Leiden criterion

---

## 30. Performance Benchmark

每个 benchmark 至少记录：

```text
total runtime
peak RSS
edge count
connected component count
number of Leiden runs
number of generated subgraphs
hierarchy node count
terminal family count
disk read
disk write
```

尤其关注：

```text
Leiden calls
subgraph constructions
```

这两个指标直接反映 V1 的核心性能瓶颈是否真正解决。

---

## 31. 性能架构要求

### Requirement A

Protein 数增加时，hierarchy 内存不得按所有 ancestor memberships 累积。

### Requirement B

Hierarchy node 数增加时，不得产生对应数量的全局 SSN copies。

### Requirement C

Worker 之间不传完整 graph。

### Requirement D

任何可能超过 RAM 的输出都必须支持 streaming。

### Requirement E

任何主要 stage 都必须支持 component-level resume。

### Requirement F

可以用 integer arrays / columnar storage 时，不使用超大 `dict[str, ...]`。

---

## 32. 科学验证体系

未来正式方法验证至少比较：

```text
OGProfiler 1
OGProfiler 2
OrthoFinder
其他适用 orthology/family inference 方法
```

评价维度：

### Family-level

- family number
- family size distribution
- single-copy family recovery
- species coverage
- large-family fragmentation

### Stability

- different random seeds
- different gamma ranges
- different search backends
- different symmetrization rules

### Evolutionary consistency

- species overlap
- gene-tree concordance
- species-tree reconciliation

### Orthology

如果有 reference/simulation：

- precision
- recall
- F1

---

## 33. Synthetic Benchmark

建议建立已知真实 genealogy 的 synthetic 数据：

```text
known ancestral family
        │
    speciation
    /      \
   A        B
             \
           duplication
           /        \
          B1        B2
```

通过模拟不同：

- sequence divergence
- duplication rate
- gene loss
- domain fusion
- family expansion

研究 hierarchical Leiden 在哪些条件下能恢复真实结构。

这不仅用于软件测试，也可直接服务未来方法学论文。

---

## 34. OGProfiler 2.0 MVP 边界

OGProfiler 2.0 第一版建议只完成：

```text
FASTA
 ↓
homology search
 ↓
normalized edges
 ↓
connected components
 ↓
hierarchical Leiden
 ↓
network families
 ↓
hierarchy
 ↓
network event annotation
```

核心结果：

```text
families.tsv
members.tsv
hierarchy.tsv
events.tsv
```

以下内容建议进入 2.1：

- gene-tree refinement
- species-tree reconciliation
-大规模 pairwise ortholog generation
- subtree-level distributed scheduling

---

## 35. 最终总体数据流

```text
                  ┌──────────────┐
                  │   Proteomes  │
                  └──────┬───────┘
                         │
                         ▼
                ┌─────────────────┐
                │ Homology Search │
                └────────┬────────┘
                         │
                         ▼
                 directional hits
                         │
                         ▼
              ┌─────────────────────┐
              │ Similarity Engine   │
              │                     │
              │ coverage            │
              │ normalization       │
              │ RBH/LRB             │
              │ symmetrization      │
              └──────────┬──────────┘
                         │
                         ▼
                    edge table
                         │
                         ▼
                ┌────────────────┐
                │   Union-Find   │
                └───────┬────────┘
                        │
                  components
                        │
          ┌─────────────┼────────────┐
          ▼             ▼            ▼
        CC1            CC2          CC3
          │             │            │
          ▼             ▼            ▼
      Hierarchical   Hierarchical  Hierarchical
        Leiden         Leiden        Leiden
          │             │            │
          └─────────────┼────────────┘
                        │
                        ▼
                Hierarchy Forest
                        │
              ┌─────────┴─────────┐
              ▼                   ▼
       Terminal Families    Network Events
                                  │
                                  ▼
                            optional phylogeny
                                  │
                                  ▼
                         Evolutionary Events
```

---

## 36. 最重要的架构原则

OGProfiler 2 后续每增加一个功能，都应先问：

> **这个信息真的必须存在于内存中的 graph object 里吗？**

职责必须严格分离：

```text
Graph            → Leiden calculation
Hierarchy table  → parent-child relationships
Parquet/table    → persistent data storage
Species bitset   → overlap/event computation
Gene tree        → phylogenetic inference
Exporter         → GML/GraphML/TSV/FASTA
```

绝不再让 `igraph.Graph` 同时承担：

- SSN
- hierarchy
- gene membership database
- evolution events
- orthology navigation
- final storage format

---

## 37. 推荐的首个技术验证原型

在正式重写整个 pipeline 之前，优先实现一个最小 prototype：

```text
core/models.py
similarity/edges.py
graph/components.py
hierarchy/leiden.py
hierarchy/engine.py
storage/hierarchy.py
```

直接读取 V1 已有 SSN 或简单 edge table：

```text
edge table
   ↓
connected components
   ↓
component-level DFS Leiden
   ↓
hierarchy.parquet
   ↓
families.tsv
```

同一个 V1 `ssn.gml` 同时输入：

```text
OGProfiler 1 hierarchy engine
OGProfiler 2 hierarchy engine
```

比较：

```text
runtime
peak RSS
Leiden call count
subgraph construction count
hierarchy depth
family count
family membership
```

如果新的 component-centric hierarchy engine 能显著降低时间和内存，并保持合理 family structure，则 OGProfiler 2 的核心技术路线成立。

---

## 38. 开发优先级

最终优先级建议明确为：

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

OGProfiler 2 的成功标准不是“代码拆得更漂亮”，而是：

1. hierarchy engine 的计算规模显著提升；
2. network hierarchy 与 evolutionary interpretation 科学边界明确；
3. 数据模型可扩展、可恢复、可测试；
4. 每一个算法改变都可以通过 regression 与 benchmark 解释；
5. 最终形成一个可以持续演化，而不是再次退化为单文件脚本的工程基础。
