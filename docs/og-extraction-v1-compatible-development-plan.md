# OG 提取：V1 策略与 V2 架构的继续开发方案

> 状态：**用户已批准；P0–P2 已完成；P3 接线及本地落盘恢复测试已完成，全局 ID 交叉验收待 P4；P4–P6 待开发/验证**。批准日期：2026-10-01（Asia/Shanghai）。
> 本文是后续开发的执行依据，不表示功能或科学验证已经完成。
> 制定时仅新增文档；后续实现进度见第 9 节及 `og-extraction-v1-reference-contract.md`。

## 1. 目标与范围

采用 **V1 决策语义 + V2 组件化实现**，将 OG 提取从终端 family 导出中分离。

目标：

1. 恢复 V1 根据事件及物种覆盖选择内部节点、收集后代蛋白的 OG 策略。
2. 保留 V2 的组件局部处理、明确数据模型、Parquet 存储、校验与 resume。
3. 对同一套层级输入，以 V1 实际函数为参考验证 OG 成员集合一致性。
4. 再在固定 SSN/hierarchy 上验证真实数据，最终单独评价 Orthobench 准确度。

本轮不扩展基因树 refinement，不同时调整 SSN、Leiden 或层级分裂规则。
V1 语义兼容与最终 OG precision 改善是两项独立验收，分别提供证据。

## 2. 已知背景与当前问题

### 2.1 上游 SSN

此前真实蛋白组验证作业 `1410705` 已通过：13 物种、128,483 蛋白；
V2 方向 SSN 与实际 OF3.1.5 原始函数在容差内一致，CLI 落盘权重匹配
完整方向 W 的 `(Wforward + Wreverse) / 2`，最终 701,947 条无向边。
这支持先固定 SSN 排查 OG 提取；它不证明 hierarchy 或 OG 结果等价。

历史结果文档：
`/Users/licongjian/Desktop/PycharmProjects/OGProfiler/docs/ssn-repair-real-proteomes-validation.md`。

### 2.2 终端 family 不等于 OG

当前 V2 的 hierarchy membership 表示：

```text
protein_id → terminal_cluster_id
```

当前 export 将终端节点直接组装成最终 OG。终止位置可能来自大小、深度、
稳定性或分裂失败，而不是生物学上的 OG 选择。因此可能出现过度拆分或大节点整体输出。

V1 的最终 OG 可以由内部节点的多个后代 family 合并而来。
本轮将这两类对象显式分开：

- **terminal family**：hierarchy 的基础蛋白分区，保留现有不变量。
- **orthogroup**：经过独立策略选择得到的成员集合，可对应内部节点。

不得通过改变终端 membership 的含义来隐式实现 OG 提取。

## 3. 参考源码与证据优先级

V1 参考：
`/Users/licongjian/Desktop/PycharmProjects/OGProfiler/legacy/OGProfiler_v1.py`。

重点函数：

- `HHN.GetEvolutionEvents`：事件判定。
- `HHN.ExtractOG`：单个物种覆盖等级下的候选选择与删除行为。
- `ExtractOGSorted`：物种覆盖等级的处理顺序。
- `GetGenesID` / `GetGenesIDs`：实际后代遍历及成员收集。
- `Processes.hnn_analysis` 所在流程：未精炼路径调用 `ExtractOGSorted(..., 'I', 0, 0)`。
- `WriteOGFiles`：成员输出及输出编号。

证据优先级：**实际源码行为 > 固定输入参考执行 > 文档文字推断**。
后续实现前记录参考文件 SHA256。参考适配器只放在测试/benchmark，
生产代码不导入整个旧脚本，不依赖其命令行、全局网络或输出副作用。

当前 V2 相关位置：

- `/Users/licongjian/Desktop/PycharmProjects/OGProfiler/src/ogprofiler/evolution/network.py`
- `/Users/licongjian/Desktop/PycharmProjects/OGProfiler/src/ogprofiler/evolution/stage.py`
- `/Users/licongjian/Desktop/PycharmProjects/OGProfiler/src/ogprofiler/output/results.py`
- `/Users/licongjian/Desktop/PycharmProjects/OGProfiler/src/ogprofiler/output/stage.py`
- `/Users/licongjian/Desktop/PycharmProjects/OGProfiler/src/ogprofiler/orthology/engine.py`
- `/Users/licongjian/Desktop/PycharmProjects/OGProfiler/src/ogprofiler/orthology/stage.py`

## 4. V1 事件兼容契约

### 4.1 新增专用事件视图

新增 `v1_event` 供 OG 选择；保留现有 `network_event` 作诊断。
现有 `legacy_event` 映射尚未完整包含 V1 整数 overlap 阈值及特殊节点行为，
不作为本轮已验证的参考输入。

设父节点物种集合为 S，两个子节点为 A、B，整数参数为 overlap_count：

| 按顺序判断的条件 | V1 事件 |
| --- | --- |
| 子节点数大于 2 | `III-3` |
| `A ∩ B = S` | `II` |
| `A ∩ B` 为空 | `I` |
| `|A ∩ B| <= overlap_count` 且 `|A ∩ B| < |S|` | `I` |
| 剩余情形中 `A = S` 或 `B = S` | `III-1` |
| 其余 | `III-2` |

V1 默认 overlap_count 为 **0**，表示物种数，而非归一化比例。
判断顺序属于算法契约，尤其应先判 `II` 再判允许重叠的 `I`。

### 4.2 特殊节点

- V1 度为 0 的层级节点，若 `n_genes == n_species`，事件赋为 `I`。
- 其他叶节点/根节点的事件初始化及未赋值状态，依据实际 V1 初始化流程核对。
- 明确区分 Python `None`、字符串 `'None'` 和“不适用”；不要依赖字符串强转隐式判定。
- 单子节点、根的无向度、已删除后代后的度数，以及同成员数父子节点应有专门 fixture。
- V2 的显式 parent/children 不应直接等同于 V1 无向图的 degree 条件。

**实现前门槛：**先建立节点形态与事件初始化的兼容映射，测试后再接选择算法。
多叉节点 `III-3` 不因当前 V2 标记为 `POLYTOMY` 就自动成为 OG。

## 5. V1 OG 选择契约

### 5.1 未精炼默认路径

复现 V1 的实际顺序：

```text
输入物种总数 N
→ GenomeNum=N, N-1, ..., 2
   → Event=I 且 genomesNum=GenomeNum 的有效候选
   → genesNum 从大到小
   → 收集实际遍历得到的成员
   → 记录被消耗的后代、跳过已消耗候选
→ GenomeNum=1：处理剩余 Event='None' 节点
→ GenomeNum=0：补充 SSN 孤立蛋白
```

默认 `max_num=0`；暂不新增候选大小过滤，避免改变 V1 默认语义。

重要边界：

1. `GenomeNum=1` 的 V1 分支只选择 `Event='None'`，没有显式过滤单物种。
   将其称为残余处理，不把它重新定义为“所有单物种节点”。
2. 选中一个节点后的成员收集、后代删除及选中节点保留，分别建模。
   “把整个 subtree 永久删除”的简化实现需先由参考测试证明等价。
3. 事件在选择前建立；不在每次逻辑删除后自动重新注释事件。
4. 原始统计属性与 active 视图下的邻接/degree 分开存储。
5. 同物种覆盖、同蛋白数的排序 tie，保留可追踪的参考节点顺序；
   如使用新的确定性 tie-break，需证明在 fixture 上不改变成员集合，并记录其契约。
6. 参考实现中的 set 输出顺序和 OG 数字编号不作为集合等价验收条件。

### 5.2 剩余蛋白与重复分配

- 选出的 OG 成员不得重复分配；如参考执行发生重叠，报告并分类，避免静默去重掩盖问题。
- 未分配蛋白单独输出 `unassigned`，并附原因或关联剩余节点。
- 暂不通过“剩余蛋白强制各自成为 OG”保证覆盖率，这会掩盖策略差异。
- SSN 孤点按 V1 单独分支补充，核查与 hierarchy singleton 的去重。
- 原始 terminal membership 仍保持全蛋白分区，与 OG 的覆盖率统计分开。

## 6. V2 实现架构

新增包，目标位置：
`/Users/licongjian/Desktop/PycharmProjects/OGProfiler/src/ogprofiler/orthogroups/`。

```text
models.py          typed records：组、选择类型、trace 状态
legacy_events.py   V1 事件兼容函数；物种 bitmap 输入
engine.py          单组件纯选择引擎，无文件/调度副作用
stage.py           输入校验、组件执行、Parquet 落盘、manifest/resume
```

建议接口职责：

```text
annotate_v1_events(component_hierarchy, species_metadata, overlap_count)
extract_component_orthogroups(hierarchy_view, terminal_membership, v1_events, config)
run_orthogroup_stage(run_root, config, command)
```

### 6.1 存储与内存

- 按组件读取 nodes/members，构建 parent、children、postorder、物种 bitmap。
- 自底向上计算 n_genes/n_species，并与节点表校验。
- 用 active/consumed/selected 状态模拟 V1 的原地删除；层级文件始终保持原样。
- 不创建全局 igraph HHN，不在内部节点持久化完整 descendant protein lists。
- 使用 DFS 区间或等价索引查后代，真正选中时才枚举成员并分批写出。
- 避免缓存每个内部节点的完整成员集合，防止深层树导致内存倍增。
- 组件内部保持确定性顺序；调度器管理组件并发，engine 不创建嵌套进程池。

### 6.2 与既有 ADR 的关系

沿用组件化 ADR：
`/Users/licongjian/Desktop/PycharmProjects/OGProfiler/docs/adr/0001-component-centric-hierarchy-architecture.md`。

V1 强制 singleton leaves 与 V2 根直接终止的差异，用**兼容视图**显式处理，
不改写现有 terminal membership 或凭空持久化新 hierarchy 拓扑。
若适配确实需要改变 ADR 的层级不变量，先提交设计决策，不与 OG engine 实现混合。

## 7. 产物与输出接线

建议增加组件级产物：

```text
orthogroups/components/component=.../
    groups.parquet
    members.parquet
    selection_trace.parquet
    unassigned.parquet
    og-manifest.json
```

组记录至少包括：

```text
component_id
local_group_id
source_cluster_id       # SSN isolate 可为空
selection_type          # EVENT_I / RESIDUAL_NONE / SSN_ISOLATE
v1_event
n_genes
n_species
membership_hash
```

成员记录：`component_id, local_group_id, protein_id`。
trace 至少记录：候选节点、处理等级/顺序、事件、选择/跳过状态、消耗来源。
大规模数据下 trace 应有明确粒度，不生成每个蛋白对的 trace。

成员 hash 使用排序后的 `(species_id, original_id)`，避免跨物种同名蛋白混淆。
全局 OG ID 在最终导出时按成员集合确定，不依赖 worker 完成顺序。

新流程：

```text
SSN → hierarchy → evolution → orthogroups → export
```

- export 的最终 OG 表及 FASTA 改为读取新的 OG membership。
- terminal family 继续作为单独诊断对象；避免同一名称代表两种语义。
- 历史 terminal 导出保留为显式诊断/兼容策略，不作为新默认 OG。
- export 类型名若仍为 TerminalFamily，应拆分或泛化，避免把内部节点 OG 强塞入旧类型。
- pairwise ortholog 的事件遍历与 OG 分组是不同产物。先明确其阶段依赖与策略版本，
  不自动把每个 OG 的所有跨物种成员对都当成 ortholog。

## 8. 配置与缓存

批准的配置方向：

```yaml
orthogroups:
  strategy: v1_compatible
  species_overlap_count: 0
  refinement: false
```

默认使用 v1_compatible；策略名称、schema、CLI 参数需集中定义及校验。
初期 refinement=true 应显式提示该策略尚未实现，避免静默忽略。

OG manifest 绑定：

- 算法版本、策略、整数 overlap 阈值、refinement 状态；
- hierarchy nodes/members、蛋白 metadata、singleton、事件输入校验和；
- V1 兼容事件算法版本；
- 输出校验和及 selected/unassigned/duplicate 等统计。

缓存迁移：

1. 新建 OG 阶段算法版本，逐组件验证 resume。
2. export 版本升级，显式绑定 OG 输入。
3. SSN/hierarchy 输入校验一致时复用，避免重新搜索或聚类造成混杂。
4. changed strategy / overlap / membership / artifact corruption 均触发对应重建。
5. 采用临时文件 + 原子替换；失败状态与完成 manifest 区分。

## 9. 开发顺序及阶段验收

### P0：冻结参考契约

- [x] 记录 V1 文件 hash，梳理初始化、degree、None 字符串及残余分支。
- [x] 制作最小 hierarchy fixtures，记录实际 V1 函数输出。
- [x] 核实新旧根/叶拓扑的兼容视图。
- [x] 明确排序 tie 和逻辑删除的可观察行为。

**退出条件：**参考输出可重复，边界无未解释分歧。

### P1：事件兼容

- [x] 实现专用 v1_event 与整数参数校验。
- [x] 覆盖 I/II/III-1/III-2/III-3、部分重叠、阈值边界、独立节点。
- [x] 将专用结果与 network_event 保持分离。

**退出条件：**同输入参考事件标签一致。

### P2：单组件 OG 引擎

- [x] 实现覆盖等级、候选排序、后代收集、逻辑消耗、残余和孤点处理。
- [x] 输出组与 trace，报告重复/未分配状态。
- [x] 测试纯函数确定性及深层/多叉结构。

**退出条件：**同层级 V1/V2 OG 成员集合一致，重复分配为零；未分配集合与参考一致或已分类。

### P3：落盘、配置与 resume

- [x] 实现组件级 Parquet、manifest、原子写入及调度接线。
- [x] 新增配置段，默认策略 v1_compatible。
- [x] 验证参数改变、输入改变、损坏、部分失败与恢复。

**退出条件：**落盘重读结果与 engine 一致；串行/并发输出成员及全局 ID 一致。

P3 已验证 local_group_id、membership_hash 和全部组件产物一致；最终全局 OG 数字 ID
依第 7 节在 P4 export 分配并验收，该项跨阶段退出条件仍待验证；本阶段不提前引入第二套全局编号。
实现/恢复契约见 `og-extraction-stage.md`。

### P4：export 与下游

- [ ] 最终 OG 输出读取新 membership；终端 family 单独保留。
- [ ] 稳定 OG ID、TSV、FASTA 与统计表相互一致。
- [ ] 更新 export cache、CLI --until-stage 与 pairwise 依赖契约。
- [ ] 更新用户文档，区分 OG 与 terminal family。

**退出条件：**实际 CLI 产物通过端到端 fixture 检查，旧导出缓存失效。

### P5：固定 hierarchy 的真实数据回归

- [ ] 固定已经验证的 SSN 与同一份 hierarchy，不重搜、不同时调分裂参数。
- [ ] V1 参考与新引擎输入同一结构，按原始蛋白成员集合比较。
- [ ] 输出首个分歧及分类，记录组大小/物种覆盖/异常大 OG。
- [ ] 检查性能与内存，证明按组件执行没有全树成员膨胀。

**退出条件：**参考策略差分通过，真实产物校验通过；此时才进入准确度评价。

### P6：Orthobench 准确度评价

- [ ] 同 SSN/hierarchy 下比较 terminal 策略与 v1_compatible，隔离提取策略贡献。
- [ ] 分别报告 precision、recall、F1、覆盖率及过度合并/拆分案例。
- [ ] 区分完整 V1 跑法与“V1 策略作用于 V2 hierarchy”的差异。
- [ ] 如仍有偏差，先按 hierarchy/提取/评估映射分类，不同时修改多个层。

**退出条件：**依据预先确定的评估指标作科学判断，不把进程成功等同于准确度达标。

## 10. 必需测试矩阵

| 场景 | 核验重点 |
| --- | --- |
| 父 I，多终端后代 | 合并为一个父节点 OG，而非每个终端各一个 |
| 父 II，子 I | 父跳过，分别选择子节点 |
| 同覆盖数嵌套候选 | 大节点优先，后代候选被消耗 |
| 覆盖数不同、大小顺序相反 | 物种覆盖优先于组大小 |
| 部分物种重叠 | 整数阈值边界及判断顺序 |
| III-3 多叉 | 不自动接受父节点，检查下层及残余 |
| degree=0 / root terminal | 特殊事件、V1/V2 拓扑兼容 |
| None 与 'None' | 初始化、分支、序列化语义一致 |
| SSN isolate 与 hierarchy singleton | 无重复补充 |
| 原始 ID 跨物种同名 | 不错误合并身份 |
| 同排序键候选 | 确定性及成员集合等价 |
| 串行/并发、resume、损坏输入 | 产物与全局编号一致，正确失效 |

差分必须使用独立 V1 参考或已冻结的期望值，避免让新实现生成自己的 expected。
第一处分歧应包含节点成员 hash、事件、候选顺序和消耗链，而不仅提供 OG 数量差。

## 11. 远端执行规则

后续涉及真实数据时遵守项目 AGENTS.md：本地开发和验证，显式部署；
重计算通过 Slurm、不可变源码快照、输入/参数 hash 和 RUN_ID/JOB_ID。
复用既有正确 SSN/hierarchy 时先核查来源和校验和，不覆盖历史科学结果。
用户只要求实现或准备脚本时，不自动启动大规模验证。

## 12. 给下一位开发者的开始指令

**从 P0 开始，不直接改 export，也不直接把 legacy_event='I' 当作已验证 OG 规则。**
先用实际 V1 函数建立最小参考测试，验证事件及逻辑删除语义，再实现单组件引擎。
保持 SSN/hierarchy 固定，逐阶段提交可独立验收的变更。
批准范围是恢复 V1 提取策略及 V2 风格实现，不包括同时优化生物学判定或调参。
