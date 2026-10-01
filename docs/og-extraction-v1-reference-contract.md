# OG 提取 P0/P1：实际 V1 参考契约

日期：2026-10-01（Asia/Shanghai）。本轮完成参考冻结、专用事件视图和单组件纯选择引擎，
尚未接阶段落盘或最终 export，SSN/hierarchy/network_event 保持原样。

## 参考与可复现测试

- 文件：`legacy/OGProfiler_v1.py`。
- SHA256：`35607fb7cdbd5f8db6c6a3bc3f22df953884b24af2df023f79307bdfcdc4ac3b`。
- `benchmarks/og_extraction/reference_v1.py` 从 AST 取原始 HHN 构造器、
  GetEvolutionEvents、ExtractOG、RunCommunityDetection 和实际 ExtractOGSorted/GetGenesID/GetGenesIDs。
  决策与遍历函数保持原样；进度条仅替换为无输出适配器。GML 写入只在临时目录发生。
- hash 改变触发参考校验失败，需重新核对和冻结。
- 13 个有效/异常 fixture 的参考事件、组、处理等级、消耗节点、剩余节点与 degree、
  未分配和重复成员保存在 `tests/fixtures/og_extraction/v1_contract.json`。
  expected 来自实际 V1 执行，不由新实现生成。组内成员排序用于消除 set 输出顺序。
- `tests/unit/test_og_v1_reference.py` 逐次参考执行与 frozen expected 比较。
- `tests/unit/test_og_v1_events.py` 将同一 fixture 转为组件 nodes/members，
  在整数 overlap=0/1/10 下与实际 V1 原始标签差分；另测试类型校验、数据损坏和 3001 节点深树。

运行：

```bash
.venv/bin/python -m pytest -q tests/unit/test_og_v1_reference.py tests/unit/test_og_v1_events.py tests/unit/test_og_v1_engine.py
```

## 已确认的事件/拓扑语义

1. 实际入口不是所有有 children 的节点：V1 只注释无向 degree>1 的节点，
   候选邻居还必须满足 `neighbor.genesNum < node.genesNum`。
2. degree=0 且 genesNum=genomesNum 时赋 I；degree=1 不赋事件，保持 Python None。
3. II 优先于 disjoint/允许重叠的 I；overlap 是整数物种数。多于两个 eligible 邻居是 III-3。
4. `ConstructedHnCC` 不创建 Event 列；任一 GetEvolutionEvents 赋值会在全局图建立该列，
   其他节点自动是 Python None。若全图没有任何赋值，后续 hnn_analysis 读取 Event 会 KeyError。
5. 组件局部兼容视图显式初始化 None，以模拟“全局列存在而该节点未赋值”的常见上下文。
   缺列错误另有 fixture，不把这个初始化称作原始失败输入上的逐位等价。
6. 未精炼入口在选择前执行 Event 字符串归一化：None→'None'。
   GetGenesID 对 Python None 和字符串 'None' 的度数>1节点遍历不同，已分别测试。
7. V1 对两基因/两物种的组件强制两个 singleton children；V2 允许根直接终止。
   本例根的 OG 事件/最终成员等价（V1 二叉根 I；终端 V2 根 degree=0 且基因数=物种数也为 I）。
   专用事件视图不持久化虚拟新层级，也不更改 ADR。

## 已确认的选择/删除语义（作为 P2 输入）

- GenomeNum 从数据集总物种数向下迭代，级别优先于 genesNum；同级 genesNum 降序，
  排序稳定、tie 保留原参考 vertex 次序。reference_order 显式保留节点输入顺序。
- EVENT_I 通过 GetGenesIDs 遍历严格小于当前原始 genesNum 的邻居；
  原事件已是字符串，只有 active degree<=1 才直接读该节点原始 geneIDs。
- 只删除遍历收集的 deletedVertex，选中节点仍保留；节点事件和原始 counts 不重新计算。
  固定 parent_II_child_I 例最终保留 parent 和两个已选 child，均由原属性解释。
- GenomeNum=1 选择 Event='None'，不检查真实物种数。多叉例中的残余 OG 都有两个物种。
- 残余读取原 geneIDs 而不是 active subtree 新统计，且不删除节点。
- GenomeNum=0 从 SSN degree=0 蛋白补充；hierarchy singleton 根若已标为 I，
  不在级别1残余分支选取。fixture 验证 isolate 最终只由 SSN 分支输出一次。
- 当前默认 max_num=0，不添加大小过滤。OG 编号和组内 set 顺序不作为集合差分条件。

## 已分类的异常，不用简化规则隐藏

| 边界 | 实际 V1 行为 | 本轮兼容处理 |
| --- | --- | --- |
| 非根 unary，degree=2 但 eligible 邻居<2 | GetEvolutionEvents IndexError | 专用事件视图抛带节点形态说明的 HierarchyError |
| 等成员 unary 根/叶，各 degree=1 | 两节点均 None；残余输出同一成员两次 | oracle 明确记录 duplicate_members；P2 必须显式报告冲突 |
| 全局 Event 从未赋值 | hnn_analysis 读取列时 KeyError | oracle 保留错误测试；组件局部视图明确 seed None 上下文 |

上述分类不是重新解释为 I、POLYTOMY，也不以静默去重通过验收。
正常严格分裂 hierarchy 的计数/parent/member 校验先于事件计算。

## P1 实现与边界

新增 `src/ogprofiler/orthogroups/models.py` 和 `legacy_events.py`。
`annotate_v1_events` 使用显式 parent/children 推导原始无向 degree，迭代 postorder
计算 bitmap 和 n_genes，对节点统计、成员唯一性、根、parent、连通性和 depth 进行校验。
仅保留 O(nodes+members) 数据，不持久化内部 descendant protein lists，也不构建全局 igraph。

返回专用 V1EventAnnotation：v1_event 为标签或 Python None，selection_event 通过显式属性
转换为选择阶段使用的标签/'None'。这与 network_event/legacy_event 分离，目前没有改其 schema/缓存。

## P2 单组件引擎

`src/ogprofiler/orthogroups/engine.py::extract_component_orthogroups` 实现 V1 未精炼 I/0/0 策略：

- 从数据集 total_species 向下处理；原始 Event 与 counts 保持固定。
- 按节点输入/reference 顺序生成候选，稳定按原始 genesNum 降序排序。
- active 无向 degree 决定遍历还是读取原始成员区间；整批覆盖等级结束后才 deactive descendants。
  选中节点仍保留，残余分支只接受归一化后的 'None'，不强制补齐未分配蛋白。
- SSN isolate 由显式参数补充，不从 hierarchy singleton 名称猜测。
- 原始成员采用一次 DFS 展开加节点区间索引；仅在真正选组时枚举，不给每个内部节点缓存后代列表。
- 返回 typed groups/trace/unassigned/remaining_cluster_ids。成员 hash 使用排序后的
  `(species_id, original_id)`，跨物种同名蛋白保持不同身份，protein_id 重编号不改变该 hash。
- 出现重复分配时抛 OrthogroupConflictError，携带诊断 result 和 duplicate_members。
  诊断候选保留实际 V1 重叠以解释分歧，但不作为正常生产 OG 返回，也不静默去重。

13 个 frozen fixture 中正常情况的 OG 成员、处理等级、剩余节点和未分配集合与 V1 一致；
等成员 unary 的重叠例得到显式冲突，诊断候选/重复集合仍匹配参考。
另验证父 I 合并多个 terminal family、已消耗候选的来源 trace、无 SSN isolate 证据的 singleton
保持 unassigned、成员 hash 身份、函数输入不变与重跑确定性，以及 3001 节点深层选择。

本轮专用测试 83 passed；完整测试集 181 passed；Ruff 和新包 Mypy 通过。
P3–P6 的配置、Parquet、resume、export、真实 hierarchy 回归及准确度仍待开发/验证。
