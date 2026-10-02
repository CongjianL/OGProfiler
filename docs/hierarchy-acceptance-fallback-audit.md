# hierarchy 候选准入与失败回退审计

日期：2026-10-02。范围：本地源码审计、P6 已保存 component 0 证据、轻量确定性复现。
后续语义与有限回退设计已形成 Proposed ADR 0002：
`adr/0002-hierarchy-admission-bounded-fallback.md`；尚未修改生产政策。
**未修改生产算法、未调阈值、未提交计算任务。**
新增 `tests/unit/test_hierarchy_acceptance_audit.py` 为当前行为的 characterization tests，
测试通过不表示准入策略已修复。

## 总结

P6 巨组不是 Leiden 没找到划分，而是 V2 的全局准入规则拒绝全部已测试划分，
随后将“搜索未接受”落为普通 terminal 完整成员；OG 的实际 V1 残余语义继承该巨组。
这是一条 **准入政策 → 失败回退语义 → 下游分组** 的链，不是蛋白丢失或 ID 映射错误。
其中严格 child size 规则是 ADR 0001 的既有设计，不是本次发现的条件写反。
若改为接纳 singleton 子群，应明确修订 ADR，而非仅声称修 bug。

## 1. 真实根节点证据

作业 `1410751`，已保存本地 `results/{historical-frozen,repaired-mean}-root/`。
两组 root 都只有一个 terminal 节点，无实际子节点；所有候选 valid=false / selected=false。

| 指标 | 历史 frozen | 修复 mean |
| --- | ---: | ---: |
| 蛋白 / 物种 | 69,701 / 12 | 69,642 / 12 |
| 候选 / Leiden 调用 | 10 / 30 | 10 / 30 |
| primary UNSTABLE 候选 | 4 | 1 |
| primary GAMMA_LIMIT 候选 | 6 | 9 |
| 搜索输出 terminal_reason | UNSTABLE | UNSTABLE |

修复 mean 候选完整摘要：

| gamma | child_count | minimum | mean ARI / stability | primary reason |
| ---: | ---: | ---: | ---: | --- |
| 0.01 | 148 | 2 | 0.834731 | UNSTABLE |
| 0.02 | 192 | 1 | 0.815493 | GAMMA_LIMIT |
| 0.04 | 247 | 1 | 0.838145 | GAMMA_LIMIT |
| 0.08 | 377 | 1 | 0.858242 | GAMMA_LIMIT |
| 0.16 | 496 | 1 | 0.904183 | GAMMA_LIMIT |
| 0.32 | 658 | 1 | 0.872424 | GAMMA_LIMIT |
| 0.64 | 895 | 1 | 0.853479 | GAMMA_LIMIT |
| 1.28 | 1,179 | 1 | 0.882533 | GAMMA_LIMIT |
| 2.56 | 1,559 | 1 | 0.893638 | GAMMA_LIMIT |
| 5.12 | 1,971 | 1 | 0.905719 | GAMMA_LIMIT |

0.16 与 5.12 的 stability 已过 0.9，仍因 minimum=1 整体被拒绝。
历史 gamma=0.08 的 stability=0.899918，刚低于阈值；后续若干 gamma 稳定性达标，
但有 singleton。两组都不应归纳为“全部划分不稳定”。

## 2. 候选准入：一个 singleton 否决整份 partition [高优先级]

`hierarchy/resolution.py:_candidate`，185–192：

```text
child_count <= 1                            → NO_SPLIT
minimum < min_child_size                    → GAMMA_LIMIT
  或 largest fraction > max_child_fraction
  或 tiny fraction > max_tiny_fragment_fraction
stability < stability_threshold             → UNSTABLE
quality < min_quality                       → LOW_QUALITY
否则 valid=True
```

`min_family_size=2` 被 CLI 映射为 min_child_size=2；任何一个 singleton 就使整份
partition 失效，而非只把该 singleton 终止成叶。默认 max_tiny_fragment_fraction=1
并未抵消这个独立 minimum gate。
轻量 fixture `[1,49,50]`：stability=1、最大子群比例0.5、tiny fraction=0.01，
仍因 minimum=1 被拒绝。若同时 stability=0.1，primary reason 仍只写 GAMMA_LIMIT。

V2 已有 singleton 叶处理能力：`engine._terminal_reason` 的第一个分支为 SINGLETON，
已接受 partition 的每个非空社区可以完整写为子节点。
目前拦截发生在创建子节点之前，是准入政策，而非存储模型不支持 singleton。

### 与 V1 的明确差异

冻结参考 `legacy/OGProfiler_v1.py` SHA256 为
`35607fb7cdbd5f8db6c6a3bc3f22df953884b24af2df023f79307bdfcdc4ac3b`。
实际 `HHN.RunCommunityDetection`：

- 多物种 >=10,000 蛋白：尝试目标20–50社区，返回结果有 >=2 社区即建立子节点。
- 多物种 3–9,999 蛋白：尝试二分，返回恰好2社区即建立子节点。
- 多物种2蛋白：直接建立两个 singleton 子节点。
- 不检查全部子群 minimum>=2，不使用三 seed mean ARI>=0.9 准入。

新增 fixture 执行实际固定 V1 函数，仅 stub 昂贵 partition 搜索，证明其接受
`1+2`、`1+3` 子群；不是声称 V1 在相同真实 SSN 必然生成这些划分。
V1 自身也可能因返回社区数不符合条件整体终止；它并非永远成功的回退。

## 3. 一个参数耦合两种语义 [中优先级]

`engine._terminal_reason`：`len(node) < 2 * min_child_size` 直接 MIN_SIZE，搜索不执行。
默认2使两/三蛋白多物种组件直接整体终止；V1 两蛋白强制分叶，三蛋白可以 `1+2`。
因此 min_family_size 当前同时表达：

1. 当前节点是否值得尝试分裂；
2. 每个输出子群的最小尺寸。

建议后续区分“递归停止尺寸”和“合法 partition 非空子群”，而非仅全局改2→1。
这会涉及 ADR 0001 的小组件根终止和 strict minimum 规则，需显式设计决策。

## 4. adaptive 不是区间内可行性的证明 [中优先级]

`search_resolution`，246–259：指数步进，只在先找到有效 coarse candidate 后
才在前一个失败点和该有效点之间做 local grid。

- 当前配置实际测试 `.01,.02,...,5.12`；下一点10.24越界，所以 **gamma_max=10 从未测试**。
- 所有 coarse 点被拒绝时，local_grid_points=5 完全未启用。
- 有效性不保证随 gamma 单调：singleton、ARI波动等会产生狭窄可行区间。
- toy fixture 使只有 gamma=10 有效，实际搜索 selected=None；另一个 fixture 的
  `(1.3,1.7)` 可行，但 coarse 仅1、2，也 selected=None。

这些是确定性搜索覆盖缺口；**未证明真实根在gamma=10或漏采区间一定存在好划分**。
因此 GAMMA_LIMIT/UNSTABLE 最多说明已测试候选未接受，不等于整个区间无可分性。
未来可单独处理端点和有预算的失败补采，不把“无限重试 Leiden”当作修复。

## 5. 稳定性是全局平均 ARI，不是逐子群置信度

`_multi_seed` 使用 seed、seed+104729、seed+209458；各对 ARI 取平均作为 stability，
选择与其他 partition ARI 支持和最高的代表（quality 用于 tie）。NMI 仅诊断。
fast 模式只有一个 seed，stability 固定1；默认CLI robust模式才执行三次。

该指标不能指出具体哪些子群不稳定；当前一个全局值否决全部社区。
代表 partition 的稳定性分数还是三者的平均，而非只计算该代表到其他结果的平均。
ARI 公式的当前实现并非已证明的错误；问题是科学准入语义和粒度。
切换fast仅绕开可靠性检查，不是经验证的准确度修复，且会改变候选/层级。

## 6. 回退将搜索失败转换成普通完整 terminal [高优先级]

serial `engine.infer_component_hierarchy` 的 selected=None 分支：

```text
split_status=TERMINAL
terminal_reason=search.terminal_reason
全部 global_ids → 当前 cluster_id
```

parallel `subtree._evaluate` 也返回无children且有terminal_reason，汇总时全体成员指向该节点。
没有“接受最佳拒绝划分”的回退，没有细分 unresolved 和已支持 terminal，
没有超大多物种失败节点的额外门槛。层级验证只验证结构/成员分区，不检验生物学置信度。
所以全部蛋白覆盖率100%与零重复仍能伴随非常低的 precision。

下游专用 V1事件：degree-zero且n_genes>n_species → None。
`orthogroups/engine.py` level=1 分支只筛选 Event='None'，无单物种条件、无大小限制。
完整真实 V1参考也执行相同行为，所以 P5/P6 parity 通过却保留巨组。
toy 6蛋白/2物种，唯一划分含singleton而失败，结构校验通过，随后生成一个完整
RESIDUAL_NONE OG；复现这条链无需真实大图。

**不建议在 OG engine 临时加大小过滤或忽略多物种 None：那会破坏已批准 V1契约，
只掩盖上游失败，且改变未分配集合。**

## 7. 诊断把不同失败原因混为一谈 [中优先级]

候选只保留首个失败条件；GAMMA_LIMIT可指child size、dominance或tiny fraction，
不一定是gamma上界耗尽。搜索汇总reason优先 UNSTABLE，只要有一个UNSTABLE，
即使其余9个因minimum=1失败，root也只写UNSTABLE。
fixture明确复现 `UNSTABLE + GAMMA_LIMIT → terminal UNSTABLE`。

建议未来保留所有 violations 与 raw measurements，并分别记录搜索覆盖/预算耗尽、
合法但可靠性不足的partition、真正NO_SPLIT，以及后续采用的fallback策略。
这属于诊断/schema改动，应绑定新算法/产物版本，不能偷偷复用旧 hierarchy cache。

## 8. V1 参数体系也不同，不能直接移植 gamma

V1 >=1000节点采用 `MappingGammaForCC` 的按大小缩放起点，再以1.5倍探索；
`BipartiteGraphs` 用二分寻找社区数范围，Leiden `n_iterations=10`，没有seed固定。
以 coefficient=1、69,642蛋白为例起点约0.00014359；V2 gamma_min=0.01，
本身已产生148社区，且不存在V1的20–50社区约束。
V2 `run_leiden` 不显式传n_iterations，不能把两套运行当作同一partition。
V1 SpecificBoard 有非严格有界搜索，直接移植会违背当前有预算执行需求。
这些是算法/配置策略差异，不是本审计已经在同SSN上复现的性能或准确度优势。
P6完整首版本 V1的源码又不同于本节固定参考，需保持版本标签清楚。

## 9. 后续设计建议及验证顺序（未执行）

1. **先明确政策/ADR**：非空singleton子群可保留并叶终止；将递归停止与
   partition准入解耦；稳定性失败、预算耗尽不直接等价为已支持的大多物种终端OG。
2. **诊断完整化**：所有失败条件、已测试gamma、端点/补采、fallback来源；
   保留质量/ARI/NMI，失效旧缓存与schema。
3. **有预算的失败策略**：端点/局部补采；若仍失败，显式 unresolved 或组件失败，
   不静默把拒绝划分当作可靠划分；如允许可靠性放宽应有单独策略/标记与证据。
4. **固定SSN，先组件0小范围Slurm验证**：比较singleton政策、搜索覆盖、
   fallback各自贡献；不同时更改OG/scorer，避免直接全量扫描参数。
5. 再做结构/成员/parity/resume与正式Orthobench，观察P/R/F1和大组污染，
   而非以root终于SPLIT作为准确度验收。

本轮验证：9个新增characterization tests + 8个既有hierarchy tests，共17 passed。
全套本地测试285 passed（显式提供实际OF/Orthobench参考源码）；Ruff与diff检查通过。
生产源码/科学配置未变，既有P5/P6快照不受影响。
