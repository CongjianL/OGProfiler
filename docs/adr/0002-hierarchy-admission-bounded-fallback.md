# ADR 0002：非空子群准入、递归停止与有预算搜索失败

- Status：Accepted（2026-10-02 用户确认实施；H1–H3 本地实现，H4/H5 待真实数据验证）
- Date：2026-10-02
- Scope：hierarchy 准入、搜索、失败状态与输入验收；OG 决策算法与 scorer 保持原样。
- Amends：ADR 0001 的 Splitting and resolution search 第4条、Terminal reasons 与根终止语义。
- Evidence：`hierarchy-acceptance-fallback-audit.md`、P6 job `1410751`。

> 2026-10-02 后续修订：ADR 0003 将新策略每次 Leiden 的默认有限迭代预算从2
> 提高到10；下文2轮为首轮实施历史记录。准入、停止尺寸、稳定性阈值与候选预算不变。

## 1. 背景与明确的决策边界

修复 SSN 上69,642蛋白的根已有148–1,971社区候选；全部因稳定性或singleton
准入被拒绝。根整体terminal随后被原样RESIDUAL_NONE输出，并非OG算法新合并。
旧 strict minimum 是 ADR 0001 的设计选择。本稿明确修订该政策，不把它描述为条件写反。

保持组件局部存储、非空分区、完整成员保留、确定性、serial/parallel等价。
不恢复全局HHN、不复制内部节点完整成员列表、不切换MCL、不更改SSN或Leiden目标函数。
不在OG引擎添加大小过滤、多物种None过滤或强制补singleton；不修改官方scorer/ID映射。

## 2. 三个独立概念

### A. partition准入：允许singleton，拒绝空子群

**结构合法性**：每个子群非空、>=2个子群、成员互斥、并集恰为父节点、ID均来自父节点。
singleton是合法子群，创建真实子节点并终止为SINGLETON；不为了拓扑兼容凭空造叶。
严禁把singleton丢弃、强塞入最大社区、或因为一个singleton拒绝全部社区。

**科学/政策准入**：结构合法后，继续检查最大子群比例、tiny fragment总比例、
稳定性与可选quality门槛；membership与quality须有效、有限。
每个不满足条件都记录在violations中，不仅保留第一个条件。
这允许明确配置的tiny fragment总比例拒绝候选，但取消隐藏的“任何singleton一票否决”。

稳定性定义本轮不变：robust三seed的mean pairwise ARI，阈值0.9；
NMI仅诊断。fast仅作为显式配置的single-seed模式，并记录stability_evaluated=false，
不得把固定1当成测得的稳定性，正式本轮验证仍用robust。

### B. 递归停止尺寸：只控制某个节点继续搜索

新增 `recursion_stop_size`，正整数，拟议默认 **1**。
节点大小 `n <= recursion_stop_size` 时不运行Leiden，形成政策终止叶：
singleton优先标SINGLETON；其他标SIZE_STOP。
该参数不约束兄弟子群的minimum，不使用 `n < 2 * min_child_size` 判断。
小节点终止不阻止其兄弟继续递归。
较大stop_size是用户显式的粒度政策，可能造成多物种粗组，不是可靠性证明。

ONE_SPECIES保持现有停止语义。MAX_DEPTH对于仍待分裂的多物种节点视为资源/策略
限制导致UNRESOLVED；NO_EDGES的多物种非singleton节点也标UNRESOLVED，而非
假定支持同一OG。上述变化纳入新版本；SINGLETON/ONE_SPECIES/SIZE_STOP只是
层级流程的政策完成，均不宣称独立的生物学正确性。

### C. 搜索失败：已测试候选未接受，不等于节点不可分

搜索结果独立记录：

| 状态 | 含义 | 节点状态 |
| --- | --- | --- |
| ACCEPTED | 至少一个候选同时满足结构和政策门槛 | SPLIT |
| REJECTED_ALL_TESTED | 按预定有限日程评估完毕，无有效候选 | UNRESOLVED |
| EVALUATION_BUDGET_EXHAUSTED | 达到候选/调用上限，尚有预定点未评估 | UNRESOLVED |
| COMPUTATION_ERROR | Leiden/输入/数值错误；不参与门槛放宽 | 组件FAILED |
| POLICY_STOP | singleton、one_species、显式size stop | TERMINAL |
| DEPTH_LIMIT / NO_EDGE_SUPPORT | 多物种节点的限制/缺乏支持 | UNRESOLVED |

all_tested_no_split作为额外诊断boolean，不升级为“整个区间不可分”。
节点UNRESOLVED是有蛋白、有完整诊断的未完成分裂，不等于缺失成员或进程异常。

## 3. 有预算的回退：补充搜索，不降低可信度门槛

拟议默认日程参数（候选数均为**每节点**，是算法配置，不是Slurm资源）：

```yaml
hierarchy:
  admission_policy: nonempty_children_v1
  recursion_stop_size: 1
  resolution_strategy: bounded_adaptive_v2
  gamma_min: 0.01
  gamma_max: 10.0
  gamma_growth: 2.0
  max_candidate_evaluations: 24
  max_coarse_candidates: 10
  rescue_grid_points: 8
  local_grid_points: 5
  stability_mode: robust
  stability_threshold: 0.9
```

gamma范围/稳定性/最大子群与tiny比例沿用本轮原设置；默认数字是设计起点，
不是已经证明改善准确度的最优值。每节点最多24个unique gamma：robust最多72次
Leiden，publication最多 `24 * publication_seeds`；底层每次Leiden也需显式固定
迭代预算，首轮沿用现有版本实际默认值并在manifest记录，不顺带改为V1的10。
候选预算与seed预算校验必须一致。无嵌套无限retry，执行错误只由组件scheduler管理重试。

确定性的阶段顺序：

1. **COARSE**：从gamma_min以gamma_growth递增，范围内至多10个点。
   首个有效coarse点出现后进入REFINE，不为了耗尽预算继续找更细partition。
2. **ENDPOINT**：coarse全失败时，若gamma_max尚未测试，精确测试它。
   预留至少一个evaluation名额；不让10.24越界导致10永远遗漏。
3. **RESCUE_GRID**：仍无有效候选，评估8个固定log-space内部点：
   `gamma_min * (gamma_max/gamma_min)^(i/9)`，i=1..8。
   按升序、去重，适用于非单调/窄可行区域；范围相等则无内部点。
4. **REFINE**：若找到有效候选，在它与最近已测试的较低gamma之间评估
   local_grid_points=5的log grid（端点已测则复用），仅在剩余预算内评估新点。
5. 在实际已测有效候选中取最低gamma。tie按固定gamma key，不依赖worker完成顺序。
6. 仍无有效候选则UNRESOLVED，保留全部violations和已执行日程，不输出“最佳拒绝”划分。

10 coarse + 1 endpoint + 8 rescue + 至多5 refine ≤24；去重不会消耗候选或seed额度。
相近gamma复用规则固定为当前rel_tol=1e-12、abs_tol=0并记录canonical key；
所有点约束在闭区间，排序确定。RESCUE阶段只完成有限预定点；发现多个可行点后
REFINE最低已测有效点，不用把可行性当成单调条件二分到任意精度。
本策略保证预算与端点覆盖，不保证找到所有狭窄可行区间。

禁止自动调低stability_threshold、切换fast、放宽最大子群、合并singleton、随机换seed
直到成功，或复制V1无严格上界的SpecificBoard搜索。若另设可靠性放宽策略，
应另立明确版本与科学验收，不隐藏在本fallback里。

### 总体执行边界

24点只保证每节点调用有界，不提供全组件/全运行耗时上界。
有限蛋白、合法真分裂、max_depth提供有限节点数，但仍可能昂贵。
正式运行另设显式component Leiden-call cap（正整数或null），用于成本限制，
默认null以免在没有估算时猜限额；命中时当前待处理节点标UNRESOLVED并停止新增调用。
component cap、耗尽节点和计数进入identity；已排队节点也应保留完整未解析成员。
Slurm时间上限由脚本显式设置；TIMEOUT由checkpoint视为执行失败，不冒充搜索结论。

## 4. UNRESOLVED 的存储、发布与恢复

每个蛋白仍恰有一个**结构叶**membership，包括TERMINAL和UNRESOLVED叶。
members.parquet的terminal_cluster_id保留既有列名，但新schema明确它指结构叶，
不能只凭该列名推断该叶已完成科学处理。
新增 node 字段 search_status、termination_kind、failure_codes、selection_phase；
UNRESOLVED叶child_count=0，resolution/quality不得伪装为已接受候选值。

分离两种manifest状态：

- artifacts_written/structural_validation_passed：结构、哈希、成员分区有效；
- hierarchy_status：RESOLVED / UNRESOLVED / FAILED。

UNRESOLVED组件保留diagnostic artifacts和完整成员，但scheduler不得计为科学DONE。
run级出现UNRESOLVED时，默认pipeline返回非零并报告节点/蛋白数。
可以显式单独导出标记清楚的hierarchy诊断；常规最终OG/score链以RESOLVED为输入门槛。

**OG的纯算法、V1事件、成员排序、RESIDUAL_NONE规则及scorer保持不变。**
新增上游hierarchy verified-input gate由stage/CLI边界复用；覆盖完整pipeline、
独立orthogroups命令及默认export入口，避免只拦pipeline而直调stage仍发布巨组。
这是新hierarchy输入验收，不是在OG引擎里过滤巨组；任何显式诊断运行均不得生成
可被常规export当作DONE复用的OG manifest。

同参数resume识别已完成但UNRESOLVED的诊断，不自动重复科学搜索；
需要用户显式retry或新policy/config后重建。错误FAILED的重试仍由scheduler有界管理。
worker与parent唯一checkpoint写入职责保持原样。

## 5. 配置与缓存迁移

- 新hierarchy算法版本、schema版本、scheduler manifest版本必须同时升级；
  serial/parallel共享同一搜索实现、输入门槛和预算计数。
- 旧 `min_family_size` / min_child_size **不静默解释为 recursion_stop_size**。
  新策略下同时提供旧键和新键视为配置错误；旧run.yaml只能显式选择legacy策略
  保持旧科学行为，或由用户明确迁移为上述新配置并记录迁移diff。
- 旧hierarchy快照仍可用于历史replay，但不被新策略resume当作有效缓存。
  若程序支持legacy策略，必须单列旧版本及门槛；常规新策略使用新identity。
- identity绑定admission政策、停止尺寸、搜索阶段参数、所有预算、seed列表、
  method/weights、Leiden迭代预算与环境版本。参数相同才复用，包括unresolved状态。
- 候选存储支持violations列表与单独结构/政策valid，保持raw metrics；
  GAMMA_LIMIT不再代理CHILD_TOO_SMALL/DOMINANT_CHILD/EXCESS_TINY_FRAGMENTS。
- 修订节点模型后，annotate与导出适配schema；更改输入验收/缓存版本，
  不更改OG决策或metric数学公式。

## 6. 分阶段实施与验收（待设计确认）

H0：冻结新语义/config/状态与ADR；保留旧政策characterization fixtures。

H1：非空singleton准入、递归stop与structural leaf验证；serial/parallel一致；
覆盖1+49+50、2/3蛋白、ONE_SPECIES、depth/no-edge失败，不丢成员。

H2：endpoint/rescue/refine日程与预算；证明每unique gamma只一次、多seed调用上限，
窄可行区间、所有拒绝、数值错误、预算耗尽分别可观察。

H3：UNRESOLVED落盘、manifest-last、checkpoint/resume、上游verified-input gate；
覆盖直接CLI绕过pipeline、stale OG/export缓存与部分失败恢复。

H4：固定P6修复mean SSN，先只在component0上通过Slurm做有来源的局部验证。
沿用seed/metric，记录root候选与真实后代，不把root出现SPLIT当作准确度结论。

H5：固定新hierarchy再执行实际V1 OG parity和两策略官方Orthobench。
报告P/R/F1/coverage、超大OG污染、unresolved覆盖率；同配置对照后才作科学判断。

设计已确认；历史数据快照保持原样，本次实现不自动提交H4/H5计算。


## 7. 本地实施记录（2026-10-02）

- 默认配置和公共 dataclass 均启用 `nonempty_children_v1 / bounded_adaptive_v2`，robust 三 seed。
- 历史回放需明确 `admission_policy: legacy_strict` 与 `resolution_strategy: adaptive` 或 `log_grid`。
  `min_family_size` 仅属 legacy；新配置中该键须为 null 或省略，停止尺寸用 `recursion_stop_size`。
- 有限 Leiden 迭代数明确为2，与已安装库的默认值一致，并计入配置身份。
- 有限组件调用上限使用确定性 DFS；即使配置 subtree_workers >1，也不独立分发预算。
  默认 null 保持子树并行，串行与并行节点/候选/成员映射等价。
- `UNRESOLVED` 保存结构叶成员及完整候选诊断；不是 DONE，不按计算失败策略自动重试。
  同配置 resume 复用诊断并返回非零；`hierarchy-all --retry-unresolved` 显式重跑。
- 新节点和候选列、schema/version、SQLite CHECK 迁移、manifest-last 和上游检查已接线。
  OG stage 与最终 consumer 共用检查；纯 OG 引擎、V1 reference、官方 scorer 在本次实现中未修改。
- hierarchy manifest 保存实际阶段参数（包含 CLI overrides）；下游使用该身份与输入/输出校验和，
  不把 prepare 时的 run.yaml 默认值当作阶段的实际参数。OG cache 绑定 hierarchy manifest。
- `export --strategy terminal` 仍是独立诊断，表格保留 unresolved 状态，manifest 标明
  `HIERARCHY_DIAGNOSTIC` 与 unresolved_nodes；不是正式 OG 完成标志。
- H4/H5 未运行；小数据契约通过不代表 component0 科学问题或 Orthobench 准确度已经解决。


### 本地验收

- 完整测试：301 passed，18 个已有科学库/历史样本警告。
- Ruff、mypy（77 source files）、git diff --check 通过。
- 新增16个本地契约测试：singleton真实子节点、独立SIZE_STOP、DEPTH_LIMIT、
  端点与非单调rescue、预算、稳定性不放宽、错误membership、配置显式迁移、
  serial/parallel一致、UNRESOLVED resume/显式retry、数据库迁移、直接OG/export gate、
  实际阶段override身份以及run.yaml变化失效。
