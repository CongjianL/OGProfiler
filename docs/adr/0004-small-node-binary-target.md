# ADR 0004：小节点binary优先与有限soft回退

> 2026-10-06：本文默认配置限制已被 [ADR0007](0007-soft42-production-default.md) 后续决定取代；下文保留历史实验范围与结果。
- Status：**Accepted for explicit opt-in（2026-10-03）；默认仍为 `kway_v1`，全数据验收待 H4/H5。**
- Date：2026-10-02。
- Amends（仅 opt-in）：ADR0001「Stable k-way splits are preserved」在小多物种节点上的适用范围，
  及ADR0002搜索日程；不修改ADR0003每次10轮、停止尺寸或稳定性定义。
- Evidence：`../binary-target-comparison.md`、`../refog-first-separation-audit.md`。

## 背景

V1对3..9999蛋白多物种节点以exactly2为目标；V2接受任意稳定非空k-way。
多路分裂的III-3事件会失去V1 I事件选择资格。固定node2694:3的11个单拷贝成员
在现有点集中选择11-way/III-3；gamma .8125的binary同时通过原门槛且产生I。
但其他binary仍为III-1，另有binary违反max_child_fraction=.95的反例。
所以binary不是准确度保证，仅增加cardinality条件也不足以保证有限搜索找到可行点。

## 拟议政策边界

1. 小多物种节点 `3 <= n < 10000` 拟引入独立可配置的binary目标政策。
   此10000边界来自V1源代码，不改变 `recursion_stop_size: 1`。
   >=10000暂保留现行k-way；不同时迁移V1的20..50大节点目标。
2. binary目标是代表性membership的child_count==2，且**全部现有结构与政策门槛仍须通过**。
   representative仍按现有多seed支持度/quality选取；若要求每seed都恰好2群，
   那是额外可靠性政策，须单独明确，不隐藏在本cardinality条件中。
3. 两蛋白节点继续现行搜索/停止语义；本稿不引入V1直接制造两个singleton的特判。
   ONE_SPECIES、singleton、max_depth、输入验收保持现状。
4. 固定method/mean权重、seed、robust mean ARI>=.9、max_child_fraction=.95、
   tiny比例与quality配置、每调用10轮。OG/scorer及事件标注公式保持原样。

## 搜索与预算必须共同修订

现有coarse/endpoint/rescue/local点对2694:3没命中可行binary；binary目标需显式
增加count-target探测，而不是继续接受首个有效k-way。

### 固定实验日程 `adr0004-hard-binary-24-v1`

本日程已确定，用于独立小组件整树实验；ADR仍为Proposed，生产默认保持原政策。
每节点共享24个unique gamma，robust固定3 seed最多72次Leiden调用；每调用10轮。
阶段上界：**COARSE 10 + ENDPOINT 1 + TARGET_PROBE 5 + RESCUE 5 + REFINE 3 = 24**。
阶段不追加额度、不补满未使用额度；所有阶段复用已测点（rel_tol=1e-12，abs_tol=0），
复用不消耗unique gamma或Leiden调用。gamma闭区间仍为[.01,10]。

1. **COARSE**：从.01开始每次乘2，最多10点，超过10则停；首个合格binary后立即进入REFINE。
2. **ENDPOINT**：coarse未找到时测10；合格则进入REFINE。
3. **TARGET_PROBE**：仍无合格binary时最多5次新点。将已测点按gamma排序，
   只考虑相邻区间，满足 `left.child_count < 2 <= right.child_count` 或任一端原始count==2。
   按 `(left.gamma, right.gamma)` 字典序取最小区间，测算术中点；插入结果后重新建立相邻区间。
   无区间或中点与已测点数值等价则停止此阶段。首个合格binary进入REFINE。
   原始count==2但政策失败时两侧区间都保持可选，低gamma优先；不把失败当作count<2。
   不假设count或稳定性单调；count>2两端内部可能存在binary但此启发式可能漏测。
4. **RESCUE**：仍无合格binary时测5个对数内点 `.01 * 1000**(i/6)`，i=1..5。
   完成此有限列表，随后进入REFINE（若已有合格点）；目标探测与rescue共享原24点，而非各24点。
5. **REFINE**：取已测合格binary的最低gamma为upper；最近的已测严格更小gamma为lower，
   无lower则使用gamma_min。固定该区间，测1/4、1/2、3/4算术内点，共最多3个新点；
   端点已测，不额外收取2点。最终选择**实际测得且合格的最低gamma**。

cardinality按现有代表性membership判断，新增 `TARGET_CHILD_COUNT`，
保留原始UNSTABLE/MAX_CHILD_FRACTION等violations；不增设每seed都恰好2群的门槛。
该日程是有限启发式，不证明未采样区间无解，不保证跨节点递归成功。
实验配置身份同时记录protocol、完整冻结HierarchyConfig和脚本hash；生产配置/缓存接线尚待实施。

## 失败与回退

只有binary目标、原门槛和结构同时通过，才可在hard binary政策下发布分裂。
预算内未找到合格binary保留UNRESOLVED及完整诊断，不冒充政策停止叶。

- `REJECTED_ALL_TESTED`：既定有限日程完成/按上述规则停止，全部已测点被拒绝。
  使用21点后失败并不是“还差3点自动补搜”；3点是仅成功后使用的refinement配额。
- `EVALUATION_BUDGET_EXHAUSTED`：另一个日程请求的新点因共享上限被阻断；
  仅恰好用满上限不自动产生此状态。若已有合格点，refinement截断仍可ACCEPTED，生产接线时须显式记录refinement截断；本次固定24点实验不触发该分支。
- `ACCEPTED`：实际合格binary发布；n=2及>=10000沿用原k-way搜索。
- max_depth/no-edge及组件总调用预算沿用现有UNRESOLVED语义；执行异常直接失败，不伪装候选拒绝。
- UNRESOLVED是结构叶但不是生物学政策停止叶，保留全部成员；组件OG发布与正式评分被门禁阻断。
  原始候选诊断仍可输出。失败节点重新运行须有新protocol/config身份，不复用旧失败缓存当作成功。
不为使binary成功而降低ARI、放宽最大子群、强拆最大群或合并多路子群。
V1搜索未得exactly2后回到普通terminal的行为不移植为新默认。

如果选择“先binary，预算失败后允许原k-way”的soft目标，需要另定义可配置政策、
selection/fallback标签与对照，明确它不保证小节点binary；不得静默回退后仍称binary完成。
当前没有自动批准这一soft回退，也没有批准将失败节点保留为普通terminal。

## 接线与验收前提

新目标、适用尺寸、阶段配额、失败/回退语义须绑定配置身份、候选trace及缓存版本；
serial/parallel必须等价。是否新增Parquet列及schema升级由具体设计决定。

先做固定小组件整体递归对照，不把当前冻结节点的一个I候选当作新整树效果。
覆盖正例2694:3、III-1反例470:0/84:5、比例拒绝844:2/500:0，以及2蛋白与10000边界。
验收resolved比例、完整成员覆盖、调用预算、同hierarchy OG parity与RefOG关系保留。
通过小组件后再显式提交不可变Slurm H4/H5验证；不得跳过UNRESOLVED输入门槛。
本次草案不自动启动全组件参数扫面、生产迁移或Slurm回归。

## 2026-10-02 修订：实验日程v2（替代上述v1作为下一轮实验日程）

Status仍为**Proposed；生产未启用**。v1规范与结果保留作为冻结对照。
依据：[失败覆盖与准入冲突审计](../binary-failure-coverage-and-schedule-v2.md)。
修订配置身份为 `adr0004-hard-binary-24-v2`，不得复用v1实验结果缓存。

### 分配与顺序

**COARSE 8 → ENDPOINT 1 → RESCUE 3 → TARGET_PROBE 9 → REFINE 3 = 24**。
前3阶段最多12点，先用全局rescue补充count证据，再用9点目标探测。
若仍失败，显式将3个REFINE名额借给目标探测（trace为`BORROWED_REFINE`）；
目标探测合计最多12点。借用后才找到合格候选时，只用剩余共享名额细化，最多3点。
不追加第二套预算；不补满未使用阶段额度；重复点不收费。

- coarse从gamma_min乘growth_factor，最多8点，超上界则停。首个合格binary跳到refine。
- 无合格点则测gamma_max；仍无合格点则测3个log内点
  `gamma_min*(gamma_max/gamma_min)**(i/4)`，i=1..3。rescue有限列表完成后选最低合格点细化。
- target每步重新排序已测点，只考察数值上仍可插入新点的相邻区间：
  1. 已存在被拒绝binary时，找下边界（左端count!=2、右端count==2）与
     上边界（左端count==2、右端count!=2）。两侧都存在时**下/上轮转，下侧先行**；
     只有一侧则探该侧。每侧选择log宽度最大区间，tie按lower/upper gamma升序。
  2. 没有上述边界时，选count从<2跨到>=2的区间。
  3. 再无则选两端count==2的区间，检查可能非单调的稳定性/准入变化。
  4. 再无则探所有其他区间，优先log宽度最大，使用几何中点，避免“下界已k>2即无探测”。
  前3类用算术中点，tie规则一致。不把UNSTABLE或比例失败视为count==1。
- refine沿用v1固定lower/upper的1/4、1/2、3/4点；借用后可截断，报告显式列出
  `refinement_truncated_nodes`。始终选择实际合格点中的最低gamma。

上/下轮转是有限覆盖政策，不保证稳定窗口被命中。
仅优先上边界会漏掉2694根的下侧合格binary；仅优先下边界又会饿死9蛋白后代上侧。
本修订避免单侧饥饿，不为单个已知gamma插入特制固定点。

### 失败语义保持与澄清

完成包括借用在内的有限日程，仍无合格点：`REJECTED_ALL_TESTED`，即使用满24点也不改名。
若显式更小共享cap提前阻断既定新点：`EVALUATION_BUDGET_EXHAUSTED`。
已有合格点而refinement被共享cap截断：仍ACCEPTED，同时记录截断，不将有效分裂丢弃。
结构成员保留、UNRESOLVED/OG门禁和异常语义沿用前文。

**不改变**gamma_min=.01、gamma_max=10、停止尺寸1、ARI .9、max_child_fraction .95、
优化10轮、seed、代表性membership、n=2/10000边界，以及OG/scorer。
下界外诊断不构成下界修订，比例反例不构成放宽准入批准。
v2整树仍0/5 fully resolved，保持实验状态；不接线生产或启动大规模Slurm验证。

## 第1项覆盖修订：自适应预算候选（2026-10-02）

仍为Proposed、生产未启用。证据及预定义A/B规范见
[固定24点覆盖对照](../binary-coverage-comparison.md)。

下一轮实验候选选择 `adr0004-coverage24-fair-v1`（A）：

1. coarse最多8点，首次原始count>=2即停。发现原始binary或已测相邻count<2→>=2区间后，
   跳过endpoint及global rescue，将未用额度交给target；否则保留endpoint+3个log rescue。
2. target沿用v2公平上下边界轮转及区间选择；在总24点内动态分配，最后3点可显式借用。
   无合格点时完成整个24点有限日程；早期成功用剩余额度refine最多3点。
3. 所有gamma由边界/中点规则生成，未加入节点专属点。不假设稳定性单调。
4. full-cap有限日程完成无合格点：REJECTED_ALL_TESTED；更小cap阻断：EVALUATION_BUDGET_EXHAUSTED；
   已合格而refine截断仍ACCEPTED，候选trace和报告标明截断。
5. A固定9蛋白及整树9蛋白节点均在第24点合格；但整树产生4蛋白UNRESOLVED，完整递归验收未过。
   对照B内部稳定区域优先未命中9蛋白窗口，不采用其优先政策。

仍保持停止尺寸1、ARI .9、比例.95、[.01,10]、10轮、seed、OG/scorer及2/10000边界。
新增4蛋白诊断中的binary均UNSTABLE；不将其改为普通terminal，不启用隐式k-way回退。
本节只修订实验覆盖日程，不批准生产政策迁移或声称整个组件已修复。

## 第2项：搜索下界实验（2026-10-02，未接受生产下界变更）

实验方案与证据见[132蛋白下界对照](../binary-lower-bound-comparison.md)。
以A日程独立对照gamma_min=.01/.001/.0001，gamma_max=10及其他准入参数保持不变。
实验各自拥有24点共享预算，配置身份含显式下界，baseline仍用冻结.01配置。

.001/.0001均可重复选中根binary（gamma .001/.0005、ARI=1、99/132），但完整树均
留下91/33蛋白比例拒绝节点，因此**不批准下界迁移，生产仍为.01**。
诊断中91/33的binary全比例失败，同时测得原门槛合格的k-way候选；这些候选本轮未准入。
A缺少已测上边界时可能只探低侧，后续须另明确补上边界的日程与共享预算，不能静默改写A身份。

本节不降低最大子群比例/稳定性门槛，不启用soft fallback，不绕过UNRESOLVED的OG门禁。
下一拓扑回退设计需覆盖新4蛋白稳定性反例及91/33比例反例，并做完整后代验收。

## 第3项：显式soft政策候选（2026-10-03，Proposed）

完整规范与证据：[上边界与soft回退实验](../binary-soft-fallback-experiment.md)。
本节显式修订小节点exactly2硬目标为**binary优先，有限搜索失败后可选择已测合格k-way**，
不把它称为hard-binary完成，不静默回退。

候选政策身份：`adr0004-soft24-upper-guard-v2`。在A共享24点内，对coarse第一次
原始count>=2但原门槛拒绝的分裂，按growth_factor补上侧，直到合格分裂、gamma_max，
或只剩3个目标探测名额。binary仍优先；失败时仅复用已测原门槛合格k-way最低gamma，
选择phase明确为FALLBACK_KWAY/原phase，分别记录binary/kway eligibility与原violations。
更小cap提前阻断不触发soft回退；无合格候选仍UNRESOLVED，OG门禁不变。

.01原下界五个小组件5/5完整resolved；18个回退均通过原门槛，成员271唯一覆盖。
同新hierarchy的实际冻结V1 OG parity五组件全部通过；局部RefOG真对保留236→694，
但不是官方全数据P/R/F1，RefOG012仍231→231，不能据此宣称全部低召回已修复。

仍保持停止尺寸1、ARI .9、比例.95、[.01,10]、10轮、seed、OG/scorer及2/10000尺寸边界。
生产默认未迁移；待政策确认及配置/cache/schema/fallback trace、serial/parallel验收后，
才进入不可变Slurm全组件与官方benchmark。本次小组件实验不是生产发布批准。


## 生产接线决议（2026-10-03）

用户在五个小组件实验验收后要求进入下一步：实施显式 opt-in，不迁移默认。
本节为当前有效决议；前面的 Proposed/hard-binary/A 日程均保留为设计历史。

- 配置 `hierarchy.topology_policy: soft_binary_24_v2`；默认 `kway_v1`。
- 生产日程完整对应 `adr0004-soft24-upper-guard-v2`，在公共 search 内执行，
  不依赖 benchmark monkeypatch。作用域3..9999；ONE_SPECIES、2蛋白、>=10000处理不变。
- 仅兼容 bounded_adaptive_v2/nonempty_children_v1，共享 cap 至多24。
  coarse 使用现有 max_coarse_candidates 与8的较小值；条件 rescue 使用现有
  rescue_grid_points 与3的较小值、固定四分对数点；REFINE最多3。
  默认配置复现冻结实验。降低这些阶段上限是显式搜索配置变更，须记录完整配置。
- binary 优先；有限日程完成后才复用已测、原门槛合格的 k-way 最低 gamma。
  更小 cap 阻断仍 EXHAUSTED，不触发回退；full-cap完成仍无解为REJECTED_ALL_TESTED。
  两者递归均保持 UNRESOLVED；已合格 binary 的 refine 截断仍 ACCEPTED。
- candidate 落盘独立记录 binary_eligible、kway_eligible、original_violations、
  selection_kind、refinement_truncated；node 同步 selection_kind/refinement_truncated。
  fallback 的 phase 保留 `FALLBACK_KWAY/`，manifest 统计 fallback 与截断次数。
- hierarchy algorithm v4、scheduler v4、Parquet schema3；topology 属于完整配置身份。
  原v3/schema2缓存失效，opt-in/default互相切换也失效；同配置新结果正常 resume。
- 停止尺寸1、ARI .9、最大子群比例 .95、gamma [.01,10]、10轮、seed、meanSSN、
  OG/scorer保持不变。原结构校验与 unresolved 发布门禁保持不变。

验收见 [生产 opt-in 回归](../soft-policy-production-regression.md)。五组件公共串行/双进程
节点、成员、候选、调用数完全一致，三张Parquet表完全一致，数值/拓扑复现原冻结实验，
V1 OG八项parity全部通过。此结果只批准生产opt-in接线，不代表全数据准确度恢复。

## H4新证据（2026-10-03，job1410810）

生产opt-in完整component0留下143个DEPTH_LIMIT/570蛋白，H5门禁阻断。
见 `../depth-limit-binary-policy-audit.md`。路径变深来自binary与fallback日程共同作用；
全部已选候选原门槛合格，但小组件5/5的成功不代表大组件完整resolved。
ADR0005提出H4-only深度42资源对照，尚未接受生产预算变更。

后续用户确认进入ADR0005的H4-only深度42实验；全局默认深度仍20，不自动进入H5。
