# H4/H5：固定 repaired-mean SSN 的 hierarchy 回归

最新状态：ADR 0003 的重跑job `1410777` 已 COMPLETED / 0:0，H4完整验收和H5 parity通过；
正式OG官方F1为40.0192%，仍低于历史完整V1的57.5075%。详见 `h4-h5-iterations10-regression.md`。

下文保留job `1410775` 的2轮历史事实：FAILED / ExitCode 2:0；H4根准入修复得到真实数据支持，
但完整resolved验收未通过，H5评分未启动。

## 固定条件

来源：P6 job `1410751` / RUN_ID `20261001T123223Z_f69c14cdfe9b_c567bab0_49330`。
仅复制其 repaired-mean 的 input、edges、components，不重跑 hits/SSN。
新配置显式迁移 ADR 0002 政策与预算，保留原 seed=42、robust 三 seed、
stability_threshold=0.9、method=rber、max_depth=20、subtree_workers=1。
不改 OG 纯算法、参考 V1 或官方 scorer。

## 作业边界

`slurm/h4_h5_hierarchy_regression.sh` 单作业、无数组，顺序执行：

1. 校验 P6 completion、固定蛋白元数据、拷贝输入后的逐文件校验和；
   保存原配置、新配置与显式迁移记录、benchmark 校验和、环境版本。
2. H4：只构建 component 0 全层级，保存根候选、实际后代、停止原因、
   unresolved 节点/蛋白数、叶大小分布、每节点预算检查、耗时和资源日志。
3. H4 resolved 后才执行 H5 其余组件；同配置复用已经验证的 component 0。
4. 若 H4/H5 hierarchy 仍有 unresolved，则记录诊断、非零退出，停止正式 OG/score。
   不强行让失败候选成为 terminal，不对缺失成员补 singleton。
5. 完整 hierarchy resolved 后冻结其输入，执行独立实际 V1 OG parity、
   正式 OG Parquet/TSV 与 terminal 诊断导出。
6. 用未修改的官方 scorer，对比 previous_repaired_mean / new_bounded_hierarchy
   两种提取策略，并保留独立 historical_full_v1 对照（不是 hierarchy-only 控制）。
   输出 P/R/F1、coverage、RefOG diagnostics 的超大组污染案例及同策略变化。
7. 再次校验固定 SSN、输入及 benchmark 未改变。

资源沿用既有 P6 包络：cu / mselab / normal，56 CPU、250G、72小时。
探索运行允许 dirty，但由标准 wrapper 保存 commit、diff、所有非忽略源码快照与内容哈希；
不把该作业称为仅运行 HEAD commit。大量结果保留远端。

## 验证

本地新增摘要/配置迁移测试与 hierarchy bounded tests：14 passed。
远端经 `dev/remote-test` 同样14 passed；正式计算走 `dev/slurm-submit`。
Ruff 和 bash 语法检查通过。作业提交/完成均不等于科学准确度改善。


## 提交记录

- Job ID：`1410775`。
- RUN_ID：`20261002T034016Z_f69c14cdfe9b_e22b9a32_98902`。
- 基础 commit：`f69c14cdfe9bee9ed1aa251e89a511382b23baca`，dirty=1；不是纯 HEAD 运行。
- source SHA256：`e22b9a32c1415fd81ae32bf5c0c876e86462352de209448dd388837de4305128`。
- 远端目录：`/home/mselab/licj/projects/running/ogprofiler-runs/20261002T034016Z_f69c14cdfe9b_e22b9a32_98902`。
- 本地 provenance：`.provenance/slurm/20261002T034016Z_f69c14cdfe9b_e22b9a32_98902/`。

本节在提交后补充；作业执行的是提交时不可变快照，之后本地文档更新不影响该作业。


## 实际结果（2026-10-02）

Slurm记录：FAILED，ExitCode `2:0`；总耗时 `00:54:56`，
2026-10-02 11:38:53–12:33:49（Asia/Shanghai）。
H4 hierarchy本身耗时约4分07秒（engine 237.89秒），峰值RSS约815MiB；
不是OOM、TIMEOUT或计算异常退出。CLI因持久化后仍存在UNRESOLVED返回2。
`h4-h5-completion.json`：h4_exit=2，h5_status=BLOCKED_BY_H4，scoring_started=false。

### 已验证的修复：根的singleton准入和稳定性门槛

| 指标 | P6 repaired-mean旧hierarchy | ADR 0002新hierarchy |
| --- | ---: | ---: |
| 根成员 / 物种 | 69,642 / 12 | 69,642 / 12 |
| 根状态 | 普通TERMINAL / UNSTABLE | SPLIT / ACCEPTED |
| 根真实子群数 | 0 | 444 |
| 选择gamma | 无 | 0.11313708498984759 |
| 最小子群 | 旧规则拒绝singleton | 1，实际合法子节点 |
| 选中稳定性 | 无 | 0.9022658582695411 ≥0.9 |
| 根候选 / 调用 | 10 / 30 | 8 / 24 |

新根gamma=0.01指标与P6一致（148群、stability=0.8347310184583124）。
gamma=0.16同样496群、stability=0.9041826997009736；新规则不再因
min_child_size=1否决该合法候选，随后local grid选择更低已测有效gamma。
根此次在coarse/local阶段成功；根的成功体现singleton准入与局部搜索，endpoint/rescue在根上尚未触发。

### 完整后代与新的碎片化风险

- 总节点81,443：SPLIT 28,176、政策TERMINAL 53,266、UNRESOLVED 1。
- 结构叶53,267，成员69,642全部唯一保留；68824蛋白进入政策停止叶，818留在unresolved叶。
- SINGLETON叶49,586，占component0成员 **71.2013%**；ONE_SPECIES叶3,680。
- 其中5,552次分裂的父节点仅2蛋白，1,975次仅3蛋白。实际最大depth为20。
- 288,402个候选、865,206次Leiden调用，严格等于3×候选数。
  每节点最多22候选（配置上限24），所有选中候选通过结构/政策/stability检查；
  所有选中稳定性的最小值0.9001357925064575。

这表明修复解决了“因一个singleton拒绝整根”的机制，但也暴露出
`recursion_stop_size=1` + 逐层寻找任何稳定分裂的科学后果：多物种小组继续拆到singleton。
稳定性只说明多seed一致，不证明分裂具有正交生物学意义。
这是**过度碎片化风险信号**，不是已获得的Orthobench准确率结论；
后续若修改停止尺寸或分裂目的，需显式修订设计，保持OG/scorer不变。

### 唯一unresolved后代的具体原因

component0 / cluster_id **67915**，depth=2，818蛋白、9物种；
路径 `0 → 64983 → 67915`。状态为REJECTED_ALL_TESTED，而非预算耗尽。
完成10个coarse、1个精确gamma_max=10端点、8个预定logspace rescue，共19候选/57调用；
全部拒绝，所以没有有效点触发local refinement。

- gamma=.01：1群，NO_SPLIT + MAX_CHILD_FRACTION；stability=1。
- gamma=.02：2群，最大子群占比0.990220 >0.95；stability=1。
- rescue gamma≈.021544：同样超最大子群比例，且stability=.333333。
- 17个候选包含UNSTABLE；其中16个满足大小比例，仅因稳定性低于.9被拒绝，
  另1个是上述rescue的同时拒绝。
- gamma=1.28是满足比例约束时最接近门槛者：49群，stability=.8973204470525067。
- 精确上界gamma=10：120群，stability=.8748040717220942，仍未过门槛。

818蛋白占component0 **1.1746%**。没有强收该划分、降低阈值、把失败叶变普通terminal，
也没有为使评分成功而补成员或发布OG。
REJECTED_ALL_TESTED只说明有限已测点均失败，不证明所有gamma不可行。

### 本地证据校验与H5状态

约5MB的component0 nodes/members/candidates/metrics与manifest、摘要和日志已取回。
`local-validation-summary.json`确认：输出校验和匹配、成员唯一覆盖、叶计数一致、
每节点预算和3-seed调用计数有效、所有选中候选稳定性门槛有效。
文件位于提交记录对应的本地provenance `results/`；大型SSN和输入继续保留远端。

**H5未执行**：其余组件层级、完整V1 OG parity、两策略官方Orthobench均未启动。
没有新的P/R/F1结果，也没有科学准确度恢复结论。默认输入gate如设计阻断发布；
本次未自动重提作业、变更参数或扩大资源。


补充轻量ID检查：取回818蛋白的原始ID（9物种）及该作业冻结的RefOG文本，
与全部RefOG（含low-certainty文本）逐ID交集为0。
这仅是局部基准成员定位，不是P/R/F1评分，也不解除完整hierarchy的输入门槛。


## 后续审计澄清（2026-10-02）

见 `recursion-stop-and-818-node-audit.md`。V1同样将多物种2蛋白直接拆为singleton；
结构叶singleton比例不是最终OG的碎片化比例，先前的风险信号尚未得到OG评分验证。
818冻结图上固定阈值.9的局部诊断显示，n_iterations=10时公共search获准一个
44群/ARI .913175的候选。这是优化预算反事实诊断，生产与job1410775事实保持原样；
尚未重新执行H4/H5或改动生产参数。

### ADR 0003 的后续配置修订

新策略与回归迁移工具已显式切换 `leiden_iterations: 10`；其余停止尺寸、稳定性
阈值、seed和固定mean SSN不变。hierarchy/scheduler版本升级至v3，旧2轮缓存失效。
上述表格仍只描述job1410775的2轮历史结果；新预算的818根搜索本地复核通过，
但尚未提交新的完整H4/H5。迁移记录会同时保留来源配置与新配置，不覆盖历史结果。

## ADR 0003 重跑的执行门槛

脚本接收 P6 来源和前次 H4 RUN_ROOT 两个显式参数；逐项验证新配置与job1410775
仅 `hierarchy.leiden_iterations: 2 -> 10` 不同，并核对双方固定输入哈希。
主H4保持serial；成功后在独立目录对component0完整后代执行两worker回放，
只改变执行并行度，保持原release尺寸。nodes/members/candidates逐表精确比较，
验收唯一完整成员覆盖、叶尺寸、所有选中候选稳定性和预算、robust调用计数与manifest。
`h4-acceptance.json`全部通过才进入H5。两次运行分别保留time和日志；对照运行
增加计算量，不改变既有Slurm资源申请。历史结果与原源快照保持不变。

## ADR0004 生产 opt-in 的后续回归准备（2026-10-03）

见 `soft-policy-production-regression.md`。脚本可选第3参数soft_binary_24_v2，
须以冻结10轮kway H4作为对照，只允许topology政策不同；旧2→10默认路径仍保留。
algorithm v4/schema3，重建hierarchy并验收fallback证据。当前只有五个小组件本地公共
串行/双进程完整回放通过，尚未提交新H4/H5，未产生新的官方P/R/F1结论。
