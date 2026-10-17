# 第2项：132蛋白组件搜索下界对照（实验预定义）

生产下界不变，ADR0004仍为Proposed。本实验显式测试不同下界，不等于接受下界迁移。

- 固定1410777 mean SSN的component84（132蛋白），完整root graph与冻结metadata。
- 固定覆盖策略A `adr0004-coverage24-fair-v1`。
- 预定义三下界 `.01`（控制）、`.001`、`.0001`，均保持gamma_max=10。
  使用整十倍级差，不将诊断已知可行gamma作为固定候选；所有节点使用同一实验下界。
- 其余参数：停止尺寸1、ARI .9、比例.95、三seed、10轮、OG/scorer保持不变。
  预算仍24 unique gamma/72调用，每个不同下界是独立配置/实验，不把三次合并称为一次24点搜索。
- 先根节点三配置各重复两次，检查原门槛与调用预算；再三配置分别从根完整递归各重复两次。
  baseline始终使用冻结原配置，实验下界仅传给hard-binary分支，避免下界变动污染控制。
- 验收成员唯一覆盖、完整后代失败节点与violations、候选预算、重复性以及OG门禁。
  若仍UNRESOLVED，保留完整成员并阻断OG，不以低gamma根成功宣称组件成功。
- 范围扩大改变搜索日程及潜在整树拓扑，需要新配置身份；不与旧失败缓存混用。
  本轮不处理稳定性/比例门槛，不引入soft k-way回退，不接线生产。

## 根节点结果（每项独立重复两次）

|gamma_min|根结果|选中gamma|ARI|最大子群比例|候选/调用数|
|---|---|---:|---:|---:|---:|
|.01|REJECTED_ALL_TESTED|—|—|—|24 / 72|
|.001|ACCEPTED|.001|1|99/132=.75|1 / 3|
|.0001|ACCEPTED|.0005|1|99/132=.75|7 / 21|

三配置均重复一致、所有点在自身显式范围内；低下界两项的根分区相同。
.001根在下界即合格，refine没有更低的范围内区间，重复gamma不收费。
.0001不是将已知.0005插入点集，而是coarse/refine按一般规则实际测得。
本局部控制支持“原下界限制了根binary搜索”，不支持“任意节点都应使用更低下界”。

## 完整递归验收

|gamma_min|节点数|候选/调用数|剩余失败节点（cluster_id:尺寸）|耗时（本机）|
|---|---:|---:|---|---:|
|.01|1|24 / 72|0:132|1.12秒|
|.001|19|178 / 534|2:91，18:33|1.48秒|
|.0001|19|214 / 642|2:91，18:33|1.97秒|

- 两个低下界的新树terminal membership相同，但各节点的选中gamma与候选日程不同。
  更低下界没有减少本例unresolved，不能从根成功推断整树改善。
- 三树均132蛋白唯一完整覆盖、validate_hierarchy通过；每节点<=24点/72调用，
  实际调用数严格为unique候选数×3；每项完整递归重复nodes/members/candidates一致。
- baseline每次仍使用原冻结配置，节点/成员/OG分区均与1410777相符。
- 三实验树均BLOCKED_UNRESOLVED；未调用实验OG/scorer，不给出RefOG recall/F1。
  Slurm尚未提交，生产串并行一致性尚未验收。

## 91/33蛋白后代复查

.001树的原24点trace：91蛋白节点10个binary、33蛋白节点9个binary，均ARI=1，
但最大子群比例分别89/91=.978021978、32/33=.969696970，全部MAX_CHILD_FRACTION。
.0001树同两节点也全部binary比例失败（分别8与3个binary）。

保持这些成员与同一SSN，以.001配置进行**独立有限诊断**：各129个全范围log点＋
65个首次原始count跨越区间的等距点，去重后各194点，总1164调用。
不把这些额外点计入原24点搜索，不借诊断绕过预算。

- 91蛋白：17个诊断binary均比例失败，无合格binary采样点；
  另有75个原门槛合格k-way点，最低测得gamma≈.04216965、4-way。
- 33蛋白：26个诊断binary均比例失败，无合格binary采样点；
  另有61个原门槛合格k-way点，最低测得gamma≈.07498942、3-way。

这些k-way被hard-binary政策排除，不是OG/scorer失败；它们只是将来显式soft政策的候选证据，
本轮未发布或递归执行这些候选。有限诊断未找到合格binary，不证明整个区间无解。

### 同时发现的覆盖限制

A在第一次raw binary后跳过global endpoint，若没有已测更高gamma，只能探低侧。
.001的91/33节点原日程最高gamma分别.016/.064，没有测试到上述高gamma合格k-way。
更极端地，下界首点即raw binary但被拒绝时，A可能只有一个已测点、无相邻区间而结束。
这不是增加下界范围本身能修复的问题；本轮为保持受控对照未静默改变A日程。
后续应显式版本化“缺失上边界补探”规则，并与soft fallback所需覆盖共同设计。

## 本轮决策

1. **不迁移生产下界**：低gamma根成功可重复，但两种低下界均未通过完整递归验收。
2. 若继续范围实验，.001已经解决本根范围限制，.0001在此例没有额外resolved收益，成本更高；
   这只是本例判断，不作全数据默认参数结论。
3. 下一步先明确补上边界的共享预算及失败状态，再设计显式标记的soft binary→k-way回退；
   比例.95、ARI .9、停止尺寸1及OG/scorer继续保持不变。
4. 后续验收需包括新4蛋白稳定性反例、新91/33比例反例与原五个完整组件，
   不只复查132蛋白根；仍unresolved时保持门禁，不进入正式评分。

## 证据与验证

1410777本地results：`binary-lower-bound-comparison.json`；派生失败输入
`binary-lower-bound-unresolved-input.json`记录母报告hash，诊断
`binary-lower-bound-failure-audit.json`记录派生输入hash、脚本hash及完整候选。
报告记录完整实验配置、生产与实验源hash；SSN输入逐项核对冻结edge manifest。

新增下界测试2个，连同既有覆盖/诊断/hierarchy/bounded-stage测试，
**56个针对性测试通过**（4.29秒）；Ruff、git diff --check通过。
生产src保持不变，未提交Git或Slurm任务。

```sh
R=.provenance/slurm/20261002T060705Z_cb0c0ba98049_1cb4fcd7_14702/results
PYTHONPATH=src:. .venv/bin/python benchmarks/og_extraction/binary_lower_bound_comparison.py \
  --result-root "$R" --output /tmp/binary-lower-bound-comparison.json
```
