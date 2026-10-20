# 第1项：固定24点的覆盖策略对照（实验预定义）

本轮只调整搜索覆盖：gamma仍为[.01,10]、停止尺寸1、ARI .9、比例.95、10轮、三seed、
mean SSN以及OG/scorer保持不变。不处理下界或比例冲突，不启用soft k-way回退。

比较三种版本：原v2；`adr0004-coverage24-fair-v1`（A）；
`adr0004-coverage24-stability-v1`（B）。A/B在相同自适应阶段分配下仅改变候选选择，
避免把阶段分配差异误算成稳定性探测优势。

## 预定义日程（运行前固定，非节点专属gamma）

- coarse最多8点，gamma_min起每次乘growth_factor；**首次原始count>=2即停**，不要求政策valid。
- 已有原始binary或相邻count<2→>=2区间：跳过endpoint与global rescue，将未使用额度交给target。
- 否则：测gamma_max，仍无合格binary时再测3个log四等分内点。
- 共享总cap24，早期成功保留最多3个refine点；未成功时借用最后3个名额继续target。
  从target预算动态分配，而不是固定阶段空置额度；重复点不收费，不重复seed调用。
- A：沿用v2的上下边界轮转（下侧先行）、最宽log区间优先；无边界则跨count区间、binary内部、
  最宽其他区间依次选择。前3类算术中点，其他类几何中点。
- B：使用相同规则，但每第3个target名额优先探测两端count==2、至少一端UNSTABLE的内部区间。
  按端点ARI差的绝对值降序，再按log宽度降序、gamma升序，测算术中点。
  无此区间则退回A的选点规则。不将ARI变化当作单调函数，不对未测点插值准入。
- 到24点后完整日程无合格binary：REJECTED_ALL_TESTED；显式更小cap阻断后续请求：
  EVALUATION_BUDGET_EXHAUSTED；已合格但refine截断仍ACCEPTED，明确记录截断。
  失败仍UNRESOLVED、完整成员保留、阻断OG评分。

固定节点以v1失败树成员为输入：2694:7（9蛋白），加2694:0（15蛋白根）和470:25（8蛋白）控制。
三版本均实际评估，记录所有点及调用数，并重复检验；不从密集诊断读取gamma。
然后A/B分别从根重跑五个完整组件，检验完整覆盖、全后代、门禁及重复性。
这是已知9蛋白反例上的开发对照，不是独立盲测；即使成功也需后续独立组件验证。
生产默认保持不变。

## 实测结果

证据：1410777本地results目录的`binary-coverage-comparison.json`，包含来源脚本hash、
v1冻结输入报告hash、精确成员、所有候选、调用数、config及完整递归结果。
每个固定节点三版本各运行两次；两个完整递归版本各重复一遍。

### 固定图对照

|固定节点|v2|A：公平边界|B：内部稳定区域|
|---|---|---|---|
|2694:7，9蛋白|24点失败|**24点合格，gamma .9584375**|24点失败|
|2694:0，15蛋白|22点，gamma .62|14点，gamma .62|14点，gamma .62|
|470:25，8蛋白|21点，gamma .97|17点，gamma .97|17点，gamma .97|

A在9蛋白第24点找到ARI=1、最大子群8/9的binary，没有剩余refine名额；
明确标记refinement截断，而非谎称细化完成。原日程在8 coarse + endpoint + 3 rescue
用掉12点；A在有count跨越区间后跳过4个全局点，给目标探测16点。
合格gamma由一般中点算法实际生成，未加入已知诊断gamma。

B同样获得16个target名额，但内部稳定区域探测分走边界名额，最终仍漏测窄窗口。
**本反例支持自适应额度分配＋公平边界覆盖，不支持把内部稳定性探测优先规则设为默认。**
这不是所有图上A优于B的结论；已知反例存在开发选择偏倚，后续仍需独立验证。

### 从根完整递归

|组件|A节点/调用数|A剩余失败尺寸|B节点/调用数|B剩余失败尺寸|
|---|---:|---|---:|---|
|2694|23 / 546|**4**|13 / 303|9|
|470|41 / 843|26|41 / 843|26|
|84|1 / 72|132|1 / 72|132|
|844|1 / 72|33|1 / 72|33|
|500|1 / 72|45|1 / 72|45|

A新树中的9蛋白node7亦选中gamma .9584375，随后8/6/5蛋白后代继续分裂，
产生新的4蛋白node11失败。根gamma仍为.62，不能仅凭固定9蛋白成功宣称组件修复完成。
两方案均0/5整组件resolved，全部BLOCKED_UNRESOLVED，未进入OG/scorer，无实验F1结论。

所有输入271蛋白完整唯一覆盖、结构验收通过；每节点<=24点/72调用，
调用数=候选数×3；重复nodes/members/candidates一致；baseline节点/成员/OG分区
再次匹配冻结1410777。没有生产串并行一致性结论。

### 新4蛋白失败复查

将A整树的2694组件单独作为冻结诊断输入（记录母报告hash），以同一mean SSN
重建node11的4蛋白完全图（6边、连通）。额外194点/582调用作为独立有限诊断，
不混入原24点搜索预算，使用既有log129＋首次count跨越区间65点日程。
证据：`binary-coverage-fair-unresolved-input.json`、`binary-coverage-fair-new-failure.json`。

诊断原始count分布：1群89点、2群62点、3群1点、4群42点。
**62个诊断binary均ARI=1/3，全部UNSTABLE，没有合格binary采样点。**
A原24点trace中的binary也均UNSTABLE（ARI约0或1/3）。
这说明新失败已表现为实测稳定性准入冲突，而不是从本次样本中发现另一个漏测合格窗口；
不证明整个区间不存在稳定binary，也不据此降低ARI阈值或把4蛋白普通终止。

## 结论与本轮边界

采用A作为**下一轮实验候选**，不采用B的稳定区域优先规则。
第1项的固定9蛋白搜索覆盖修复有证据；**完整递归验收未通过**，生产默认保持原样。
不继续对已知9蛋白gamma调参，不静默放宽新4蛋白节点的稳定性门槛。
新4蛋白冲突与既有范围/比例冲突都要保留UNRESOLVED，后续另行明确拓扑回退政策。
本轮不执行搜索下界修改、比例门槛修改、soft k-way回退、生产接线或Slurm提交。

```sh
R=.provenance/slurm/20261002T060705Z_cb0c0ba98049_1cb4fcd7_14702/results
PYTHONPATH=src:. .venv/bin/python benchmarks/og_extraction/binary_coverage_comparison.py \
  --result-root "$R" --output /tmp/binary-coverage-comparison.json
```

## 验证状态

新增覆盖策略测试7个，连同既有诊断、hierarchy和bounded-stage测试，
**54个针对性测试通过**（3.36秒）；Ruff及git diff --check通过。
固定图/整树重复性、冻结baseline parity及派生失败报告的母报告hash已验证。
生产`src/`保持不变，未提交Git或Slurm任务。
