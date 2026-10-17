# ADR0004：共享24点日程的小组件完整递归对照

## 范围与冻结条件

2026-10-02，本地独立实验 `adr0004-hard-binary-24-v1`，完整日程见
[ADR0004](adr/0004-small-node-binary-target.md)。仍为Proposed；未切换生产默认。
这是对现有ADR0001小节点k-way准入的**显式实验修订**，不是已接受的生产变更。

使用Slurm 1410777不可变输入的五个小组件，SSN仍为方向权重mean。
method=RBER、seed=42、robust三seed、mean ARI>=.9、最大子群比例<=.95、
每调用10轮、停止尺寸1、max_depth20全部沿用冻结manifest。
仅小多物种3..9999节点替换搜索日程与exactly2准入；2蛋白及>=10000沿用生产搜索。
使用现有DFS、结构验收与纯OG引擎，不修改生产源文件、OG或scorer。
实验adapter仅在该进程的串行上下文中替换搜索入口，不是生产/并行接线。

证据目录（未提交的大体积诊断）：
`.provenance/slurm/20261002T060705Z_cb0c0ba98049_1cb4fcd7_14702/results/binary-recursive-comparison.json`。
报告包含完整config、生产源hash、实验脚本hash、所有节点/候选/成员、调用数、失败与门禁状态。
组件edge输入逐项核对冻结manifest的sha256；未重建SSN。

## 完整递归结果

|组件|蛋白数|baseline节点/调用数|hard-binary节点/调用数|失败节点（新树cluster_id:尺寸）|最大候选数/节点|
|---|---:|---:|---:|---|---:|
|2694|15|18 / 102|13 / 342|7:9|21|
|470|46|62 / 627|27 / 639|1:26，25:8|21|
|84|132|184 / 1989|1 / 48|0:132|16|
|844|33|37 / 129|1 / 63|0:33|21|
|500|45|57 / 423|1 / 63|0:45|21|

- baseline五组件均resolved；节点拓扑/状态/选中gamma、terminal membership、最终OG分区
  与冻结1410777一致。quality比较容忍1e-12数值差异，忽略非存储event字段；不宣称Parquet字节一致。
- hard-binary五组件均包含UNRESOLVED，**0/5整组件resolved**。
- 两方案均覆盖各自输入的全部271蛋白，成员唯一且无遗漏；每份结果通过validate_hierarchy。
- hard-binary分别再次完整串行重跑，nodes/members/candidates完全一致。
  这证明本地重复确定性，不等同于尚未实现的串并行一致性验收。
- 实际Leiden调用数严格等于候选数×3；最大21点/63调用，全部在24点/72调用共享上限内。
- 本机首轮hierarchy耗时baseline约.17/.81/2.77/.44/.75秒，hard-binary约.50/.81/.77/.29/.37秒。
  失败树提前停止，耗时减少不是完整解析效率改善；小输入本地时长不外推Slurm大组件。

## 失败的具体含义

全部为`REJECTED_ALL_TESTED`，没有共享预算阻断新请求的`EVALUATION_BUDGET_EXHAUSTED`。
失败日程没有REFINE，最多用21点；剩余3点是成功后refinement配额，不自动补搜。

1. **2694**：新根在gamma .6200000000000001选中合格binary。
   上游拓扑已变化，原2694:3的11蛋白节点不再是验收对象。
   新9蛋白后代测试21点，原始count包含1/3/5/9，没有测试到2；保留UNRESOLVED。
2. **470**：根仍是gamma .08的binary，反例并未消失。
   新26蛋白后代的5个binary均违反MAX_CHILD_FRACTION；另一个8蛋白后代20点未测到binary。
   rescue与target探测有重复点，所以该后代仅20个unique gamma。
3. **84**：根在最小gamma .01已经是3-way，测试点的原始count均>=3。
   没有符合规范的count跨越/二群端点区间，因此TARGET_PROBE阶段为空；
   10 coarse + endpoint + 5 rescue共16点完成后失败。
   这只说明本有限日程未命中binary，不证明整个[.01,10]内无binary，更不改变下界。
4. **844**：3个binary候选均违反MAX_CHILD_FRACTION（32/33）。
5. **500**：3个binary候选均违反MAX_CHILD_FRACTION（44/45）。

所有UNRESOLVED保留完整成员，报告为`BLOCKED_UNRESOLVED`，没有将结构叶转为普通terminal，
没有k-way回退、降低ARI或放宽最大子群。实验明确执行门禁后才调用纯OG引擎；
五份hard-binary都未进入OG/正式评分，故不报告其RefOG recall/F1，也不以失败叶分区替代正式OG。
baseline OG与已通过固定hierarchy V1 parity的冻结分区一致；本次没有新的实验hierarchy OG parity成功结论。

## 下一步判断

本日程不满足完整递归验收，暂不生产接线或启动全组件H4/H5。
先分别审计9/8蛋白新后代的未采样区间、84根的下界/非单调count覆盖，
以及26/33/45蛋白binary的比例冲突；区分覆盖不足与准入冲突。
若重新分配失败后的3点，须新protocol明确共享预算与区间优先级；
若引入soft binary→k-way回退，须另行修订ADR并标记fallback，不仍称hard-binary完成。
停止尺寸1、ARI .9、OG/scorer继续保持不变。

## 复现与测试

```sh
PYTHONPATH=src:. .venv/bin/python benchmarks/og_extraction/binary_recursive_audit.py \
  --result-root .provenance/slurm/20261002T060705Z_cb0c0ba98049_1cb4fcd7_14702/results \
  --output /tmp/binary-recursive-comparison.json
PYTHONPATH=src:. .venv/bin/python -m pytest -q \
  tests/unit/test_binary_recursive_audit.py \
  tests/unit/test_binary_target_audit.py tests/unit/test_refog_split_audit.py
```

12个诊断测试通过，新增6个测试覆盖：count目标中点探测、原门槛拒绝保留、
共享预算阻断、尺寸边界、失败整树成员保留、24点用满但仍ACCEPTED。

另外现有hierarchy单元与bounded-stage集成测试28个通过；合计40个针对性测试通过，
新增实验代码与测试Ruff及git diff --check通过。未提交Git改动，未提交Slurm任务。
