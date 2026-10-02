# 失败节点搜索覆盖与准入冲突：审计及24点日程v2

## 冻结证据与诊断预算

固定`adr0004-hard-binary-24-v1`失败树，按其terminal membership重建6个失败节点的精确成员，
以同一mean SSN induced graph诊断；逐项验证edge manifest hash。
全部子图仍连通，不是缺边/断连造成的终止。method/seed/10轮/原准入门槛保持不变。

诊断是**独立有限覆盖审计**，不是对24点生产搜索暗中加预算：
每节点129个[.01,10]对数点；有count<2→>=2区间则额外65个该最早相邻区间的等距点；
如果最小gamma已经count>2，再33个[.00001,.01]点（显式BELOW_MIN_DIAGNOSTIC_ONLY），
重复点复用。总共1131点、3393次robust Leiden调用；小图最大132蛋白，本地执行。
不把未采样区间判为无解；测得连续多个点也不证明它们之间全部成立。

证据位于1410777本地结果目录：`binary-failure-coverage.json`，记录输入比较报告hash、
诊断脚本hash、精确成员、所有gamma/代表性membership/ARI/violations与调用数。
原递归报告不覆盖，修订日程输出独立`binary-recursive-comparison-v2.json`。

## 逐节点结论

|固定失败节点|n|诊断binary数|范围内合格binary数|判断|
|---|---:|---:|---:|---|
|2694:7|9|38|1|存在漏测合格点；其余37个binary均UNSTABLE|
|470:1|26|44|0|所测binary均MAX_CHILD_FRACTION|
|470:25|8|26|26|存在漏测合格点|
|84:0|132|7|0|7个合格binary全在gamma_min下方|
|844:0|33|6|0|所测binary均MAX_CHILD_FRACTION|
|500:0|45|59|0|所测binary均MAX_CHILD_FRACTION|

- **9蛋白**：gamma .9584374999999999合格，ARI=1、最大子群8/9。
  周边大量binary ARI=.25。旧日程仅探到.94后就转rescue，未进入窄窗口。
  这是搜索覆盖失败，不是binary目标和原准入必然冲突。
- **8蛋白**：26个合格采样点从.9690625到.976875，ARI=1、最大子群7/8。
  旧五步探测停在.98，漏测这些点。
- **26/33/45蛋白**：binary分别为25/26=.961538、32/33=.969697、44/45=.977778，
  均超过.95；所测点不是靠增加同类探测次数就可接受。
  这不是整个区间无合格binary的证明，只是全部诊断binary的实际准入冲突。
- **132蛋白**：范围内全部采样count>=3；下界外7个合格采样gamma
  .0004869675.. .0017782794、最大子群比例.75。
  说明binary在另一搜索范围存在，不证明当前整个范围无binary。
  本轮不将下界改成小值，也不把下界外候选计入验收。

## 修订日程

完整规范见[ADR0004 v2](adr/0004-small-node-binary-target.md)。共享预算仍为24点/72调用：

`COARSE 8 → ENDPOINT 1 → RESCUE 3 → TARGET_PROBE 9 → REFINE 3`

仍失败时显式将3个refine名额借给target，总目标探测最多12点。
减少高gamma coarse重复count信息；global rescue提前供区间选择使用；
失败binary的上下边界轮转探测，而非持续选择最低gamma区间。
没有target边界时以最宽log区间的几何中点补覆盖。
未加入任何已知节点专属gamma，没有降低阈值或隐式k-way回退。

单侧“只追上边界”会漏掉2694根下方的合格binary，因此最终规范使用双侧轮转。
双侧轮转仍是有限启发式，可能在狭窄稳定窗口上花费不足。

## 修订后完整递归结果

同5个完整组件，从根重新构建两种hierarchy；每个hard-binary再独立串行重复一次。

|组件|v1失败节点尺寸|v2失败节点尺寸|v2节点数|v2候选/调用数|
|---|---|---|---:|---:|
|2694|9|9|13|125 / 375|
|470|26，8|26|41|334 / 1002|
|84|132|132|1|24 / 72|
|844|33|33|1|24 / 72|
|500|45|45|1|24 / 72|

- 470的新树node25（8蛋白）在gamma .97选中合格binary，后代完整解析；该失败已消除。
- 2694根仍选gamma .6200000000000001。其9蛋白后代使用24点，测到多个binary但均UNSTABLE；
  最后测.95875为3-way，仍漏掉已知合格.9584375。本轮不继续针对这个数字调参。
- 84/844/500及470:1仍失败，分别与范围/比例冲突证据一致；不宣称因此证明无解。
- 5组件仍**0/5 fully resolved**，均保持BLOCKED_UNRESOLVED，未进入OG/scorer；无RefOG F1结论。
- 271蛋白覆盖完整且唯一，所有结果通过结构验收；每节点<=24点/72调用，调用数=候选数×3。
- 重复运行nodes/members/candidates一致；baseline节点、成员、OG分区再次匹配冻结1410777。
  尚未实现生产并行接线，不将串行重复性作为串并行一致性。
- 本机v2单次hierarchy约.64/1.14/1.22/.34/.63秒，不外推大组件成本。

## 后续边界

v2是明确版本化的实验日程，不满足完整递归验收，暂不接线生产或跑大组件H4/H5。
9蛋白问题是窗口覆盖，26/33/45是实测比例冲突，132是实测搜索范围冲突；
单纯重分配24点并未同时解决三类问题。
下一政策决策须分别明确：进一步有限覆盖如何分配、是否显式修订gamma下界、
以及是否引入标记清晰的soft binary→k-way回退。
停止尺寸1、ARI .9、最大子群比例.95、OG/scorer继续保持不变。

## 复现

```sh
R=.provenance/slurm/20261002T060705Z_cb0c0ba98049_1cb4fcd7_14702/results
PYTHONPATH=src:. .venv/bin/python benchmarks/og_extraction/binary_failure_audit.py \
  --result-root "$R" --comparison "$R/binary-recursive-comparison.json" \
  --output /tmp/binary-failure-coverage.json
PYTHONPATH=src:. .venv/bin/python benchmarks/og_extraction/binary_recursive_audit.py \
  --protocol adr0004-hard-binary-24-v2 --result-root "$R" \
  --output /tmp/binary-recursive-comparison-v2.json
```

## 验证状态

新增诊断测试2个、v2日程测试5个；结合既有binary/refog诊断与hierarchy/bounded-stage测试，
**47个针对性测试通过**（3.76秒）。Ruff及git diff --check通过。
生产`src/`没有改动；未提交Git或Slurm任务。
