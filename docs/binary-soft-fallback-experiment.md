# 第3项：上边界补探与显式soft binary回退（实验预定义）

Status：Proposed，仅独立实验；生产与既有hard-binary版本保持不变。
这是小节点拓扑政策的显式修订，不仍称hard-binary完成。

## 两个分开的实验身份

- `adr0004-coverage24-upper-guard-v1`：hard binary +上边界补探。
- `adr0004-soft24-upper-guard-v1`：相同日程，但完成hard binary搜索后允许标记k-way回退。
- 对照原A，各在gamma_min=.01/.001进行五组件完整递归；范围与回退是两个独立因素。
  每个实验单独24点，不把不同实验调用相加称一次预算。

## 共享24点日程

沿用A的early coarse、条件global点、公平边界和失败时借用3个refine名额。
新增：coarse停止在被拒绝raw binary、且没有已测更高点时，在target前最多3次
`min(previous_gamma*growth_factor, gamma_max)`上探，trace=`UPPER_GUARD`。
若已找到合格binary则停止；若已找到原门槛合格k-way则停止上探，但仍进入binary目标搜索。
全部上探从同一24点中扣除；不追加一轮预算，不在失败后偷偷调用原搜索。

## 回退语义

- 合格binary优先，仍按实际测得最低gamma选择。
- 完成有限binary搜索、无合格binary时，soft版本仅复用已测点中的**原门槛合格k-way**，
  选择最低gamma；不追加Leiden调用，不放宽ARI、比例、quality或结构门槛。
- 显式更小cap阻断日程时不触发回退，以EVALUATION_BUDGET_EXHAUSTED保持失败。
- 回退phase=`FALLBACK_KWAY/<原评估phase>`；报告分别列binary/kway eligibility、
  原门槛violations、selected及fallback节点。不得用同一个valid标志掩盖准入政策改变。
- soft选中k-way后属于SPLIT，继续完整后代递归；没有合格binary或k-way仍UNRESOLVED，
  完整成员保留，阻断OG。NO_EDGE_SUPPORT/DEPTH_LIMIT等失败不通过该回退变成终止叶。
- 无失败组件才调用原纯OG引擎；其公式与scorer不改。解析成功不是准确度改善证据。

验收：独立重复、唯一成员覆盖、<=24点/72调用、fallback追溯、所有后代resolved及OG门禁。
本实验不批准生产迁移、缓存/schema接线或Slurm大规模运行。

## v1结果与v2补探规范（第二轮运行前固定）

v1 soft在两个下界下均4/5 resolved；45蛋白的首次coarse split是原门槛不合格3-way，
因此“仅raw binary触发上探”漏掉其上侧。完整候选仍保留，未修改原参数或v1身份。

新增明确身份 `adr0004-coverage24-upper-guard-v2`、`adr0004-soft24-upper-guard-v2`：

- 上探触发改为coarse最后原始count>=2且**原门槛不合格**，不限binary。
- 上探直到原门槛合格的分裂、gamma_max，或共享24点只剩3点时停止；
  最后3点留给目标探测（不是额外3点）。更小显式cap也按同一保留规则处理。
- 上探可找到合格binary而直接成功；或找到合格k-way仅作为潜在soft候选，仍先搜索binary。
- 沿用其他规则、失败状态、候选复用和明确fallback标记；不为45蛋白插入特殊gamma。
- 固定原下界.01，从根重跑全部五组件；不同时再降低下界。

## 完整递归结果（2026-10-03，北京时间）

固定gamma_min=.01的v2三策略对照，原五组件从根完整重建，每种配置独立重复一次。
生产下界仍.01，没有同时降低下界、停止尺寸、ARI或最大子群比例。

|组件|A调用数/状态|仅v2上探调用数/状态|soft-v2调用数|fallback节点数|soft-v2 OG数|
|---|---|---|---:|---:|---:|
|2694|546 / unresolved 4|561 / unresolved 4|594|1|3|
|470|843 / unresolved 26|843 / unresolved 26|1656|1|6|
|84|72 / unresolved 132|72 / unresolved 132|4713|10|36|
|844|72 / unresolved 33|72 / unresolved 33|900|5|14|
|500|72 / unresolved 45|72 / unresolved 45|1707|1|10|

- soft-v2 **5/5完整resolved**，A与仅补上边界的hard版本仍0/5。
  这明确把解析改善归因于**显式拓扑准入政策改变**，不宣称仅靠搜索覆盖修复hard binary。
- 18个fallback节点全部有`FALLBACK_KWAY/<phase>`标签，所选候选均k-way eligible、
  binary ineligible，原门槛violations为空；未追加调用、降低门槛或将失败普通终止。
- 271蛋白完整唯一hierarchy覆盖，每节点<=24点/72调用；soft-v2单次五组件合计9570调用=3190候选×3（重复运行与对照调用另计）。
  完整树全部结构验收通过，重复nodes/members/candidates一致。
- 同一原OG引擎输出69个OG、271个唯一成员、unassigned=0；最大OG尺寸按组件为13/15/17/15/13。
  不把这些小组件尺寸结论外推成全数据大组污染验收。
- 原baseline节点、成员、OG分区再次与1410777冻结结果一致。
- v1在两下界均4/5 resolved；45蛋白根最高已测gamma仅.04/.032。
  v2按通用rejected-split条件上探，本例实际测到.64合格6-way，再在binary日程完成后显式回退。
  未插入45蛋白专属gamma，旧v1结果保持原文件。

## 同新hierarchy的实际V1 OG parity

使用已有冻结V1 AST执行适配器，而不是以生产事件/OG结果构造oracle，
对soft-v2五份新hierarchy分别执行原V1事件与ExtractOGSorted。

**5/5 MATCH**：raw events、ordered members、member multiset、unassigned、duplicates、
remaining nodes、selection sources、consumption八项全部一致，无重复成员或首次分歧。
这是同hierarchy的OG parity；不表示本拓扑与V1原始binary搜索拓扑相同。

## 局部RefOG关系诊断（不是官方全数据P/R/F1）

只计五个选中SSN组件内部、同组件的confident RefOG真基因对；不计跨组件或未选组件。

|RefOG|局部真对数|冻结baseline保留|soft-v2保留|
|---|---:|---:|---:|
|001|105|0|78|
|012|1035|231|231|
|035|435|3|97|
|036|990|2|172|
|047|528|0|116|
|合计|3093|236|694|

RefOG012仍未改善，是“binary或完整解析不自动保证全部召回”的实测反例。
这些组件按已知问题选择，属于开发对照，存在选择偏倚；不报告官方P/R/F1、全数据coverage
或benchmark superiority。OG/scorer本身未修改。

## 证据与后续门槛

1410777本地results：

- `binary-soft-fallback-comparison.json`：旧v1三策略×两下界实验，保留不覆盖。
- `binary-soft-fallback-expanded-comparison.json`：v2原下界三策略完整递归，
  保存源hash/源内容快照、Git commit与tracked diff、完整config及候选eligible/fallback明细。
- `binary-soft-fallback-expanded-v1-parity.json`：冻结V1 reference hash、母报告hash及八项检查。
- `binary-soft-fallback-expanded-refog-pairs.json`：母报告/truth hash及严格局部真对关系计数。

本轮满足小组件完整递归与同hierarchy V1 OG parity，但仍**Proposed**。
下一步需要显式确认soft政策，再接线配置身份、缓存/schema/fallback trace与生产串并行实现；
完成相同小组件的serial/parallel验收后，再通过不可变Slurm H4/H5做全组件与官方Orthobench。
尚未通过生产串并行、全component0、全数据accuracy或大组污染验收。
生产src未改，未Git提交或Slurm提交。

```sh
R=.provenance/slurm/20261002T060705Z_cb0c0ba98049_1cb4fcd7_14702/results
PYTHONPATH=src:. .venv/bin/python benchmarks/og_extraction/binary_soft_fallback_comparison.py \
  --expanded --result-root "$R" --output /tmp/binary-soft-fallback-expanded.json
```

## 测试状态

61个针对性测试通过（新增soft/guard测试5个），Ruff与git diff --check通过。
未提交Git改动，未推送或提交Slurm。
