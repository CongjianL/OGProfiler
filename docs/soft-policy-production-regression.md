# ADR0004 soft policy：生产 opt-in 接线验收

日期：2026-10-03。默认保持 `kway_v1`，本轮不提交 Slurm 作业。

## 启用及身份

```yaml
hierarchy:
  topology_policy: soft_binary_24_v2
```

也可对 hierarchy / hierarchy-all 使用
`--set hierarchy.topology_policy=soft_binary_24_v2`。
默认停止尺寸1、ARI .9、比例 .95、gamma [.01,10]、seed42、robust3、优化10轮不变。
生产模块 `src/ogprofiler/hierarchy/soft_search.py` 实施冻结 soft-v2 日程，
串行和 spawn workers 通过同一个公共 search 调用，不使用审计 monkeypatch。

algorithm `hierarchical-leiden-v4`、scheduler v4、schema3；完整配置身份含 topology。
旧schema2/v3失效；opt-in/default切换失效；同配置新结果 resume 通过。
新增原门槛 eligibility/violations、selection_kind 与 refinement_truncated 落盘。
manifest 有 fallback_count/refinement_truncated_count。回退不额外增加调用。

## 五组件完整递归

固定 job1410777 meanSSN 小组件输入；公共串行与真实 spawn 双进程分别完整递归并落盘。
未覆盖源结果，独立输出位于本地 `.provenance/production-soft-optin-20261003-verified/`。
报告为 `report.json`，完整配置、输入与源hash、Git dirty状态为 `provenance.json`。

|组件|蛋白|节点|每次运行调用数|fallback|串行秒|双进程秒|
|---|---:|---:|---:|---:|---:|---:|
|2694|15|28|594|1|0.729|1.341|
|470|46|82|1656|1|1.884|1.649|
|84|132|238|4713|10|7.720|5.198|
|844|33|51|900|5|1.671|2.074|
|500|45|83|1707|1|1.937|1.995|

五组件全部：
- unresolved=0；271蛋白唯一完整覆盖，无新增或丢失成员。
- 每节点最多24 unique gamma，Leiden调用严格等于候选数×3。
- 串并行 nodes/members/candidates 及调用数完全一致；三张Parquet表完全一致。
- 对冻结实验原字段逐项比较：候选、拓扑与终端成员完全复现。
- 实际冻结V1 OG适配器 events/ordered_members/member_multiset/unassigned/duplicates/
  remaining/selection_sources/consumption 八项全部通过。
- 全部18次fallback带显式kind/phase及原准入证据。

小组件的进程启动开销明显，耗时仅为本机开发诊断，不代表大型组件加速比例。
此批是开发选中的低召回反例，不是官方全数据P/R/F1。

## 测试

新增生产测试覆盖：默认/opt-in、非法策略与>24预算、原UNSTABLE/比例失败保留、
full-cap rejection与小cap exhaustion、binary优先、refine截断、2蛋白原k-way、
Parquet诊断、配置切换失效与schema失效、soft H4公共CLI串并行验收。
全量本地测试342 passed、8 skipped；跳过项依赖外部OF源码或官方scorer路径。
Ruff、diff whitespace与Slurm shell/内嵌Python语法检查通过。既有OG/scorer代码无变更。

## 后续 H4/H5（准备，尚未执行）

现有 `slurm/h4_h5_hierarchy_regression.sh` 增加可选第3参数 `soft_binary_24_v2`。
使用P6原始meanSSN来源与job1410777的冻结10轮H4结果作为对照；迁移工具显式写入
soft配置，脚本验证新旧配置仅 topology 政策不同，逐项核对固定输入hash。
默认第3参数仍kway_v1，保留原2→10迭代预算回归路径。

资源请求保持原56CPU/250G/72小时，无array。通过dev/slurm-submit不可变快照执行。
H4验收完整component0及双进程回放：成员、预算、原门槛fallback证据、调用数、耗时。
H4通过才执行其余组件H5、冻结hierarchy、V1 OG parity及官方Orthobench。
若出现unresolved，保留诊断并阻止OG/评分；不调低停止尺寸、ARI或比例门槛。
