# ADR 0003：显式提高有限 Leiden 优化迭代预算

- Status：Accepted（2026-10-02，用户确认优先修订优化预算，停止尺寸与稳定性阈值不变）
- Amends：ADR 0002 的每次 Leiden 迭代预算；准入、搜索日程、失败与发布语义不变。
- Evidence：`../recursion-stop-and-818-node-audit.md`，H4/H5 job `1410775`。

## 决策

新策略 `nonempty_children_v1 / bounded_adaptive_v2` 的每次 Leiden 调用统一使用
有限 `leiden_iterations: 10` 默认值，替代2。不是无限收敛模式，也不是失败后追加
第二轮搜索。CLI 显式覆盖与已保存配置中的显式数值继续生效。

保持 `recursion_stop_size: 1`、`stability_threshold: 0.9`、robust 三 seed、
`max_candidate_evaluations: 24`、gamma 日程、最大子群比例、目标函数与 mean SSN
权重不变。每节点 robust 仍至多72次调用，每次至多10轮优化；component call cap
仍按调用数而非优化轮数计费。该预算不保证收敛、可分性或全运行墙钟时间。

## 依据与成本

818 蛋白节点在2轮下是19个候选全部被拒绝，不是候选额度耗尽。固定图、seed、
阈值的公共搜索对照在10轮下接受 gamma 1.0763474115247547、44子群、mean ARI
0.9131750774646644：11候选、33调用，而2轮是19候选、57调用。
本地对照耗时约3.85秒与1.96秒，调用减少不等于耗时减少。
这只是该节点的证据，不证明完整 component0 已解决或 Orthobench 准确度改善。

暂不增加失败后第二轮回退：旧19点加新搜索可能超过24候选，而且同gamma不同
优化预算应视为不同评估。未来回退须另修订版本、共享总候选/调用预算，并将
迭代数纳入候选键；本轮保持现有单轮日程与失败状态。

## 兼容与缓存

- 新 YAML/CLI 配置省略迭代数时默认10。legacy_strict 配置省略该键时保持历史2；
  显式2或10均保留，覆盖后的最终政策决定省略值。
- 公共 `HierarchyConfig()` 默认10；直接构造 legacy 回放配置应显式传入
  `leiden_iterations=2`。低层 `run_leiden`/counter 的历史默认2保留，生产路径
  总是从 HierarchyConfig 显式传入预算。
- 历史 run.yaml 的显式2不静默重写；回归迁移工具显式改为10并记录来源配置。
- OFAT 搜索策略对照显式固定10轮，避免切换 legacy 时同时改变优化预算；
  matrix 配置哈希通过同一配置加载接口计算。
- hierarchy 算法版本升级至 `hierarchical-leiden-v3`，scheduler manifest 至v3；
  参数身份继续绑定迭代预算与环境。Parquet schema保持2：列语义没有改变。
- OG/scorer、停止尺寸、ARI阈值保持原样；旧诊断不按新版本直接resume。

## 验收

通过配置默认/显式覆盖/legacy回放、迁移预算、串并行等价、实际阶段预算
落盘及变更预算失效的本地契约测试。全 component0 的H4、通过门槛后的H5仍须
显式提交新的不可变 Slurm 快照验证；本次实现不自动提交。

### 本地实施验收（2026-10-02）

- 完整测试307 passed，18个已有警告；Ruff、mypy（77文件）、diff检查通过。
- 默认10、旧配置显式2、CLI覆盖与省略值的legacy兼容测试通过；历史API回放fixture显式固定2。
- 迭代预算变更导致resume失效，实际manifest记录10；v2缓存与v3隔离，schema仍为2。
- 冻结818图在新公共默认配置下复现11候选/33调用、44子群、ARI .9131750774646644。
- 尚未执行新预算的完整H4/H5，未提交Git改动或推送。

### 后续真实数据验收

commit `cb0c0ba` 的clean不可变快照通过job `1410777` 完成H4/H5。
component0及全运行unresolved为0，完整覆盖、预算、稳定性、串并行和固定hierarchy
V1 OG parity通过。新正式OG官方P/R/F1为91.7440% / 25.5911% / 40.0192%；
工程修复得到支持，但仍未达到历史完整V1的57.5075% F1。
详见 `../h4-h5-iterations10-regression.md`。停止尺寸、阈值和OG/scorer保持不变；
当前无未解决节点，不启动额外回退或重提作业。
