# ADR0007：将已验证的 soft42 设为 V2 默认配置

- Status：Accepted（2026-10-06，用户明确要求）。
- Supersedes：ADR0004/0005中的仅opt-in及全局kway/depth20限制，和ADR0006本轮维持默认的范围限制。
  历史实验和图2对照保持原始含义，不回写历史结果。

## 决定

公共V2配置入口 DEFAULT_CONFIG / load_config 的完整默认配置与已冻结的
`presets/embleya-soft42.yaml`一致。实际差异仅两项：

- hierarchy.topology_policy：kway_v1 → soft_binary_24_v2。
- hierarchy.max_depth：20 → 42。

其余配置全部保持：seed42、10轮、nonempty_children_v1、bounded_adaptive_v2、
ARI .9、比例 .95、停止尺寸1、每节点24候选、v1_compatible OG、refinement=false、
search.threads=8、runtime.workers=1、subtree_workers=1。不把集群56-worker运行覆盖
混入生产默认，也不引入实验图评分模型。

## 依据与局限

冻结Embleya13基因组相对OF参考：P95.7133/R69.0770/pairF1 80.2424%。
Open Orthobench job1411546：P92.6961/R33.8220/F1 49.5608%；此前默认kway20
job1410777为P91.7440/R25.5911/F1 40.0192%。完整树/OG parity/输入与官方评分检查通过。
这是用户选择的默认策略，不是通用最优证明；soft42的Orthobench F1仍低于完整V1
57.5075%，召回不足仍存在。不同参考体系的指标不直接横比。

## 兼容性

新建CLI运行、无配置及稀疏YAML的未指定字段继承soft42。完整显式配置与CLI override
优先级不变。已有run.yaml若省略topology/depth也会继承新默认；历史重放须显式固定
旧字段，不能把继承新默认的运行当作原配置重放。旧策略可以通过两项显式覆盖继续使用：

```sh
--set hierarchy.topology_policy=kway_v1 --set hierarchy.max_depth=20
```

不批量修改已有run.yaml或远端运行、不删除缓存、不自动重算历史数据。
底层实验构造器 HierarchyConfig()/ResolutionSearchConfig()保留原兼容默认；
公共V2生产入口通过hierarchy_config(load_config()["hierarchy"])得到新配置。
这些低层构造器的历史测试不等同于V2生产配置入口。

## 验收

增加整个默认配置与冻结preset逐项相等的回归测试；检查partial YAML继承及旧策略
显式覆盖；旧depth实验测试显式固定depth20，保留原实验语义。运行完整轻量测试，
不为默认配置发布自动追加科学计算。

历史gamma_strategy矩阵中的legacy log_grid分支需显式固定kway以满足兼容约束；
该分支比较的是整套搜索政策（strategy/admission/topology），不是单独gamma字段。
保持显式10轮及既有矩阵数量，不运行新矩阵实验。legacy_strict稀疏配置也须
显式指定topology_policy=kway_v1；不静默改写不兼容的科学参数。

验收结果：完整本地416 passed / 8 skipped，5项既有warning；远端通过标准wrapper
进行无参数默认配置与冻结preset逐项相等的短验证。历史数据、图2和Slurm运行不重算。
