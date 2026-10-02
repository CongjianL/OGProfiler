# P5：固定真实 hierarchy 回归

状态：**P5 固定 hierarchy 真实策略差分与产物验收通过**（2026-10-01）。
性能/内存结论限于本次数据及下述计量范围，不表示 Orthobench 准确度已改善。

## 本次提交记录

- Job：`1410721`
- RUN_ID：`20261001T061912Z_f69c14cdfe9b_344c4dd4_17518`
- 基础 commit：`f69c14cdfe9bee9ed1aa251e89a511382b23baca`；dirty=1，精确源码快照和 diff 已保存。
- 源码 SHA256：`344c4dd4efb9d9be39380527f087fbf4bc33e2691ca12bc59ff793a5f9b45eab`
- 远端运行目录：`/home/mselab/licj/projects/running/ogprofiler-runs/20261001T061912Z_f69c14cdfe9b_344c4dd4_17518`
- 本地全套：270 passed；远端小型 fixture：14 passed；Ruff、76 文件 mypy、bash 语法检查通过。

本节在提交后补充，作业源码快照包含提交前版本的本说明；科学执行脚本/工具与快照一致。

## 真实结果验收

作业 `1410721`：COMPLETED，ExitCode `0:0`，耗时 `03:23:02`。
调度记录为 2026-10-01 14:17:50–17:40:52（Asia/Shanghai）。
`p5-audit/report.json` 的 `passed=true`、`artifacts_passed=true`、
`frozen_inputs_unchanged=true`、`artifact_error=null`。
本地重读全部 12,616 条组件诊断，八项检查每项均为 12,616/12,616：
events、ordered_members、member_multiset、unassigned、duplicates、remaining、
selection_sources、consumption；失败组件 0，首个分歧 null。

| 指标 | 结果 |
| --- | ---: |
| 输入蛋白 | 128,483 |
| SSN 组件 | 12,616 |
| 最大组件蛋白 | 34,982 |
| 层级节点（含 5,005 个 OG 输入适配的临时孤点根） | 19,353 |
| 持久化 terminal families | 16,582 |
| 最终 OG | 16,178 |
| 已分配 / 未分配 | 128,483 / 0 |
| 重复分配 | 0 |
| 最大 OG | 2,481 |
| >=1,000 蛋白的 OG | 5 |

五个大 OG 的大小为 1,366、1,378、1,411、1,568、2,481。它们与同层级 V1
参考策略一致；大小本身不证明生物学正确，也不作为删除或调参依据。
完整大小/物种覆盖分布在 report.json。

按组件的原始蛋白成员比较通过后，生产 OG stage 与 export 实际执行成功；
重读组件 membership hash、实际 TSV 原始成员集合及 unassigned 集合均与
纯 engine 结果一致，冻结的 SSN/输入/层级/配置校验和保持不变。
这验证的是 **同一份 V2 hierarchy 输入下的 V1 未精炼 OG 语义兼容**，
不是 V1/V2 独立构建 hierarchy 的等价，也不是完整历史 V1 程序的结果等价。

### 性能与内存观察

| 阶段 | Wall time | 最大 RSS（KiB） |
| --- | --- | ---: |
| components | 22:58.46 | 317,408 |
| hierarchy-all | 53:47.44 | 359,868 |
| annotate-network | 27:42.68 | 138,248 |
| P5 比较 + OG stage/export/校验 | 1:38:17 | 799,932 |

12,616 个组件的参考计时合计 31.939 秒（含参考视图构造/GML），
纯 engine 计时合计 9.404 秒；均开启 tracemalloc，不作独立公平性能 benchmark。
最大组件为 component 0：34,982 蛋白、1,817 节点；参考 Python 峰值
18,010,379 bytes（17.18 MiB），engine 峰值 9,881,182 bytes（9.42 MiB）。
最大 engine Python 峰值也发生在该组件。生产 engine 不存储内部节点完整
descendant lists，本次按组件检查未观察到全 HHN 成员复制；该证据只支持本数据规模。
整体 P5 进程最大 RSS 799,932 KiB（781.18 MiB）包含参考、全局 metadata、
stage/export 和原生库，不能用 9.42 MiB 替代其总内存。

纯策略比较耗时小，而端到端落盘/校验约 98 分钟；GNU time 记录 CPU 使用率
19%，存在明显非纯 engine 开销。具体 I/O、重复 checksum、metadata 读取与
export 成本仍需独立 profiling；本次不更改参数或缓存逻辑掩盖该问题。

小型验收文件已取回本地 `.provenance/slurm/20261001T061912Z_f69c14cdfe9b_344c4dd4_17518/results/`：
report.json、components.jsonl、fixed-origin.json、statistics.tsv、四份 time、
local-validation-summary.json、retrieved-files.sha256。大型数据与 hierarchy/OG/TSV 留在远端。

## 固定输入

使用作业 `1410705` 验证通过的 `v2-run`，SSN 算法
`of315-nbs-lrb-directional-mean-v4`，701,947 条无向 mean 边。
原目录只有 prepare/search/edges 产物，hierarchy 与 components 目录为空。
用户于 2026-10-01 同意：从该 SSN 构建一次 hierarchy 后冻结，再比较 OG。

`slurm/p5_fixed_hierarchy.sh VERIFIED_SSN_RUN` 在独立 Slurm run 下复制
input/edges/run.yaml，核验原 verification 和 SSN manifest 校验和；
不复制旧 checkpoint，不运行 prepare/search/edges。
依次执行 components、hierarchy-all、annotate-network、固定 hierarchy 审计。
科学参数直接读取原 run.yaml，无 overrides；仍为 seed=42、robust、adaptive、
max_depth=20、subtree_workers=1、runtime.workers=1。
资源沿用 SSN 验证脚本：cu/mselab/normal，56 CPU、128 GB、24 小时，无 array。
独立目录与源码快照由 `dev/slurm-submit` 创建，并记录 dirty diff、内容 hash 和 job ID。

## 回归契约

工具：`benchmarks/og_extraction/fixed_hierarchy.py`。

- 每次只读取一个组件的 hierarchy，保留原 Parquet 节点行序。
- 独立 postorder 适配器仅为 V1 参考构造完整 descendant gene/species 属性；
  V2 纯 engine 不使用这些列表、不构建全局 HHN。
- V1 参考仍执行固定 SHA 的原函数。复用已加载的 AST 函数减少组件间解析开销；
  参考后处理用 Counter 统计重复，未更改原 V1 函数。
- 两方使用同一蛋白 token；比较集合使用 `(species_id, original_id)`，避免同名冲突。
- V1 OG 函数只查询 SSN degree-zero。参考 SSN 视图精确保留显式 isolates，
  非孤点以链表示；这是 OG 输入适配，不是重建/比较 SSN 权重。
- 比较事件、顺序成员、组集合多重性、选择 source、后代消耗、剩余节点、
  未分配与重复。重复分类为 REFERENCE_OVERLAP，并使验收失败；不静默去重。
- reference/输入形态异常写入组件分类；首个分歧进入汇总。
- 仅纯策略比较全通过后运行 OG stage/export，重读组件 membership hash、
  实际 members.tsv 原始成员集合及 unassigned.tsv，与纯 engine 比较。
- 冻结 input/edges/components/hierarchy 与 run.yaml 的 SHA256，处理后再次核验。

输出：`p5-audit/report.json`、`components.jsonl`、`frozen-inputs.json`；
真实 OG/TSV 位于同一独立 run 的 `fixed-run/orthogroups`、`fixed-run/results`。
汇总包含 OG 大小/物种覆盖分布、最大 OG、预先设定的 >=1000 蛋白大 OG 数。
该阈值只用于诊断，不过滤候选或改变科学参数。

## 性能证据边界

每组件分别计量参考（含适配/GML）与纯 engine 时间和 tracemalloc Python 峰值；
作业阶段另有 GNU time 最大 RSS。进程 RSS 包含 V1 oracle、全局 metadata 与 export，
不是独立 engine native RSS；tracemalloc 也不覆盖原生库全部内存。
据真实结果报告组件规模/峰值关系，不把 Python 峰值直接当作全部内存证明，
也不把 export 全局排名内存当作组件纯 engine 内存。

## 验收

本地新增测试覆盖 13 个冻结参考 fixture 及实际 Parquet/TSV 流程。
真实回归必须 `report.json.passed=true`，Slurm COMPLETED 本身不等于通过。
若 hierarchy 构建或 network annotation 失败，先检查对应 stderr/manifest；
保留新 run 与固定 SSN，局部诊断，不重搜或盲目重复提交。
P6 Orthobench 准确度评价另行进行。
