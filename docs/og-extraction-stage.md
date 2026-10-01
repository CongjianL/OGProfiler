# P3：V1-compatible OG 组件落盘与恢复

## 阶段接口

```bash
ogprofiler orthogroups --run RUN_ROOT
ogprofiler orthogroups --run RUN_ROOT --set runtime.workers=2
ogprofiler run --out RUN_ROOT --from-stage orthogroups --until-stage orthogroups
```

阶段顺序现在为 `SSN → hierarchy → annotate-network → orthogroups → export`。
**本阶段只发布组件 OG 数据；最终 export 当前仍读取 terminal family，P4 才切换。**
`run --until-stage orthogroups` 可用于单独验收，避免把旧 terminal 导出当成新 OG。

配置默认值（支持 YAML 和 `--set`；省略 `--config` 时读取已有 `run.yaml`）：

```yaml
orthogroups:
  strategy: v1_compatible
  species_overlap_count: 0
  refinement: false
runtime:
  workers: 1
  component_retries: 1
```

整数 overlap 参数拒绝负值、浮点数和布尔值；当前只有 `v1_compatible` 策略。
`refinement=true` 明确报错，避免默默采用未精炼策略。
配置及单组件纯引擎继续分离，worker 内部没有嵌套进程池。

## 输入与 isolate 契约

- `input/proteins.parquet`：`protein_id, species_id, original_id`。
- `input/species.parquet`：完整数据集物种数，不使用组件物种数替代。
- `components/index.parquet`：完整组件成员分区。
- `components/singleton_terminal_families.parquet`：显式 SSN isolate 证据，允许空表。
- 非孤点组件的 `hierarchy/components/component=.../{nodes,members}.parquet`。
- 有 hierarchy 的组件绑定 `evolution/components/component=.../events.parquet` 校验和。
  该表是上游诊断依赖，选择过程重新计算专用 v1_event，不使用 network_event/legacy_event
  作为 V1 选择标签，也不改变这些标签的 schema。

检查蛋白分区、species metadata、跨组件原始身份唯一性、singleton 来源；
单组件继续检查 hierarchy 的根/parent/depth/计数/terminal membership。
没有 hierarchy 的真实 SSN isolate 使用临时单节点输入视图，按 P2 的 SSN_ISOLATE 分支选择。
该视图不会写入 hierarchy；普通 hierarchy singleton 不自动获得 isolate 身份。

## 产物

每个组件目录 `orthogroups/components/component=00000000/`：

| 文件 | 内容 |
| --- | --- |
| `groups.parquet` | component/local_group/source_cluster、选择类型、v1_event、处理等级、基因/物种数、membership_hash |
| `members.parquet` | component_id、local_group_id、protein_id |
| `selection_trace.parquet` | trace_order、候选 cluster、等级、selection_event、状态、消耗来源 |
| `unassigned.parquet` | protein_id、terminal_cluster_id、原因；即便无未分配也发布有 schema 的空表 |
| `v1_events.parquet` | 原始 V1 事件、reference_order、原始无向度、eligible 子节点数；Python None 为 Parquet null |
| `og-manifest.json` | 仅完成产物；算法/schema/事件版本、配置、输入/输出 hash、统计、remaining_cluster_ids |
| `og-failure.json` | 失败诊断；异常类型/说明、失败时身份，重复分配时附蛋白及 duplicate 统计 |

成功重建会清理该组件的 failure marker。异常 OG 候选不作为正常结果发布。
每次 Parquet 写入最多缓存 65,536 行；不把所有成员展开成一次性 row list。
worker 读取过滤后的组件 index/protein/isolate 数据，以及该组件的 nodes/members。
调度进程读取并校验全局轻量 metadata/index，不构造全局 igraph 或完整 hierarchy。

`membership_hash` 使用 P2 的排序 `(species_id, original_id)` 契约，local_group_id 使用确定性选择顺序。
串行和并发结果中两者及全部 Parquet 行一致；**最终全局 OG ID 的分配和导出仍在 P4**。
本阶段没有重新定义终端 membership，也没有调整 SSN、Leiden、分裂或 pairwise 策略。

## 恢复与失败

- 算法 `component-orthogroups-parquet-v1`，schema 1，绑定纯引擎和专用事件算法版本。
- resume 同时核对身份、DONE manifest、五份产物校验和和 Arrow schema。
- overlap 或全局 protein/species/index/isolate 文件变化保守地重建全部组件；
  单组件 nodes/members/events 变化只重建该组件。
- workers/retries/命令变化不改变科学身份，不强制重建已验证产物。
- 输出损坏、遗漏、manifest 损坏或版本变化触发重建，不凭文件存在宣称完成。
- 临时目录写完后逐文件 `os.replace`，最后原子发布 completion manifest。
  **不是整目录事务**；中途失败会留下部分数据文件，但没有可复用的 DONE marker。
- 运行前移除旧的阶段总 manifest；组件失败不阻止其他组件完成。
  总 `orthogroups/og-manifest.json` 记录 DONE/FAILED、组件清单及 manifest hash。
  读取时应以清单为准，忽略历史残留、已不属于当前 index 的目录。
- scheduler 独占 SQLite 转换；worker 只写自身组件产物。
  stale RUNNING 自动恢复；即便数据库状态落后，完整验证的产物也可恢复为 DONE。
- `component_retries` 指本次调用内的额外尝试次数；跨调用保留 attempts 审计记录，
  但不因历史失败次数耗尽而永久阻止用户恢复已修复的组件。
- 成功组件在后续恢复中复用；FAILED 总 manifest 及 SQLite 失败在 `status` 中明确显示。
- 上游全局预检失败不发布完成 manifest；损坏输入通过 CLI 错误返回。
- 同一 workspace 只运行一个 orchestrator，执行期间输入应保持固定；本阶段不实现多写者锁。

## 验收范围

冻结的 13 份 V1 fixture 经 Parquet 落盘后，成员集合、处理等级和 unassigned 匹配独立 oracle；
参考中的重复分配场景被明确归入 FAILED，记录 duplicate 数而非静默去重。
另覆盖配置、所有产物/manifest 损坏、输入 checksum 改变、输入 metadata 错配、
临时写入和部分发布失败、同次 retry、跨调用恢复、stale RUNNING、串行/双 worker 及 CLI stage range。
这些是小型本地工程验收；未执行 Slurm/真实 hierarchy 回归或 Orthobench 准确度评价。

本地验证记录：新增 P3 测试 **44 passed**，完整测试集 **225 passed**（14 条已有 SciPy
拟合 covariance warning）；Ruff、新 orthogroups 包 Mypy、`git diff --check` 均通过。
完整测试使用本地已固定 OF3.1.5 参考源码：

```bash
OF_SOURCE_ROOT=/private/tmp/ogp-ssn-audit-of315/src/orthofinder .venv/bin/python -m pytest -q
.venv/bin/mypy src/ogprofiler/orthogroups --follow-imports=silent
```

开发环境的已安装 wheel 可能仍为旧版本；未重新安装时，源码 CLI smoke 使用：

```bash
PYTHONPATH=src .venv/bin/python -m ogprofiler run --out RUN_ROOT \
  --from-stage annotate-network --until-stage orthogroups --dry-run
```
