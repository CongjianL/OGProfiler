# P4：最终 OG 导出与独立终端诊断

默认 `ogprofiler export --run RUN_ROOT` 读取已完成的 `orthogroups` membership，
而不是把 hierarchy terminal family 当成最终 OG。内部节点可合并多个终端 family；
terminal membership 仍是原来的全蛋白分区，既不修改也不重新聚类。

```bash
ogprofiler run --out RUN_ROOT --from-stage annotate-network --until-stage export
ogprofiler export --run RUN_ROOT
ogprofiler export --run RUN_ROOT --family-fasta OG000000123
ogprofiler export --run RUN_ROOT --all-family-fasta
```

`export` 支持 `--config` / `--set`，默认读取 `run.yaml`。实际 OG 策略配置与完成产物不一致时，
先重新执行 `orthogroups`。源数据过期或 OG 产物损坏时明确报错，不隐式回退到 terminal 导出。

## 默认结果

| `results/` 文件 | 语义 |
| --- | --- |
| `families.tsv` | 最终 OG；保留历史 `family_id` 列名，但行对象是独立 ExportOrthogroup |
| `members.tsv` | 仅已分配蛋白的 OG membership，附 species_id / original_id |
| `terminal_families.tsv` | TF 前缀的终端 family 诊断；不是最终 OG |
| `terminal_members.tsv` | 全蛋白 terminal partition；包括真实 SSN isolates |
| `unassigned.tsv` | 未分配蛋白、原始身份、terminal_cluster_id 和原因；不强制生成 singleton OG |
| `statistics.tsv` | OG/terminal 数、input/assigned/unassigned 蛋白数、重复分配数 |
| `hierarchy.tsv` | 原始 hierarchy 结构及统计；isolate 的诊断行不写回 hierarchy store |
| `events.tsv` | 原有 network_event 诊断，独立于 v1_event |
| `export-manifest.json` | 当前导出清单、版本、输入/输出 hash、统计及每个 OG 的成员 hash |

OG 行字段：

```text
family_id component_id local_group_id source_cluster_id selection_type
v1_event processing_level n_genes n_species membership_hash
```

SSN_ISOLATE 的 source_cluster_id 为空。Python None 的 v1_event 写为 TSV 空值；
选择 trace 中的字符串 `None` 仍保留在 P3 Parquet，二者不混淆。

## 稳定全局 ID

`canonical-orthogroup-membership-rank-v1` 按每组排序后的 `(species_id, original_id)`
成员元组排序，再赋 `OG000000000` 等全局编号。

- 与 worker 完成顺序、组件遍历顺序、protein_id 重编号和 component_id 重编号无关。
- 不按 transient source_cluster_id 排名；成员集合改变仍可能改变后续编号。
- 终端诊断另用 TF 前缀；OG 与 terminal 不共享对象类型或 membership 含义。
- 组大小、物种数、membership_hash 从实际成员重新核查。
- assigned 和 unassigned 必须互斥并恰好覆盖 prepared 蛋白；重复或缺失为错误。

验证串行/双 worker 的全部交换表和全局 OG ID 一致，以及 protein/component 重编号后 ID 和
membership_hash 不变，完成 P3 的跨阶段全局 ID 验收。
12 份正常冻结 V1 fixture 的最终成员、等级和全局 ID 符合参考及排序契约；
剩余等成员 unary fixture 属于已冻结的重复分配异常，FAILED 状态阻止导出。
它不是合法 V2 split；现有 network annotation 对单子节点不支持，本轮没有调整该上游规则。

## FASTA

默认不生成每组 FASTA，只有 `--family-fasta` 或 `--all-family-fasta` 才生成。
主标识使用唯一的 prepared protein ID，原始身份保留为 header 属性，例如：

```text
>OGP2P000000000012 original_id=geneA protein_id=12 species_id=3
```

跨物种同名蛋白仍是不同条目，不再产生重复 FASTA 主标识。
选择集合、TSV member 数、FASTA 序列和 manifest 一致；未分配蛋白不进入 OG FASTA。
上次 manifest 管理、这次不再选择的 FASTA 移入 `results/previous-fasta/<unique-id>/`，
避免把历史文件误当成当前输出，同时保留历史内容。其他用户文件不纳入清理。

## 显式历史终端策略

```bash
ogprofiler export --run RUN_ROOT --strategy terminal
```

该策略写到 **`results/terminal-families/`**，不覆盖默认 OG 结果。
保留原四表和 `canonical-terminal-membership-rank-v1` 兼容编号（历史 OG 前缀只在这个独立命名空间中）。
新默认结果的终端诊断采用 TF 前缀。终端策略对比只是工程接口，本轮未启动科学准确度评价。

## 消费者验证、缓存和发布

新增 read-only `verified_orthogroup_inputs` 消费者门禁：总 manifest 必须 DONE，组件清单必须匹配
当前 index；组件 manifest hash、当前算法/schema/config、当前上游校验和、五份 Parquet
校验和及 Arrow schema 全部核对。忽略清单外历史组件目录，不凭 glob 的目录存在判断当前产物。

导出版本升级为 `orthogroup-and-terminal-diagnostic-export-v2`；旧 terminal export cache 失效。
新 cache 绑定 OG manifest、所有 OG 产物及上游输入、prepared FASTA、策略与 FASTA 参数；
验证 OG 当前有效性先于导出缓存复用。文件损坏会重写导出；上游损坏/改变须先恢复 OG 阶段。

交换表和 FASTA 采用 temporary + os.replace；重建时移除 completion manifest，最后原子发布新清单。
这不是整目录事务；中途写入失败可能留下部分已替换表，但没有可复用的完成标记。
当前导出仍沿用全局 canonical ranking 的 O(proteins + hierarchy nodes) 内存规模，
不保存每个内部节点的完整 descendant list；超大数据性能验收属于 P5，未在此宣称达标。

## Pairwise、phylogeny 与 benchmark

- pairwise 算法版本 `hierarchy-cross-child-v2`，仍从 hierarchy terminal partition 与
  **network_event** 的 cross-child 规则产生候选；OG grouping_dependency 明确为 `none`。
  修改/损坏 OG 产物不影响其独立缓存；每个 OG 的跨物种全配对不是默认 ortholog。
- `run ... --until-stage export` 仅在 `output.emit_pairwise_orthologs=true` 时额外调度
  `orthologs`；dry-run 同样展示该显式 opt-in。单独 `orthologs` 保持独立可调用。
- 可选 gene-tree refinement 继续消费 terminal 诊断。默认导出后的显式 family 参数使用 TF ID；
  历史 terminal 子目录保留其兼容编号。版本为 `selected-terminal-family-phylogeny-v2`，
  不把 OG source_cluster 的内部节点偷换成原 terminal refinement 对象，也不实现 OG refinement。
- benchmark 的最终 grouping 指标读取 OG 表；hierarchy/network-event 指标读取 terminal 诊断。
  pairwise 的 ID 映射使用完整 terminal membership，避免遗漏 unassigned 蛋白；记录额外输入 checksum。
  本轮只修复数据依赖接线，未变更评分公式/ground-truth 政策，也未评价准确度。

组件 GraphML 仍为独立 opt-in，来源和语义保持不变。

## 本地验收记录

P4 新增 **31 项集成测试**，并将原 13 份冻结 V1 落盘测试延伸到最终 OG 成员和全局 ID；
完整测试集 **256 passed**，14 条已有 SciPy covariance warning。
Mypy：76 个源文件通过；本轮修改文件 Ruff 及 `git diff --check` 通过。
还以 `PYTHONPATH=src` 调用真实 CLI 子进程，验证双 worker、overlap 配置向 export 传递及产物一致性。

本轮只做小型本地 fixture 验收，没有启动 Slurm、真实 hierarchy 回归或 Orthobench 评价。
P3 提交为 `4e97a63`；P4 为该提交后的独立工作区改动。
