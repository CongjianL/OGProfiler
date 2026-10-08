# P6：Orthobench 策略贡献评价

状态：**作业 `1410751` 评分与校验完成；OG 语义兼容通过，但准确度恢复目标未达到**。
验收日期：2026-10-02（Asia/Shanghai）。

## P6 结果与科学判断

Slurm COMPLETED、ExitCode `0:0`，耗时 `10:43:11`；
2026-10-01 20:31:02 开始，2026-10-02 07:14:13 结束（Asia/Shanghai）。
`p6-completion.json` 确认 evaluation_completed、benchmark_unchanged、
strategy_parity_passed 均为 true。五份评分均通过原始 ID 映射与官方 reader 校验。
两组独立 V1 parity 分别覆盖 30,166 / 30,165 个组件，失败组件为 0，
Parquet/TSV 与冻结输入校验通过。

以下均为 **百分数**，来自实际官方函数未舍入值，不是 best-group diagnostics：

| 对照 | Precision % | Recall % | F1 % | Coverage % | 组数 |
| --- | ---: | ---: | ---: | ---: | ---: |
| historical_frozen terminal | 0.058863 | 80.972900 | 0.117641 | 100 | 33,743 |
| historical_frozen v1_compatible | 0.058919 | 81.049326 | 0.117752 | 100 | 33,094 |
| repaired_mean terminal | 0.058506 | 80.436875 | 0.116928 | 100 | 33,927 |
| repaired_mean v1_compatible | 0.058822 | 80.870489 | 0.117558 | 100 | 33,283 |
| historical_full_v1 | 78.121549 | 45.501087 | 57.507516 | 100 | 58,218 |

每组 RefOG 原始 universe 覆盖率也是 100%，但这不代表成员分组正确。
组内 v1_compatible 的变化（**百分点**）：

- historical_frozen：P +0.000055524，R +0.076426167，F1 +0.000110967。
- repaired_mean：P +0.000315206，R +0.433613897，F1 +0.000629955。

提取策略有极小幅数值改善，但 precision/F1 仍处于崩溃水平，不能称为准确度修复。
70 行 confident best-group diagnostics：历史组 1 个 RefOG 的 best-F1 改善、69 个不变；
修复组 3 个改善、67 个不变，两组均无 best-F1 下降；这不是官方 pairwise F1 的替代。
完整历史 V1 对照与旧汇总 78.1/45.5/57.5 一致，但其源码版本、搜索和层级构建不同，
因此该巨大差距仍需分层定位，不能直接归因为当前 OG 代码未复现参考函数。

### 已定位的大组案例：hierarchy 拒绝候选后整体终止

两组 component 0 都只有一个 root 节点，child_count=0：

| 项目 | historical_frozen | repaired_mean |
| --- | ---: | ---: |
| root 蛋白数 / 物种数 | 69,701 / 12 | 69,642 / 12 |
| root terminal_reason | UNSTABLE | UNSTABLE |
| 候选数 / Leiden 调用数 | 10 / 30 | 10 / 30 |
| 有效候选数 | 0 | 0 |

实际 terminal 输出已经将该 root 整体作为一组。
v1_compatible 的 OG000000000 则为 source_cluster_id=0、selection_type=RESIDUAL_NONE、
processing_level=1，完整保留相同成员数；不是 OG 提取阶段新合并出的巨组。
由于 degree-zero root 的 n_genes > n_species，参考事件为 None，实际 V1 残余分支
同样选择整组；两组 parity 证实了这一行为。
巨组在至少 21 个 RefOG 的 fragment 示例中出现，是 12 个 RefOG 的 best group。
例如 repaired RefOG004 的 best group 有 69,630 个额外成员，RefOG025 有 69,629 个。
这给出明确过度合并案例；尚未单独计算该巨组占全局官方 FP 的比例。

关键不是“Leiden 没产生划分”，而是 **产生的全部候选被准入规则拒绝**：

- repaired gamma=0.01：148 个子群，min_child_size=2，max_child_fraction=0.105913；
  stability=0.834731 < 0.9，因此 UNSTABLE。
- repaired gamma=0.02 起的候选包含 singleton；配置 min_family_size=2，
  `minimum < min_child_size` 使整个划分被标为 GAMMA_LIMIT。
- repaired gamma=5.12：1,971 个子群，stability=0.905719 已过 0.9，
  但 minimum=1，仍因 GAMMA_LIMIT 拒绝。
- root 终止理由采用集合优先级：只要任何候选为 UNSTABLE，就汇总为 UNSTABLE，
  从而掩盖了同时存在的最小子群规则拒绝；不是全部候选均仅因稳定性失败。

对应代码为 `hierarchy/resolution.py` 的 `_candidate` 先检查 min_child_size 等约束，
再检查 stability；无有效候选时根直接终止。根输出为 family/OG 后导致大量混合成员。
降低 stability 阈值不是本次已经验证的修复，也不能据此跳过 singleton 准入与失败回退审计。

**下一步定位层：hierarchy 候选准入、singleton 子群处理、无有效候选时的回退，
再对照 V1 层级构建规则。保持 OG engine 和 scorer 固定，不同时改多层。**
2026-10-02 已完成该专项源码/fixture 审计，见 `hierarchy-acceptance-fallback-audit.md`；
发现准入政策差异、失败回退语义及搜索覆盖缺口，生产规则尚未调整。
映射/提取兼容检查已经通过；这不证明当前 hierarchy 与 V1 的生物学分裂等价。
本次没有调参、重提任务、修改生产算法或宣称准确度达标。

### 本地证据

小型结果保存在 `.provenance/slurm/20261001T123223Z_f69c14cdfe9b_c567bab0_49330/results/`：
五份 official stdout/metrics、五份 70 行 diagnostics、两份 paired delta、评分 report、
两组 parity 汇总、completion、component 0 的 nodes/candidates/metrics、
local-validation-summary.json 与 retrieved-files.sha256。
大型成员、预测和完整 hierarchy 保持远端。

## 提交记录

- RUN_ID：`20261001T123223Z_f69c14cdfe9b_c567bab0_49330`
- 基础 commit：`f69c14cdfe9bee9ed1aa251e89a511382b23baca`，dirty=1，精确 diff/源码快照已保存。
- 源码 SHA256：`c567bab0daca12541c3eb4a9201fb8706e594ea0e7c57b284ddb61252bc9bb41`
- 远端 run：`/home/mselab/licj/projects/running/ogprofiler-runs/20261001T123223Z_f69c14cdfe9b_c567bab0_49330`
- 本地：276 tests passed；远端小型 P5/P6 测试：20 passed；Ruff、76 文件 mypy、bash 语法检查通过。

本节在提交后补充；作业源码快照保存了本说明的提交前版本，执行代码未改动。

## 预先固定的设计

P5 数据为 13 个实际蛋白组，不是带 RefOG 真值的 Orthobench。
P6 改用此前 `v2_modified_rep1` 的 Open Orthobench：12 物种、251,378 蛋白，
70 个 RefOG、1,945 个原始真值成员；输入 digest
`380c85d9548c607f5df6daf656d74b21597fe215ef78d4f8ff296776bb1fdd07`。

| 对照组 | SSN/hierarchy | 两种策略 |
| --- | --- | --- |
| historical_frozen | 复制历史 forward SSN 与已有 hierarchy，冻结校验 | terminal / v1_compatible |
| repaired_mean | 复用相同 hits，当前 SSN 修复算法 + mean；按原 hierarchy 参数构建一次后冻结 | terminal / v1_compatible |
| historical_full_v1 | 冻结的完整首版本独立跑法 job 1409360 | 原始完整结果单列 |

第一组不重建 SSN/hierarchy；第二组不重新搜索、不调 hierarchy 参数。
第二组只按已经批准的修复算法重新生成 SSN，明确将 forward 改为 mean，
沿用历史 seed=42、robust/adaptive、max_depth=20、workers=56、subtree_workers=1。
每组内部策略比较共用同一个 hierarchy；跨组变化属于上游 SSN/hierarchy 总效应，
不混称为 OG 提取贡献。

使用 `fixed_hierarchy` 工具对两组分别执行实际 V1 参考差分与生产 OG 落盘校验，
全部通过后才评分。terminal 导出写入独立 `results/terminal-families` 命名空间。
没有同时改变 SSN、层级规则、OG 策略参数或 scorer，也不运行 parameter sweep。

历史完整 V1 使用 `OGProfiler1First/rep1_retry2` 的 `groups_rep1.tsv` 与
RUN_MANIFEST（版本 commit `729675fee0dabc02f63ae5c2cda7571748751deb`）。
它不是 P5 固定参考 SHA 所代表的同一源码版本，其 search、SSN 和 hierarchy
也不同；只作为完整历史结果参照，不作提取策略归因或偷偷替代固定参考。
历史序列/ID 适配已恢复 original ID，输入 digest 再校验，原结果保持只读。

## 固定评估口径与失败门槛

1. 官方 primary：原 `BENCHMARKS/benchmark.py`，默认 smart reader、默认
   `q_even=True`、逐 RefOG low-certainty exclusion。保留实际 CLI stdout/stderr。
   scorer SHA256：`81eb1e660c17819549b07eea8a54b4fb42a89180cafeb4569d92195d282f5e6f`。
2. 调用同一实际官方函数取未舍入 P/R/F1，并核对 CLI 一位小数结果；
   JSON 为 0–1 比例，不与历史表 0–100 百分数混用。
3. coverage 单独报告全输入已分配比例与原始 RefOG universe 覆盖率。
   未分配蛋白保持未分配，不补单蛋白 OG，也不只评分 1,945 个真值蛋白子集。
4. 原始 ID 须在完整官方输入内、无重复分配；正式计算前确认 prepared metadata
   与原始输入 ID 集合相同。官方 smart reader 必须与预测组集合完全一致，
   有任何 ID 丢失/误解析即失败。官方 stdout 的 ERROR 即失败，即使退出码为 0。
5. 附加 diagnostics：按 confident RefOG 的 best-group precision/recall/F1、
   碎片数量、额外基因、missing，以及过合并/拆分类别；**不是官方 pairwise 指标**。
   低确定性基因逐 RefOG 排除；保留少量额外成员/fragment ID 示例供复查。
6. 输出两组各 70 行 paired diagnostics delta 和各自官方 P/R/F1/覆盖率 signed delta。
   不只展示 F1，不隐藏 precision–recall 权衡，不作单次固定 seed 的显著性结论。

科学判断：只用每组内 `v1_compatible - terminal` 判断提取策略贡献；
若差分与映射检查通过但 accuracy 较低，先按层级/参考策略本身分析，
不把它称为未定位的代码错误，不同时修改多个层。
`evaluation_completed=true` 表示评分与校验完成，不等于 accuracy 达标。

## 执行与 provenance

工具：`benchmarks/og_extraction/orthobench.py`；脚本：`slurm/p6_orthobench.sh`。
通过 `dev/remote-test` 小型测试和 `dev/slurm-submit` 单作业执行，无 array。
资源沿用现有 Orthobench modified 正式脚本：56 CPU、250 GB、72 小时；
partition/account/qos 沿用本项目 cu/mselab/normal。与 P5 的 128 GB/24h 不同，
因为这次是 251,378 蛋白与 56 worker 的 Orthobench campaign，非重复 P5 数据。
不修改 shared environment 或历史结果，不覆盖已有科学目录。

独立 run 保存 benchmark.py、官方 Input/RefOG 文本快照、输入 SHA、原始路径、
历史 V1 groups/manifest 副本、dirty source snapshot、两组 Parquet/TSV、
两份 parity report、五份官方评分、paired delta、环境与阶段 time。
只有所有评分和 parity 通过且 benchmark 快照未变时才写 `p6-completion.json`。
失败时先查明确阶段日志，保留成功的固定 hierarchy，不盲目重提全套任务。
