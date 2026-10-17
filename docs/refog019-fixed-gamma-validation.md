# RefOG019 gamma .953125 固定点验证

2026-10-04。固定来源job1410868 mean SSN、component3976根（12蛋白）；
V1实际gamma/membership来自job1411083原始V1调用trace。

## 方案与边界

- 唯一待测gamma **.953125**，仅诊断调用的搜索边界设为同一点；原run.yaml不改。
- 通过生产 `search_resolution` 调用，原soft topology、代表性membership选择、
  robust3 seed、每调用10轮、seed42、稳定性.9、比例.95、停止尺寸1及OG/scorer保持不变。
- seed按生产规则为42、104771、209500。观测代理只记录实际分区，直接委托生产run_leiden，
  不额外运行另一套seed，不用V1seed-free结果替代代表分区。
- 仅1个unique gamma、3次Leiden调用；非生产日程新增点，不向原24点缓存追加候选。
- 保存各seed的成员/子群大小/quality、两两ARI、代表分区、原始violations、
  与V1实际membership的ARI及事件。冻结表、配置及图输入hash前后核对。
- 同时核对V1run和V2冻结输入checksum完全相同，igraph/leidenalg版本一致。
- 如果点通过binary+原门槛，定位为本固定点覆盖遗漏；如果UNSTABLE，
  定位为该点稳定性冲突，未采样并不代表补点能修复。
  结论不外推未采样区间，也不代表完整后代验收。
- 输出在独立不可变Slurm run；1CPU/8G/1h，与前一小诊断资源相同，无array。

工具：`benchmarks/og_extraction/refog019_fixed_gamma.py`；
脚本：`slurm/refog019_fixed_gamma.sh`。
本地9项相关fixture测试通过（含固定点成功与UNSTABLE反例、配置不变、三seed及ARI计算），
ruff/bash语法及diff检查通过。集群fixture验证后提交。

## 执行身份

- Slurm job **1411086**。
- RUN_ID `20261004T044604Z_27a2533e1334_7496e3b5_90608`。
- HEAD `27a2533e13340b9192e39251dc18db3af028fc7c`，dirty=1；不是clean HEAD结果。
- 精确源码snapshot SHA256 `7496e3b55e3dfaad4cafe0a2114d9abf5ea8e374dc03d95463fd0d561793031c`。
- 本地及远端fixture测试均9 passed。

## 实际结果：稳定性冲突，不是单纯覆盖遗漏

job1411086 **COMPLETED / 0:0**，Slurm elapsed **4秒**。
报告：`.provenance/slurm/20261004T044604Z_27a2533e1334_7496e3b5_90608/results/report.json`。
输入/冻结产物checksum前后不变，mean SSN、版本和原准入均核对通过。
原始V1输入checksum与当前根的冻结输入完全相同。

| seed | 子群尺寸 | 隔离的singleton原始ID | 与V1实际分区的ARI |
|---|---|---|---:|
| 42 | 11+1 | ENSDARP00000146398 | -.09090909090909119 |
| 104771 | 11+1 | ENSDARP00000066396 | 1 |
| 209500 | 11+1 | ENSDARP00000066396 | 1 |

两两ARI：

- 42 vs104771：**-.09090909090909119**。
- 42 vs209500：**-.09090909090909119**。
- 104771 vs209500：**1**。
- 原生产mean ARI：**.27272727272727254**，低于.9。

不是某个seed未产生binary：三个都恰好2群，但隔离的是不同斑马鱼蛋白，
成员边界不同。不能仅比较child_count或11+1尺寸就认为重复划分一致。

生产代表分区已选两个相同seed支持的分区，且**与V1实际分区ARI=1**。
代表分区最大子群比例11/12=.916667，通过.95；结构有效、tiny比例有效。
唯一原始拒绝原因是 **UNSTABLE**，binary_eligible=false、kway_eligible=false。
原事件代数为III-1，根本身并不直接获得I选择资格；此处未运行其完整后代或OG。

## 结论与范围限制

1. 该gamma确实未进入原24点日程，但它不是原V2门槛下的漏掉的合格binary。
2. 即使单独补入准确的V1 gamma，原稳定性门槛仍会拒绝；
   诊断分类为 `STABILITY_CONFLICT_AT_FIXED_POINT`。
3. V1实际分区并非V2代表选择算法找不到；它已成为代表，只是整组三seed的
   mean ARI没有通过。把代表机制换成多数投票也不改变这一既定mean ARI准入条件。
4. 此结果仅排除“.953125补点即可通过”的解释；不证明整个gamma区间无稳定binary。
5. 停止尺寸1、ARI.9、比例.95、seed规则、优化轮数、搜索日程和OG/scorer保持不变。
   没有将本点写回候选缓存、没有启用额外预算或发布失败节点OG。
6. 若继续，应先审计小节点中单个重复物种拷贝交换对ARI的影响及V1单次随机分区的
   可重复性；任何稳定性定义/拓扑政策调整仍需独立证据和显式ADR修订，
   不是本固定点验证自动批准的修改。

本地及远端相关fixture **9 passed**，ruff/bash/diff检查通过。代码和文档尚未提交。
