# ADR 0003：H4/H5 迭代预算10的真实数据验收

Date：2026-10-02（Asia/Shanghai）。

## 来源与执行状态

- Job：`1410777`，Slurm `COMPLETED / 0:0`，总耗时7:01:53。
- 开始14:05:42，结束21:07:35；集群时钟与本地提交时钟有小幅差异，保留各自原始记录。
- RUN_ID：`20261002T060705Z_cb0c0ba98049_1cb4fcd7_14702`。
- Commit：`cb0c0ba98049ed9a485e17ce9f420fd7c7ffac70`，dirty=0。
- Source SHA256：`1cb4fcd745e670e85b7389bbde979f0c39c8f400e4c9a5923918aa970fe87d52`。
- 不可变源码快照、原资源56 CPU/250G/72小时；未重新搜索hits或构建SSN。

`iteration-only-control.json`确认与job1410775固定输入哈希相同，配置仅改变
`hierarchy.leiden_iterations: 2 -> 10`。停止尺寸1、mean ARI阈值.9、seed42及
robust三seed、每节点24候选、RBER、mean边权全部保持原样。

## H4：完整component0，而非只检查旧818节点

`h4-acceptance.json`全部9项验收为true：resolved、成员唯一完整覆盖、叶尺寸、
候选预算、robust调用计数、选中候选政策、两个manifest和串并行一致性。

| 指标 | job1410775：2轮 | job1410777：10轮 |
| --- | ---: | ---: |
| 蛋白数 | 69,642 | 69,642 |
| unresolved节点/蛋白 | 1 / 818 | **0 / 0** |
| 根gamma | .11313708498984759 | .09513656920021767 |
| 根子群数 | 444 | 470 |
| SPLIT节点 | 28,176 | 29,418 |
| 政策TERMINAL节点 | 53,266 | 53,751 |
| Leiden调用 | 865,206 | 896,919 |
| 主H4墙钟时间 | 4:06.96 | 12:50.00 |

根和上游拓扑已改变，旧cluster_id和818边界不是本次树中稳定标识。
此次结论是完整component0已resolved，不是将旧818节点静默标为停止叶。
10轮的调用数略增、耗时约3.12倍；更充分优化有实际计算成本。

独立两worker回放耗时10:26.69，nodes/members/candidates三张表精确相等。
两份结果分别有有效manifest；这不是只对根候选或成员数的弱一致性检查。
取回的component0四项产物哈希均与manifest匹配；本地独立复核83,169节点、
298,973候选，每节点实际最多14点（上限24），选中最小ARI .9000082112697011。
主H4 GNU time peak RSS约772MiB；两worker日志的RSS不作为所有worker内存总和。

## H5：全层级完成，固定hierarchy下OG parity通过

- scheduler：10,712个非singleton组件已resolved，其中component0复用；
  新执行10,711、failed=0、unresolved=0；另19,453个singleton组件。
- 总hierarchy Leiden调用3,186,366，候选预算通过；其余层级耗时1:00:49。
- 独立实际V1提取对照覆盖30,165组件，failed_components=0、first_divergence=null；
  产物验收与冻结输入校验通过。parity阶段耗时3:47:30。
- 251,378蛋白全部分配一次，duplicates=0、unassigned=0。
- V1-compatible OG为95,895组，其中69,967个singleton；最大OG600，>=1000组为0。
  已消除旧69,642蛋白粗根直接输出的大组病态；不代表所有合并/拆分均正确。
- 固定SSN与官方benchmark哈希复核未变；正式smart reader精确识别预测成员。

## 官方Orthobench结果

以下为百分数；coverage与RefOG raw coverage全部100%。terminal仅为结构叶诊断，
V1-compatible才是当前正式OG提取策略。

| 结果 | Precision | Recall | F1 |
| --- | ---: | ---: | ---: |
| 旧P6 repaired-mean / V1-compatible | 0.0588% | 80.8705% | 0.1176% |
| 新hierarchy / terminal诊断 | 62.4810% | 1.0759% | 2.1153% |
| **新hierarchy / V1-compatible** | **91.7440%** | **25.5911%** | **40.0192%** |
| 历史完整V1 | 78.1215% | 45.5011% | 57.5075% |

同策略相较旧P6，F1提高39.9016个百分点；相较历史完整V1仍低17.4883个百分点，
主要表现为更高precision与更低recall。全蛋白覆盖不等于正确成对关系召回。

控制边界：2→10与job1410775的比较支持hierarchy执行修复，但该旧作业未完成H5。
评分的旧P6对照还使用旧准入/搜索政策，因此不能把39.90个百分点全部归因于迭代预算。
历史完整V1是独立全pipeline基线，不是仅切换hierarchy的一因素对照。

## 剩余问题与下一步

70个RefOG的best-group诊断：新策略SPLIT37、MERGE_AND_SPLIT23、EXACT8、OVERMERGE2。
这是独立诊断分类，不是官方pairwise分数；低recall与分散关系相符，但尚未定位首次分散阶段。
例如RefOG001的15成员分成15个预测组，best recall=1/15，成员并未丢失。

**ADR 0003 的完整resolved与工程契约验收通过；准确度达到历史V1的目标尚未完成。**
无需为本作业设计搜索失败回退，因为当前unresolved为0。
下一步按低recall RefOG审计首次分散位置：SSN连通分量 → hierarchy候选/拓扑 →
V1事件及OG选择。固定hierarchy parity已经通过，先定位上游差异，不直接修改OG/scorer，
不从结构singleton比例推导停止尺寸或稳定性阈值应降低。
binary与k-way、V1大/小组件分裂目标须另做明确设计与控制实验。

## 证据位置

远端原始结果保留在该RUN_ID目录。本地仅取回component0约数MB的诊断Parquet、
验收/配置/资源/评分与RefOG摘要：

`/Users/licongjian/Desktop/PycharmProjects/OGProfiler/.provenance/slurm/20261002T060705Z_cb0c0ba98049_1cb4fcd7_14702/results/`

其中 `h5-report-compact.json`为只移除逐根列表的报告副本，不重新计算评分。
源码与已提交快照保持原样，本次没有重提交Slurm或修改科学参数。
