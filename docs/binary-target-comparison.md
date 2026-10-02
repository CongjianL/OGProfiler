# V1小节点binary目标与V2 k-way的固定节点对照

Date：2026-10-02。来源job1410777；生产代码、配置及OG/scorer保持原样。

## 源码语义

V1 RunCommunityDetection：3..9999蛋白多物种节点要求exactly2；>=10000使用20..50
目标并接受>=2；多物种2蛋白直接拆成singleton；其余停止。
节点<1000时BipartiteGraphs从gamma=.5在[0,1]做社区数方向的中点搜索，最多1000步；
更大节点另有SpecificBoard括区步骤，其中扩张没有严格上界。Leiden10轮、无固定seed，
没有现行ARI=.9与最大子群比例=.95门槛。

因此，“V1小节点binary目标”不等于“现行V2仅增加child_count==2就精确复现V1”。

## 两个独立有限诊断

选择六个固定小节点：2694:0、2694:3、470:0、84:5、844:2、500:0（最多46蛋白）。
取回五组件约608KiB边分区，逐文件校验与原hierarchy输入manifest一致；诱导图保留
原global ID顺序、RBER目标和mean边权。冻结选中gamma的群数与ARI全部复现。

1. **共同候选bank**：10个coarse、精确端点、8个rescue、冻结local点；去重后每节点
   22点/66调用。比较同一批候选的k>=2准入与k==2目标，除此以外原门槛不变。
   这是固定bank，不假称它是生产adaptive访问日程。
2. **V1-style目标探测**：从.5在[0,1]有限中点探测社区数，找到第一个代表性binary
   就记录其原门槛结果；使用固定三seed及V2代表性选择，最多24点/72调用。
   六例实际1..7点。没有复制V1不固定seed的1000步完整搜索，也不是全V1 replay。

总计462调用；每个诊断独立记账，不叠入生产单节点24点预算或解释为生产fallback。
每调用10轮、robust seeds42/104771/209500、mean ARI>=.9、stop_size1、最大子群.95不变。
目标探测的[0,1]是显式诊断搜索差异；六例实际点均落在生产[.01,10]范围内。

## 结果

| 节点 | 蛋白 | 冻结V2 | 共同bank合格binary | V1-style首个binary | 原门槛 / 事件 |
| --- | ---: | --- | --- | --- | --- |
| 2694:0 | 15 | 4-way / III-3 | 无 | gamma .625 | **通过，ARI1；III-1** |
| 2694:3 | 11 | 11-way / III-3 | 无 | gamma .8125 | **通过，ARI1；I** |
| 470:0 | 46 | binary / III-1 | gamma .08 | gamma .125 | **通过，ARI1；III-1** |
| 84:5 | 25 | binary / III-1 | gamma .3805462768 | gamma .5 | **通过，ARI1；III-1** |
| 844:2 | 26 | 26-way / III-3 | 无 | gamma .6875 | **比例拒绝，25/26=.961538；III-1** |
| 500:0 | 45 | 4-way / III-3 | 无 | gamma .0234375 | **比例拒绝，44/45=.977778；III-1** |

2694:0在bank gamma .64确实有binary，但mean ARI=.285714而被拒绝；
很近的.625却稳定1，说明稳定性可行性不能被当作单调函数。
2694:3的bank完全没有代表性binary，却在.8125有可行候选。这支持改进目标搜索覆盖，
不支持通过放宽稳定性恢复它。

对于2694:3的11个单拷贝成员，真实binary候选把1+10个物种分为不相交集合，
其V1事件确为I；这比此前仅按species bitmap的假设分组更强，但仍只是**冻结节点**证据。
新binary根及其后代不一定保留同一个node3，所以没有宣称RefOG001整体已恢复。

## 设计结论

1. binary目标有明确正例：现有SSN与原门槛下存在候选，能恢复该节点I资格。
2. 单纯增加k==2条件不够：相同bank漏掉正例，且可能将已resolved树变成UNRESOLVED。
3. binary不保证I：三个通过门槛的binary仍为III-1，必须保留这些反例。
4. 另两个binary即使ARI1也违反.95比例；本次不放宽条件、不强收、不升级为terminal。
5. 当前V1与V2还存在搜索域、seed、可靠性门槛、size2特判与失败停止差异；
   本次隔离了一步目标/点覆盖问题，没有将所有策略差异混成“binary一定更准”。

据此形成 **ADR0004 Proposed**。原生产k-way政策仍生效；需明确完整目标日程、共享配额、
失败状态及是否存在显式soft回退后，再进行小组件整树及后续Slurm验证。
不修改停止尺寸、ARI阈值、OG/scorer；不把本次单节点事件结果当作全pipeline评分。

## 复查

工具：`benchmarks/og_extraction/binary_target_audit.py`；三个单元测试覆盖bank去重/预算、
species disjoint binary I、overlap binary III-1及单拷贝k-way III-3。
与既有RefOG/hierarchy契约一起共20项针对性本地测试。

原始候选曲线与每个拒绝原因：
`/Users/licongjian/Desktop/PycharmProjects/OGProfiler/.provenance/slurm/20261002T060705Z_cb0c0ba98049_1cb4fcd7_14702/results/binary-target-comparison.json`

该报告保存实际冻结trace、全部候选指标、预算和审计脚本SHA256。无需远端推理或新Slurm
作业；本次工具、测试、设计草案未提交Git。
