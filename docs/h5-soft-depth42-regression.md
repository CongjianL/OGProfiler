# H5 soft+depth42 全组件实验验收

日期：2026-10-04（Asia/Shanghai）。job1410868，源码27a2533，dirty=0。
RUN_ID：20261003T120507Z_27a2533e1334_72b113e6_93925。
Slurm COMPLETED/0:0；2026-10-03 20:03:46至2026-10-04 01:59:55，5:56:09。
资源56CPU/250G/72小时；未重新计算hits/SSN，OG/scorer算法保持原样。

## 工程验收与复用

- job1410845的component0已验收并复制，scheduler DONE/attempts=0，未重算；H5后hash未变。
- 新调度10,711个非singleton组件，跳过component0；failed=0、unresolved=0。
  总10,712个非singleton组件，另19,453个singleton组件。
- 完整hierarchy resolved：368,278节点，170,911分裂、197,367政策叶；最大深度33。
  每节点预算通过，总调用6,700,458（包含复用component0的历史调用，不等于本job新增调用）。
- 6,384次FALLBACK_KWAY、300次refinement截断；部分根为ONE_SPECIES政策叶，
  summary的root_split=false不是未解析或验收失败。
- resolved gate通过，10,712份组件manifest冻结，输入与benchmark哈希复核未改变。
- 固定新hierarchy的实际冻结V1 OG parity覆盖30,165组件：失败0、first_divergence=null，
  产物验收通过。证明同hierarchy下提取策略一致，不证明与V1完整pipeline拓扑一致。
- 251,378蛋白完整分配、coverage与RefOG raw coverage均100%；67,916个OG，
  singleton OG41,276，最大OG600，>=1000的大组为0。无69,642蛋白粗根污染重现。
- hierarchy-all墙时55:22.74；parity 3:45:57；官方评分1:20.30。

## 官方Orthobench结果（百分数）

|条件|Precision|Recall|F1|
|---|---:|---:|---:|
|job1410777：kway、10轮、depth20、V1-compatible|91.7440%|25.5911%|40.0192%|
|job1410868：soft、10轮、depth42、V1-compatible|92.6961%|33.8220%|49.5608%|
|历史完整V1（独立全pipeline）|78.1215%|45.5011%|57.5075%|

相对同meanSSN的job1410777，P +0.9521、R +8.2309、F1 +9.5416个百分点。
这是soft topology/search日程及显式深度预算的**联合修订**对照，不将全部改善归因于binary
或深度单因素。job1410810仅改变topology但留下DEPTH_LIMIT，未完成评分。
相对历史完整V1，新结果F1仍低7.9467个百分点、R低11.6791个百分点；准确度目标尚未达到。
历史V1是独立pipeline，不是固定SSN的一因素对照。

terminal仅作结构叶诊断：P62.5849%、R1.0901%、F1 2.1429%；不等同正式OG策略。
官方scorer SHA256仍81eb1e660c17819549b07eea8a54b4fb42a89180cafeb4569d92195d282f5e6f，
mapping/official reader检查通过，不调整scorer、ID映射或不确定成员口径。

## 碎片化与污染诊断（非官方pairwise指标）

70个RefOG的best-group confident recall：25个提高、45个持平、0个下降。
RefOG-预测组交叠片段总数809→488（是逐RefOG交叠计数，非全局独立组数）。
例如RefOG001：best recall 1/15→13/15、片段15→3；RefOG070：1/13→12/13。
但RefOG012保持15/46、片段6→6，仍是不能被binary优先完全修复的反例。

分类仍SPLIT37、MERGE_AND_SPLIT22、EXACT8、OVERMERGE3；旧kway分别37/23/8/2。
OVERMERGE由2增至3，不能宣称污染全面改善；最大OG600和>=1000为0也不是无污染证明。
缺失confident成员总数为0；coverage恢复不代表pairwise recall恢复。

## 科学结论与下一步

全组件工程契约通过；在同输入、同提取策略/官方scorer下，联合政策的聚合F1/R实测改善，
但仍低于历史完整V1，不能发布“准确度已完全修复”的结论。脚本的
accuracy_improvement_claimed=false为预设保守标记；上述增量来自读取冻结官方分数，
不是修改报告标记或重新评分。

继续定位未改善及低召回RefOG首次在SSN、hierarchy或事件选择阶段拆散的位置，重点保留
RefOG012等负例，并检查新增OVERMERGE。暂不降低停止尺寸、ARI或比例门槛，不改OG/scorer。
本轮未提交新作业，未迁移全局默认topology/max_depth。

## 本地证据

仅取回小型摘要、日志与RefOG表，原始结果/SSN继续留远端。
`.provenance/slurm/20261003T120507Z_27a2533e1334_72b113e6_93925/results/`中包括：
- h5-metrics/report.json：官方报告（内置对照为旧P6 repaired-mean，不能当作1410777对照）。
- baseline1410777/report.json：另取回的同meanSSN kway报告。
- soft-vs-kway-deltas.json：本地直接相减的增量及逐RefOG诊断。
- h5-parity/report.json、scheduler-manifest.json、component0-reuse-check.json、h5-hierarchy-freeze.json。

后续只读定位见 [首次偏离审计](h5-refog-first-deviation-audit.md)：召回分解、候选准入冲突及全部新增污染组路径。
