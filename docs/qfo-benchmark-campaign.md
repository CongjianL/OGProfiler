# QFO比较分析：soft42生产默认（2026-10-06）

用户明确要求继续使用 /home/mselab/licj/project_data/OGProfiler/QFO 做benchmark。
本地源码权威，用户目录只读；所有全数据处理及计算通过Slurm。

## 当前执行范围：bacteria 独立数据集

2026-10-06 用户决定停止 all，改用
`/home/mselab/licj/project_data/OGProfiler/QFO/bacteria` 的 23 个物种。
本节取代下文历史 all 计划；旧快照与结果保留，仅作归档。

job **1411823** 经 sacct 核验为 **COMPLETED / 0:0**（2分18秒），
审计、完整 prepare、3×20 search smoke 与两个方法环境检查均退出 0。
`subset_checks_passed=false` 表示目录差异，不是作业失败；
该作业尚未运行完整 V2/V1/OF3，也未产出正式 QFO 分数。

新入口：`slurm/qfo_bacteria_preflight.sh BACTERIA_DIR OF_ENV V1_ENV [PREVIOUS_JOB]`。
资源保持 **4 CPU / 16 GB / 2小时**。通过标准提交器冻结源码和 provenance，
记录 supersedes_job=1411823 与用户指定的范围变更。
输入审计使用 `--collection bacteria`，只读取所给细菌目录，
不打开 all/eukaryota，也不做跨目录子集判断；单独生成内容摘要与 ID 映射。
输出 `preflight/input/bacteria`、`prepared-bacteria`，报告采用
`collection`、`n_species`、`n_proteins`、`dataset_sha256`，
`subset_audit_completed=false`、`subset_checks_passed=null`（未进行此项审计）。
使用生产 soft42 默认完成全细菌 prepare，并以确定性 3×20 子样本验证 search。
预检完成后再依据实际规模准备独立完整 V2/V1/OF3 运行；
所有比较限于相同冻结细菌输入，保留原有官方 QFO 评分门槛。

本地新增回归测试：细菌目录独立存在，跨目录比较函数若被调用即失败；
确认输入只快照到 bacteria，报告不含 all 专属统计。

## 历史 all 计划与已检查事实（已被上节取代）

- all 78 FASTA，约558MB；bacteria 23，约39MB；eukaryota 48，约491MB。
- all是主输入，另外两组独立审计与all的字节、ID及序列差异；其余文件不臆断域标签。
- 头部为sp|ACCESSION|或tr|ACCESSION|；数据release与canonical/isoform范围未有元数据佐证。
- 目录内未发现非FASTA的评分/参考文件，既有benchmark未发现B2结果。
- 现有ogprofiler、ogprofiler-v1-first、orthofinder-3.1.5环境已检查存在。
- 原项目QFO_READINESS约束仍有效：官方QFO评分须有验证的orthology表示及评分包，
  不把OG组内全对展开后冒充ortholog，不把network_event当phylo_event。

## 分阶段执行

1. 本轮预检：4CPU/16GB/2小时、单作业，不是全量方法运行。为输入建立内容快照，
   逐文件SHA256、蛋白/残基数量、ID/UniProt accession唯一性、域目录差异分类。不把all以外序列并入主输入。
   同时保留original ID与accession双射，不裁剪isoform、不改非法字符为X。
2. 用当前soft42默认做78物种完整prepare输入验证，记录timing；仅3物种各20条
   确定性序列做DIAMOND search smoke。不把该子样本作为正式分数。
3. 预检成功后根据蛋白总量、资源和格式检查准备独立完整V2/V1/OF3运行，
   科学默认固定。既有Orthobench资源不是QFO资源充足性的证明；全运行提交前
   明确资源/时限并检查环境与源版本。此预检不隐式启动后续作业。
4. 全量比较：OG coverage、组大小/单例/巨大组、OF3参考一致性（pairwise及B-cubed
   如适用，明确是方法一致性而非生物学真值）、split/merge诊断、域分层和资源描述。
   OF3 unassigned遵循assigned主评估+未分配单独报告；不把未分配当负例。
5. 正式QFO生物学指标在匹配参考版本/评分包和验证orthology导出后独立开展，
   与OG划分一致性报告分开，保留当前READINESS状态，不发布伪官方分数。

## 可追溯性与文件

入口：slurm/qfo_preflight.sh QFO_ROOT OF_ENV V1_ENV。基础设施值来自项目配置/脚本参数。
用dev/slurm-submit创建独立源码快照，记录commit、dirty、source hash、job/run及参数。
所有步骤通过OGProfiler2_benchmark/workflows/run_timed.sh保存命令、时间、stderr、
退出状态及资源；原始输入和已存在结果不覆盖。presets不是环境安装包，PYTHONPATH
固定快照src，避免旧安装包覆盖新默认。环境声明与软件版本记录于run。

输出包括 preflight/qfo-input-audit.json、input-manifest.tsv、input-sha256.json、
id-map.tsv、subset-comparison.json、input/all快照、soft42-default.yaml、prepared-all与smoke日志。
qfo-preflight-completion.json仅表示预检完成；不代表正式benchmark已完成。
后续将紧凑报告接入OGProfiler2_benchmark/05_metrics/qfo，原始大结果留远端。

## job1411817失败与修正

sacct：FAILED/1:0，00:00:53。audit步骤51.32秒，峰值321464KB，退出1；
不是OOM、TIMEOUT或节点故障。仅输入审计运行，prepare/search/方法全量都未开始。
失败于23个bacteria同名文件的字节hash均不同，脚本把目录误设为严格子集。
轻量样本还发现eukaryota人类文件记录列表与all不同；不将该目录当作同一版本。

修订：all主输入仍严格冻结、ID/accession唯一和格式检查仍执行。其它目录逐文件
区分BYTE_IDENTICAL、FORMAT_OR_DESCRIPTION_ONLY、ID_DIFFERENCE、
SEQUENCE_SET_DIFFERENCE及NOT_IN_ALL；记录新增/缺失accession、序列变化、ID变化
数量和有限例子。不裁剪差异、不替换主数据、不把集合不一致标成通过。
input_ready只表示all可用；subset_checks_passed与subset_audit_completed独立。
域分层仅使用目录物种列表对all序列归类，不将独立子目录序列混入主评估。
后续若使用子目录做独立运行，需独立快照与dataset digest；不得声称与all严格同输入。

run_timed原报告git_commit=UNKNOWN（快照无.git、无旧.dev-deploy-meta），
修订为优先读取标准wrapper注入DEV_GIT_COMMIT；原作业根provenance完整，不回写历史记录。
重试保持4CPU/16GB/2h，无资源增加；通过新源码、新RUN_ID执行，并记录retry_of_job=1411817。
原失败目录、冻结输入和id-map保留。本地测试覆盖格式/顺序变化、真实序列差异、
ID变更、主输入不被混入及快照提交号传递。

真实细菌样本UP000000425_122586本地新逻辑验证：两目录各2001蛋白，accession集合
与序列完全一致，但3个original ID变化（P0A0Z2、Q9JXS2、Q9K0J8），分类ID_DIFFERENCE。
此证据仅适用于该样本，其余文件由重试审计逐一核验；不外推所有细菌差异性质。
样本hash与诊断归档 docs/qfo-preflight-1411817-sample-difference.json。
修订后完整测试423 passed / 8 skipped（5项既有warning）。

## bacteria 预检结果与方法验证（2026-10-07）

job1411846：COMPLETED / 0:0，27秒，源05a5239，所有6项计时步骤退出0。
独立细菌输入23物种、82507蛋白、26792356残基；SHA256：
`703718a40b464a4193a79e5070777a63d768a39ccaceda2f4c7cfd4885c55fe4`。
完整prepare 3.16秒、峰值151648KB；小样本search 8.26秒。
这些数值仅用于预检，未表示完整方法计算所需时间/内存。

方法入口：`slurm/qfo_bacteria_smoke.sh`（4CPU/16GB/2小时）与
`slurm/qfo_bacteria_full.sh`（32CPU/128GB/72小时）。完整资源档沿用历史
V1首次版本基准资源档，是新完整对照的显式请求，不是提高预检资源；
单作业按V2/OF3/V1顺序运行，无数组、无参数搜索，三种方法不同时争用资源。
参数均为 `PREFLIGHT_RUN EXPECTED_DIGEST OF_ENV V1_ENV`。
方法smoke验证前不提交full；smoke沿用冻结3×20，绝非正式分数。

`benchmarks.qfo.campaign`校验来源完成状态、bacteria标签、冻结摘要，
完整模式核对文件hash与蛋白/物种数量，独立复制同一组输入。
历史V1通用适配器拒绝UniProt竖线ID，因此专门生成
`OGPV1_species|QFO_species_record`，用双射表恢复完整原始ID；
不改变序列、物种成员或选择蛋白。回归测试涵盖原始ID恢复、
V1现有converter契约、被修改输入拒绝及all来源拒绝。
V1使用原首次版本SHA256 65ec43d269b410956fadc4f8215a3e0a0772eac8bcc1060ab011b090823cca41，
环境前缀由已配置环境查询；沿用已验证的DIAMOND makedb兼容启动器。
其CLI没有随机种子参数，明确记录unavailable，不把计时元数据seed42当成V1已控种子。
V1参数固定为历史B1配置（evalue1e-5、lrb、NBS、rber、gamma1、so0）。
V2加载冻结生产soft42配置，仅覆盖运行线程/worker数量，不调科学参数。
OF3使用diamond，`-og`停止于OG划分，不启动树推断；这是OG对照而非正式orthology评价。
所有方法及转换均经过run_timed，环境显式导出、源码快照和独立结果目录保留。
OF3原始assigned分组路径单独保存；标准converter添加的unassigned单例仅用于
完整覆盖验证，不用作OF3 assigned pairwise参考标签。
`qfo-methods-completion.json`只记录方法执行，comparison_evaluated=false；
实际一致性评分、split/merge与资源汇总留待完整输出检查后进行。

## job1411854小样本失败诊断

2026-10-07核验：FAILED / 2:0，25秒。V2 prepare/search/edges/components均成功，
60蛋白全部为单例，scheduler完成0个非单例组件且failed=0。
network annotation误将“无非单例层级输出”视为缺失输入；OF3/V1尚未启动。
这是全单例边界代码错误，不是资源问题。保留样本与科学参数，资源不增加。
修正仅接受scheduler RESOLVED、零失败/未解决、单例数与index/proteins完全匹配的
全单例数据集；生成零组件/零事件注释manifest。非单例遗漏或输入不一致仍报错。
新增正例及缺scheduler、非单例index、缺单例表、缺蛋白负例回归。
重试用新源码快照与RUN_ID，第五参数记录retry_of_job=1411854；
完整方法运行仍需等待三方法smoke验证成功。

## job1411870：原版 V1 的全单例边界与验证集调整

sacct FAILED / 1:0，49秒；V2端到端23.19秒退出0，OF3 12.11秒及转换退出0。
V1 7.71秒退出1，ExtractOG访问空层级图的Event属性时KeyError。
ssn.gml为60个节点/0条边，hm/hmm为空；所有物种对均无非自身有效hits，
OF3也将60条全部列为unassigned。未产生V1分组，三方法完成标记未生成。
这不是资源失败，不能把该作业说成三方法验证成功。
保留历史V1源hash，不添加算法修补或把单例文件冒充V1成功输出。

新增 `slurm/qfo_bacteria_smoke_proteomes.sh`：相同冻结输入的前三物种完整蛋白组，
模式smoke-proteomes，独立dataset digest并记录完整23物种digest；
4CPU/16GB/2h不变，retry_of_job=1411870。所有选中FASTA按字节复制，
不合成序列、不按参考标签挑选同源家族。这是明确扩大验证集范围，
不是调科学参数；已有3×20全单例结果继续保留。
成功后还需检查有效跨物种连接与三方法分组/ID完整性，
再提交23物种完整对照。回归覆盖前三物种选择、完整序列复制及独立digest。

## job1411874验证通过，启动完整细菌对照

2026-10-07：COMPLETED / 0:0，22分44秒；3物种、4449蛋白。
V2完整pipeline退出0，19分40.68秒；OF3退出0，2分02.29秒；
原版V1退出0，55.06秒；两项分组转换也退出0。
V2有803个非单例组件全部成功，failed=unresolved=0；V1 SSN有2385条边。
输出审计：V2 3004组、V1 2993组、OF3含补充单例的覆盖分区2899组；
三方法均4449行、每条输入蛋白恰好一次，无未知ID，无遗漏/重复。
V1 ID映射双射及原始UniProt ID恢复通过。
OF3原生assigned蛋白2382、unassigned2067，完整覆盖分区不冒充assigned参考。
这些是格式、覆盖和执行验证，并非方法精确率/召回率或优劣结论。

完整对照以预检1411846冻结的23物种82507蛋白为唯一输入，
SHA256 703718a40b464a4193a79e5070777a63d768a39ccaceda2f4c7cfd4885c55fe4。
使用已准备的full入口，32CPU/128GB/72小时；单作业顺序V2/OF3/V1，
完整数据重做search，不复用3物种样本结果；source/input/results均独立。
第五参数1411874记录为validated_smoke_job，不误记为失败重试。
完整作业结束后先核验所有输出，再执行assigned主评估、独立unassigned报告、
pairwise/B-cubed一致性、split/merge及资源汇总，仍不发布官方QFO分数。
