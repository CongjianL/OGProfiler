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

## 按用户要求：三方法并行执行

2026-10-07用户明确要求V2/OF3/历史V1并行。原顺序job1411889在V2
DIAMOND阶段（prepare已完成、搜索未完成、OF3/V1未启动）被精准取消，
原snapshot/provenance/输入/不完整日志全部保留，不改写运行中的源码。
新增三个入口：qfo_bacteria_v2.sh、qfo_bacteria_of3.sh、qfo_bacteria_v1.sh；
每个32CPU/128GB/72h，总并发申请96CPU/384GB，由Slurm实际调度。
三作业无相互依赖，共享同一冻结数据来源和source commit，但各自有独立
RUN_ID、源码快照、输入副本、结果与timing，避免同目录争用。
科学参数、搜索线程、worker数量不变；未完成的搜索重新执行，不当作有效缓存。
原顺序入口保留用于历史重放，此轮不再使用。选择器QFO_METHOD默认all以兼容
原smoke，独立入口仅启用一个方法；未知方法执行前报错。
第5参数validated_smoke_job=1411874，第6参数supersedes_sequential_job=1411889。
方法路由桩测试覆盖v2/of3/v1/all及错误值，确认没有误触发其它方法。
三方法输出核验结束后再统一比较，不用任一方法完成状态代替全campaign成功。
并行资源时间不相加成campaign墙钟时长；报告需分别保留各方法时间/内存，
另报告三个作业最早开始至最晚完成的campaign跨度。

## 三个完整方法完成，进入独立评估

2026-10-07核验1411893/1411894/1411895均COMPLETED / 0:0；
完整方法计时分别V2 3:19:20、OF3 43:43.31、V1 30:12.37。
这是方法运行时间，不是评分完成或优劣证据；OF3只运行到OG，不含树推断。
新增独立evaluate入口：4CPU/16GB/2h，经run_timed保存评分资源与版本。
检查三方法full完成标记、冻结输入逐文件一致、覆盖唯一、V1双射；
OF3原生assigned主参考与完整覆盖分区分开，unassigned独立汇总。
复用已有线性列联表compare计算pairwise/B-cubed、split/merge与top30诊断，
不枚举蛋白对、不把OF3视作生物学真值、不输出官方QFO指标。
组件诊断明确使用V2 retained graph，包括V1相对该图的诊断，不误称V1图的上界。
输出summary.json、metrics.tsv及产物SHA256；大结果树仍留远端。
回归测试验证assigned限定评分、unassigned混合报告、双方法已知答案和
被修改输入、缺失/重复预测、未知/重复OF3 ID的拒绝。

## job1411946评分完成

COMPLETED / 0:0，14秒。已取回紧凑summary.json与metrics.tsv，归档于
OGProfiler2_benchmark/05_metrics/qfo/bacteria_soft42_1411946，REPORT.md解释范围与局限。
三方法完整唯一覆盖82507蛋白；assigned主评分62466蛋白/9828 OG。
V2 pairP/R/F1=76.5618/54.2763/63.5211%，V1=46.2698/70.5801/55.8961%；
B-cubed F1 V2=75.4443%，V1=77.2541%。结论依指标不同，不宣称全面优胜。
V2图组件仅丢23参考对，但最终丢176130对；需固定现有树诊断内部树边界与OG提取，
未据此改变生产默认或添加过滤规则。完整方法及评分不再重复提交。

## 固定图/树的边界与选择诊断

按用户明确要求：不调gamma、不改变soft42默认、不重做search/edges/hierarchy。
新增qfo_bacteria_fixed_tree.sh：4CPU/16GB/2h，读取1411893固定V2图/树及
1411894 OF3参考，使用已有run_audit与pair_f1_oracle接口。
比较terminal cut、生产OG、全部节点完整cut oracle、事件资格约束cut oracle；
已有run_audit校验所有输入与节点/成员哈希、完整覆盖、fractional全局最优性。
标签oracle只保留独立DIAGNOSTIC文件，绝不写回生产结果。

另用节点参考标签计数及(reference,production OG)联合计数计算pair LCA损失，
精确分解当前丢失对为跨组件、结构叶内部、纯assigned参考LCA、混合assigned LCA。
纯LCA表示树内有可用完整祖先边界；混合LCA的整个节点合并需夹带其他assigned家族。
这两项是可用机会/表示限制诊断，不构成生产政策的因果归因；unassigned不参与纯度定义。
图root总能提高recall但代价是误合并，因此不把recall差或oracle差硬拆为因果贡献。
完整cut家族与生产active-view提取不等价，oracle不是任意OG输出的上界。
计数算法不枚举真实蛋白对；小型随机树用显式pair枚举核验所有分类之和，
同时覆盖纯LCA、混合LCA、结构叶内部损失及非法重复成员。
再次核验固定图与树哈希后输出boundary-summary.json；大per-OG及oracle cuts留远端。

## job1411949固定树诊断完成

COMPLETED / 0:0，3:45；哈希核验通过，结果归档于
OGProfiler2_benchmark/05_metrics/qfo/bacteria_fixed_tree_1411949。
生产/全节点cut oracle/资格约束cut oracle的pairF1=63.5211/81.5008/66.0049%。
资格约束oracle的recall52.6112%低于生产54.2763%，F1提升来自减少误合并；
不将两个F1差当作树/选择的因果贡献。
当前丢失176130对分为跨组件23、结构叶内0、纯assigned LCA52121、混合LCA123986；
component0占全部损失73.5803%。任意节点完全匹配5723OG，资格节点4739OG。
下一步固定树关联纯LCA与事件资格/生产active-view轨迹，核查资格与选择路径，
同时诊断混合LCA的家族交错；不扩大gamma、不改默认、不自动重跑方法。

## 纯LCA资格与生产路径追踪

新增pure_lca_trace只读入口，读取1411949已核验的report与boundary-summary、
固定层级节点/成员，以及生产og-manifest、v1_events、selection_trace、groups。
按节点聚合纯LCA丢失对，严格核对总数等于52121的既有诊断（动态读取基准）。
资格规则与生产相同：None，或多物种I。首先区别资格排除与资格通过；
后者再分已选但参考家族仍碎片化、未选且有消耗轨迹、没有可见选择轨迹。
已选与后续消耗可能同时发生，因此已选优先；不把所有资格通过的损失误称未提取。
保留trace_order/processing_level/consumed_by及源OG规模，便于核查active-view。
所有既有树/图和生产输入输出哈希再次检查，大逐节点记录与哈希清单留远端。
使用独立4CPU/16GB/2h只读Slurm作业，不重做方法、不改默认、不扩大gamma。

## job1411967纯LCA生产路径完成

COMPLETED / 0:0，44秒，固定输入核验通过。纯LCA损失总数52121，涉及2360节点，
全部被事件资格规则排除；资格通过但未提取与已选仍碎片化两类均为0。
事件II/III-1/III-2/III-3分别贡献101/18576/10200/23244对。
结论限于纯assigned LCA子集，不推广至全部损失；节点可嵌套，不可按纯度直接全放行。
下一步需固定图树采集资格排除节点及混合节点对照，确定标签无关资格/选择实验，
完整划分同时评估召回收益与误合并代价，暂不修改默认或扩大gamma。

## 固定树资格结构证据对照

qualification_structure以不含参考标签的函数计算各内部节点结构特征：
直接子树间/内部边权密度（缺失边按0，内部无可比对或权重0时比值置空）、
节点外部边权占incident strength比例、最小子树规模比例与子树物种重叠。
采用LCA边权累积和后序传播，不逐个节点扫描全图、不展开蛋白对。
参考标签只在计算特征后分组：纯LCA有损失、纯节点无损失、混合节点、全unassigned。
比较所有资格排除内部II/III节点，不仅挑已观察损失或top案例；核对纯LCA52121对。
事件、log2基因数、log2物种数、component0/其他分层，报告分位数及缺失比例。
嵌套节点并非独立样本，不做假独立显著性结论；不根据纯度或经验分位数拟合生产阈值。
结果出来后再选择单一明确的标签无关资格对照；不直接开放全部II/III，
不以低覆盖/高集中度做过滤，不扩大gamma或修改生产默认。
独立只读Slurm资源4CPU/16GB/2h，使用原固定图树与1411949基准报告。

## job1411988结果与资格对照设计

COMPLETED / 0:0，1:21，哈希检查通过。2360纯损失/1589纯无损失/6907混合/1全unassigned节点。
61可匹配分层覆盖纯损失2356节点、51951对；跨/内部密度52层区间重叠，纯中位数仅18层更高。
III-3大损失层与多个III-1层方向不同，证据不支持通用单指标阈值或全面开放II/III。
报告归档bacteria_structure_1411988，提出尚未实现的MSA-v1资格实验：
二叉多物种排除节点仅在两子树相互连接密度严格超过各自到所有祖先兄弟分块密度时附加资格。
参数/标签无关，原资格与active-view排序不变；先验证生产基线回放，再评估完整划分净TP/FP变化。
这是可证伪假设，效果未知，不修改默认，不自动提交策略实验。

## MSA-v1实验实现

benchmark-only适配层暂时替换共享提取器的内存annotation视图，把通过MSA的排除节点
作为I候选；原事件文件不改，恢复函数绑定，生产代码与默认配置未编辑。
MSA只读取图树和原事件：非根二叉多物种II/III节点的兄弟密度必须严格超过
双方到每层祖先兄弟分块的密度，正边权、严格平局排除，不读OF3标签。
稀疏子树/祖先兄弟权重以edge*height累积，小随机树显式密度核验；根/原资格不变。
逐组件先回放原提取，成员集合必须与冻结生产划分完全相同，事件必须与存档相同。
然后独立提取MSA完整划分；重叠/未覆盖组件记录失败，存在失败则不输出全量分数，
不偷偷修补或退回生产结果。完整通过后才加载参考标签评分，报告新增/移除TP/FP。
哈希检查、候选原事件及资格证据、实验members均存独立输出；不修改树/gamma/默认。
资源4CPU/16GB/2h，使用同一冻结V2/OF3/1411949基准，作为开发集实验而非官方QFO准确率。

## job1412011 MSA完整对照结果

COMPLETED / 0:0，2:12；生产划分与事件回放、完整覆盖、冻结输入均通过。
5473/7522节点通过资格：III-1 4177、III-2 1263、II33、III-3 0。
P/R/F1由76.5618/54.2763/63.5211%变为27.3530/74.7020/40.0436%；
新增TP78681、FP700248，无移除TP/FP。最大组440assigned蛋白/73参考家族/94839FP。
MSA假设未带来净改进，不采用默认，也不继续其阈值搜索；生产参数保持不变。
admitted_selected字段报表缺陷忽略，修正在8723180，不影响成员与评分。
下一步如继续，应固定树定位新增误合并的source_cluster与active-view消耗路径，
对照新增TP无新增FP的组，区分资格与提取范围机制；本轮不自动提交作业。

## MSA新增误合并来源与active-view路径

只读取1412011保存的candidates、members、summary与input-hashes，不重新计算MSA资格。
调用同一共享提取器回放，完整分区和component:local_group ID必须与既有实验一致。
逐有新增TP/FP的组计算source_cluster、原事件、资格证据、处理层次和consumed_by轨迹，
比较原完整子树与实际输出的inside/outside/omitted成员，区分完整来源节点放行与范围变化。
分别报告实际输出与假设完整source的新增TP/FP，不将后者作为可同时实现的划分。
逐组新增TP/FP求和必须严格匹配1412011的78681/700248；归档top误合并与无新增FP的TP对照。
新增active-view诊断不修改默认、图树或资格规则，不自动尝试新评分/阈值。
独立只读回放Slurm资源4CPU/16GB/2h，全部输入前后哈希核验；大轨迹留远端。

## job1412049路径定位完成

COMPLETED / 0:0，41秒；1412011完整分区/组ID回放与哈希检查通过。
全部新增TP78681/FP700248来自1205新资格完整来源节点输出，无新增对来自scope变化。
III-1 872组/50403TP/431416FP；III-2 326组/26617TP/235692FP；II7组/1661TP/33140FP。
最大误合并0:0对应原cluster22807，原III-1，完整458=实际458，新增FP94075，消耗786后代。
清洁收益21394/10040/48967也都是完整来源提取，说明事件/路径本身不区分污染。
结果支持定位资格放行污染祖先节点，而非越出原子树的active-view扩张；范围限于新增关系。
下一步如继续应检验内部家族边界与标签无关完整合并/切分目标，不修补active-view范围或继续MSA阈值。
默认、图树与gamma保持不变，本轮仅归档，不提交新实验。

## 完整来源节点内部合并/切分目标诊断

固定1412049的1205个有新增对且互不重叠的完整来源节点，复用graph_cut接口。
节点诱导子图上固定resolution1的strength-null目标：sum Win/W-(strength/(2W))^2。
根合并得分0；同一节点子树全部节点可选，exact DP输出完整无重叠cut，平局保留父节点。
零边权保持显式0目标；不寻找gamma或拟合阈值。外部边不进入该局部诱导目标，局限明确。
所有图树全蛋白先算cut，再读OF3 assigned标签；事后分污染/清洁收益及事件，
报告cut gain、切分数量、相对完整合并移除TP/FP与相对生产保留新增TP/FP。
此为目标鉴别力诊断，不是生产策略、祖先推断或完整数据替代划分。
独立4CPU/16GB/2h只读Slurm，输入前后哈希校验，大case和DIAGNOSTIC cuts留远端。

## job1412069内部目标结果

COMPLETED /0:0，22秒；原合并计数、完整节点覆盖与输入哈希核验通过。
765污染组678拆分，440清洁收益组235也拆分；移除合并FP629443，但丢TP57493。
保留相对生产新增TP29375/FP79040（仅1205节点，非全量方法评价）。
清洁III-2节点10040 cut gain0.50193仍丢333TP，说明高增益不专属于污染边界。
目标有污染结构敏感性但家族边界鉴别不足，不接默认、不搜索resolution/gamma。
下一步如继续应固定cut检查物种结构与局部null期望的增益来源，不增加过滤或修改提取范围。

## 固定cut物种组成与null期望诊断

读取1412069保存的cuts/cases与哈希清单，不重算cut、不调resolution。
使用完整蛋白物种元数据和诱导边权，按子群计算物种计数与incident strength。
跨cut子群的null期望Wij=si*sj/(2W)，再用各子群同物种强度分解same/cross-species期望，
对实际跨子群边权作同样分解，严格核对sum(expected-observed)/W等于存档cut gain。
保留负贡献与零边权，不把单个正物种贡献视为因果机制。
分别报告清洁/污染各事件的物种贡献分布、子群物种重叠和10040/21394/22807/14737案例。
标签只用既有cohort，不参与特征或资格；默认配置、原图树与提取范围保持不变。
独立只读Slurm4CPU/16GB/2h，结果小摘要取回，大逐节点统计留远端。

## job1412073结果

COMPLETED /0:0，19秒，固定cut gain与输入核验通过。
清洁拆分225/235与污染646/678节点跨物种贡献为正；该信号不专属于污染。
清洁21394的一组26基因覆盖23物种，不支持简单按物种分群解释。
固定原cut继续物种对条件化null敏感性诊断，不优化cut、不调resolution、不改默认。

### 下一步：固定切分的物种对条件 null 诊断

在 job1412069 保存的全部 1205 个 source cuts 上，新增显式 opt-in
`--condition-species-pairs`；原诊断入口默认行为及生产 soft42 均保持。
每个 source induced graph 内，对无向物种块 `(s,t)` 保留总边权 `W_st`，
并保留每个切分组在此块上的端点强度 `k_i,s,t`。跨物种块使用二部配置期望：
`E_ij,st=(k_i,s,t*k_j,t,s+k_j,s,t*k_i,t,s)/W_st`；同物种块使用
`E_ij,ss=k_i,s,s*k_j,s,s/(2W_ss)`。零权块贡献零；不添加伪计数。
以 `sum_i<j,st(E_ij,st-observed_ij,st)/W` 重新评价既有切分，保留正负增益。
同时验证各物种块的端点强度守恒及原无条件增益与保存值一致。

输出逐节点、逐物种块的期望与实际跨组边权，以及 clean/polluted、原事件分层的
增益变化和符号翻转。`1e-10` 仅用于报告数值正负容差，并非资格过滤阈值。
参考标签仅用于事后队列汇总，不进入特征或 null。该对照是 null 敏感性诊断，
不优化切分，不产生生产完整分区，不把正增益解释为生物学真值。
Slurm 延用 4 CPU /16G /2h 的单作业资源，无参数矩阵。

实现验收：本地完整轻量测试 484 passed /8 skipped /5 项既有 warning；
远端合成短验证 12 passed；ruff、diff whitespace、Slurm bash 语法检查通过。
提交前确认冻结输入可读，正式诊断使用 clean committed source snapshot。

### job1412081 完成：物种对条件 null 的固定切分结果

COMPLETED /0:0，43 秒。RUN_ID=20261008T051058Z_6910e701d7c1_8012db3f_2327，
source=6910e70，dirty=0。固定 1205 个 source cuts，哈希与原增益一致性、
物种块强度守恒通过；无生产改动。
已切分清洁节点 235 中 133（56.60%）由原正增益变非正，
污染节点 678 中 147（21.68%）翻转。比例是节点统计，不是加权 TP/FP 修复率。
清洁 21394 增益 .134134→.078713，10040 .501930→.221222；
污染 22807 .774820→.657379，14737 .528123→.351742。
条件化具有差异化评分影响，但仍有 102 个清洁已切分节点正增益，
尚未建立新的完整最优切分或验证生产收益。默认 soft42 保持。
紧凑归档：OGProfiler2_benchmark/05_metrics/qfo/bacteria_species_pair_null_1412081/。
下一步候选：固定树、完整可加和条件目标的动态规划切分对照，
不是事后符号过滤或分辨率搜索。

### 固定树的完整物种对条件 merge/cut 目标

新增 benchmark-only `species_pair_cut`，不是对旧 cuts 的符号过滤。
在每个原 source 的完整诱导图上固定所有 `W_st` 和蛋白端点强度，向上聚合
`k_v,s,t`；节点 keep 分数为
`q(v)=[Win(v)-sum_s k_v,s,s²/(4W_ss)-sum_s<t k_v,s,t*k_v,t,s/W_st]/W`。
分母在 source 内固定，不随候选子节点重新估计，否则失去可加和与同尺度比较。
零权块省略；总权零时所有节点分数零。每个内部节点递推
`best(v)=max(q(v), sum_child best(child))`，平分保留父节点，结构叶不可再分。
归一化分数绝对值不超过 1e-12 视作浮点零，不设置可搜索的生物学门槛。

独立以所选完整分区的跨组块期望减实际边权重算目标，核验与 DP 分数一致。
重算无条件目标的最优 cut，并验证新 DP 分数不低于该可行 cut 在新目标下的分数。
全部切分完成后才加载参考标签，输出按 cohort/event 的完整 merge/cut TP/FP 损失、
相对无条件最优 cut 的 TP/FP 变化，以及对原生产组的保留新增 TP/FP。
每个 source 完整覆盖、source 间不重叠、输入哈希保持；范围仍为 1205 个变更 clades，
不是整套 82507 蛋白的替代生产分区。不改变资格或 active-view 范围。

合成测试包含穷举所有完整 cuts 与 DP 最优性对照、独立目标恒等式、单物种退化、
物种分组的零增益父节点 tie、混合物种边界、缩放、零权及输入边验证。
单 Slurm 作业沿用 4 CPU /16G /2h，无参数矩阵；生产 soft42 不变。

本轮实现验收：本地完整轻量测试 489 passed /8 skipped /5 项既有 warning；
远端合成测试 15 passed，ruff、diff whitespace 和 Slurm bash 语法通过。

### job1412089 完成：条件化完整 DP 切分并未解决关键清洁节点

COMPLETED /0:0，39 秒，RUN_ID=20261008T052424Z_5bb4e523a0ca_ab885050_14156，
source=5bb4e52，dirty=0。1205 个 source 的完整切分与独立目标恒等式、覆盖、
输入哈希、旧 cut 可行值比较检查通过。仍非整套生产替代分区。
clean 440 中切分节点235→109，移除merge TP12448→7873；
polluted 765 中切分678→557，移除TP45045→41027、FP629443→620867。
相对原无条件最优 cut 合计多保留8593 TP、8576 FP，存在权衡。
相对原生产新增且保留的 pairs 为33919 TP、84498 FP，是另一统计量。
关键清洁source21394由3→5 groups，移除TP305→352，新TP89→42；
10040由6→7 groups，移除TP333→342，新TP54→45。
10040的新最优总目标与旧 cut 在新目标下重评相同，需检查局部 tie 路径，
暂不把相同总分直接当作数值误差证明。21394则有真实报告增益变化但更低TP。
优先审计固定新旧 cuts 的局部增益与数值平分，再讨论真实边界目标；
默认soft42保持，无阈值、gamma或分辨率搜索，本轮不追加作业。
紧凑归档：OGProfiler2_benchmark/05_metrics/qfo/bacteria_species_pair_cut_1412089/。

### 10040 同分异切与 21394 正增益边界的定点审计

新增 opt-in `--tie-replay-dir`，只读取 component0/source10040、21394。
先逐蛋白核验当前浮点条件目标分区与冻结 job1412089 分区一致；
所有树图输入仍用原 source、原诱导权重，不改变任何 score 或默认路径。
精确审计使用 `Fraction.from_float` 保留存储 binary64 边权的完整值，
逐节点重算内部权、物种块端点强度与原可加和目标，精确判断 keep/split 的
差值符号和等号。不使用 epsilon、分辨率调参或参考标签决定平分。

输出浮点与精确 DP 的局部 keep/split 最优值、差值、choice_changed、
两套 active 路径及 exact tie 标志；精确最终分区独立核验跨组 null deficit 恒等式，
并验证不劣于原浮点 cut 在精确目标下的值。精确分区是诊断反事实，非生产修正。
逐节点报告孩子最优 cut 的同/跨物种边界实际权重、期望与 signed gain；
参考标签随后仅用于 keep 与孩子最优 cut 的 TP/FP 变化及精确最终 cut 的指标。
原浮点选择函数与 job1412089 快照不修改。双目标定点 Slurm 沿用4 CPU /16G /2h。

验收：本地完整测试492 passed /8 skipped /5项既有warning；
远端定点合成测试11 passed；ruff、diff whitespace、Slurm bash语法通过。

### job1412103 完成：数值tie与真实边界增益已区分

COMPLETED /0:0，13秒；RUN_ID=20261008T055443Z_53fa3e8ce4f6_3daf636e_31238，
source53fa3e8，dirty=0。10040/21394浮点分区逐蛋白回放、精确目标恒等式、
覆盖与输入哈希通过。10040活跃10052浮点split−keep=2.7756e-17，
精确差0，parent-on-tie恢复9 TP，7→6组、TP93→102，精确目标不变；
但相对merge仍lost333，未解决全部过切。非活跃10046的tie不影响实际分区。
21394精确分区与浮点相同，5组TP278 lost352；活跃21396精确正增益
.0001226236115487494，跨物种边界actual20.061502195660537 vs
expected20.1140384784857，keep TP300 vs childcut253，真实损失47 TP。
非活跃21399虽有数值tie，对该例输出无贡献。嵌套节点局部TP损失不累加。
下一步分开修复实验数值tie与研究21396真实边界目标，避免可调epsilon
抹去真实小增益。默认soft42保持；紧凑报告归档bacteria_conditioned_tie_1412103。

### 修复实验条件 DP 的精确平分语义，并分解21396物种块

实验 `species_pair_cut` 新的规范路径使用精确有理数 keep 分数传给原完整树 DP；
仅在严格数学相等时保留父节点，不进行 1e-12 分数置零或 epsilon 比较。
所有权重以 `Fraction.from_float` 表示存储边权，分母固定 source；
LCA 内部权与端点强度使用精确累加、后序聚合，避免逐节点重扫全部边。
最终跨组 null deficit 独立精确重算，恒等式用等号核验。

历史浮点分支保留为显式 `exact=False`，定点审计和既有 benchmark 入口继续显式
回放历史分支；修复回归入口 `--exact-conditioned --raw-replay-dir` 使用新路径。
逐蛋白核验全部1205 source 的旧浮点 cut 与 job1412089 一致，再比较精确 cut 的
TP/FP 与旧 cut，历史快照与文件不改写。生产 src/config/defaults 不修改。
新增合成回归构造 `nextafter(1,0)` 边权：真实增益低于1e-12仍应split，
并验证非零keep分数上的真tie、边顺序不变性、零权、审计一致性及目标恒等式。

单 Slurm 作业内先进行全部1205 source的精确修复回归，成功后再定点审计
10040/21394，并为21396的孩子最优边界输出各物种块的 source W_st、端点强度、
观察跨组边权、期望及 signed gain。仍固定现有诱导图和树；分解用于定位真实
小正增益，不添加阈值或修改资格。资源沿用4 CPU /16G /2h，无参数矩阵。

修复验收：本地完整测试496 passed /8 skipped /5项既有warning；
远端轻量合成测试15 passed，ruff、diff whitespace、Slurm语法通过。

### job1412104 完成：精确修复回归通过，21396增益来自唯一物种对

COMPLETED /0:0，69秒；RUN_ID=20261008T061431Z_ff314c0b1b9e_737ab568_40915，
sourceff314c0，dirty=0。全部1205 source旧浮点分区回放、精确目标恒等式、
覆盖、哈希通过。相对旧浮点条件DP：clean净+312 TP/0 FP，polluted净+739 TP/+732 FP，
合计+1051 TP/+732 FP。根切分109/557不变；不归因于纯tie，更不宣布生产改善。
10040恢复9 TP到102、6组；21394仍278 TP、5组、lost352。
21396三子切分为24+1+1（两个单例31672、31901）。255物种块中254零贡献，
唯一正块(1,11) expected .5885919586270754 vs actual .5360556758019136。
核心—31672的期望 .4673649989732148，核心—31901的期望 .12122695965386063；
actual仅有块总量，尚不能据此分配逐子组对实际边权。此前对子组实际强弱的
口头解释证据不足，报告已明确更正，需读取原始边核验。
node外species1端点强度约 .2660848122536832已由块总量确认，
具体单例内外连接方向仍待查。树为三子同时切分，未提供核心+单例中间节点。
下一步定位原始边流与完整切分表示贡献，不增加门槛。零贡献不是缺边证据。
soft42默认保持。
紧凑归档：bacteria_exact_dp_1412104；本轮无追加作业。

### 原始冻结边流与三子完整切分表示的核验

新增只读 `--singleton-flow`，只在已有定点 replay 上启用。
扫描冻结 component0 的所有 retained weighted edge 记录，捕获31672/31901的
全部incident edges（含source外、component内），保留u/v、实际权重、原蛋白ID与物种ID。
这里的“原边”指冻结图的原始记录，不是预过滤相似性hits，不补造缺失边。
逐单例区分：node内直接子组、node外但source内、source外；若两个单例直接相连，
该无向边会分别出现在两个单例的incident报告中，不把两份报告相加当作总边权。
source内incident记录与source诱导图逐边/逐权精确一致才继续。

对各直接子组对、各物种块独立重算actual和固定source null期望及signed deficit，
核验三组总目标与已有精确局部增益一致。列出三子集合全部五种集合分区，
按原树节点蛋白集合逐一核验可表示性；部分子组合并若无原节点，明确标记为
仅解释性反事实，不作为新生产cut或新选择策略。
参考标签只在结构/边流/候选分区固定后事后计分。
保持现有树图、精确DP与历史浮点回放、资格、active-view、soft42；无门槛搜索。
Slurm单定点作业沿用4 CPU /16G /2h。

验收：本地完整测试499 passed /8 skipped /5项既有warning；
远端合成测试9 passed，ruff、diff whitespace、Slurm语法通过。
