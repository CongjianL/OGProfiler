# QFO比较分析：soft42生产默认（2026-10-06）

用户明确要求继续使用 /home/mselab/licj/project_data/OGProfiler/QFO 做benchmark。
本地源码权威，用户目录只读；所有全数据处理及计算通过Slurm。

## 已检查事实

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
