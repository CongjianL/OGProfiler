# 图2：默认 V2 与 soft42 评估更新

2026-10-06；用户要求将两项现有结果接入OGProfiler2_benchmark图2。本轮不重新调参或提交科学作业。

## 来源与口径

- 当前默认V2：job1410777，RUN_ID 20261002T060705Z_cb0c0ba98049_1cb4fcd7_14702。
  有效hierarchy配置逐项与当前DEFAULT_CONFIG一致（kway、depth20、10轮）。
- soft42：job1411546，RUN_ID 20261005T232747Z_8bcccc9d3302_512d821e_28231。
- 两项都为固定mean-SSN的全组件hierarchy/OG结果，v1_compatible提取；不是重跑搜索。
- 原图OGProfiler2=早期V2，保留并标注V2 legacy；另保留V1 first和四个竞争方法。
  本地原图含V1 first而远端旧图尚未包含，因此同步后两端统一为八个条件。

|条件|官方P %|官方R %|官方F1 %|raw median split excess|top10 FP占比|
|---|---:|---:|---:|---:|---:|
|V2 default|91.7440|25.5911|40.0192|7|0.6333915|
|V2 soft42|92.6961|33.8220|49.5608|4|0.5988063|

soft42相对默认V2：P+0.9521、R+8.2309、F1+9.5416个百分点。F1仍低于完整V1。
这些变化联合涉及topology/search schedule和depth，不作单因素因果归因。

四面板同步：A官方P/R，B每个方法70 raw RefOG best-group F1，C raw split excess和
best-group contamination，D官方归一化FP的top1/5/10集中度。raw指标沿用原图实现，
与前次报告中排除低置信成员的confident诊断分开；D仍按官方低置信排除。
新FP分解计算的P/R/F1与官方报告逐项差异小于1e-10。原六个条件的源数据行保持原样。
原Stage2A统计/CI和报告保留历史含义；本轮不声称已扩展这些推断统计或运行时对照。

## 文件与复现

本地与远端benchmark目录的 `07_figures/Fig2_candidate/` 更新PNG、PDF、图注与data。
图尺寸183x160mm（原145mm高度增加以容纳八个标签），600dpi，已视觉检查。
新评估表另存 `05_metrics/orthobench/v2_configuration_update/`，不覆盖历史B1分析表。
脚本 `OGProfiler2_benchmark/scripts/statistics/update_fig2_v2_configs.py` 对每个输入校验
评分验收、预测哈希、完整覆盖、配置类型；元数据JSON记录job/run及输入输出哈希。
预测组只取回两份约5MB文本；远端大树/SSN和原始实验不移动。

执行：
```sh
.venv/bin/python OGProfiler2_benchmark/scripts/statistics/update_fig2_v2_configs.py \
  --base PRE_UPDATE_FIGURE_DATA --default-run LOCAL_JOB1410777 \
  --soft42-run LOCAL_JOB1411546 --out UPDATED_FIGURE_DATA --metrics-out UPDATED_METRICS
Rscript OGProfiler2_benchmark/scripts/plotting/plot_stage2a_fig2.R UPDATED_FIGURE_DATA OUTPUT
```

本地原图保存在 `.benchmark_cache/fig2-update/local-before/`；远端原图在
`07_figures/Fig2_candidate_archive_before_v2_configs_1411546/`保存完整副本后再同步。
新增单元测试3项通过；全部四张源表原始行逐项相等，RefOG长表420→560行。
图2更新溯源另以Git归档：`docs/fig2-v2-config-update-provenance.json`。

## 交付验收

- 完整本地测试：413 passed / 8 skipped（5项既有warning）。
- 复现脚本重新生成的四张数据表逐字节相同，provenance JSON语义相同。
- 远端旧图副本创建成功；新版PNG/PDF、图注和六份data文件全部sha256sum校验OK。
- 图与新指标目录已同步回OGProfiler2_benchmark；未覆盖原B1官方/统计表或原始运行。
