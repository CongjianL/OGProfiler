# H5 soft+depth42：剩余低召回及新增污染首次偏离审计

Date: 2026-10-04。当前冻结 job **1410868**，对照同 mean SSN 的 kway job **1410777**。

## 范围与可复现性

只读审计，不运行 Leiden、SSN 构建、OG 选择或官方 scorer；停止尺寸1、ARI .9、
深度42实验配置及 OG/scorer 全部保持原样，生产默认不迁移。
本地只取回67个相关组件的诊断表、原始ID/组件映射、最终members与官方诊断；
大型 SSN 及全组件 hierarchy 保留远端。hierarchy/OG表按manifest逐项核对hash。
对照时校验蛋白与组件映射字节相同、RefOG/低可信文件hash相同；
比较成员clade和子分区，绝不比较两棵树的numeric ID是否相同。

覆盖全部70 RefOG、1,895可信蛋白、39,242真值无序对。按官方q_even的
`1/(n_confident-1)`权重重建recall **.3382199840850459**，与官方值误差小于1e-12。
诊断原始对数、best-group recall及官方加权recall分别报告，不混用。

本地证据目录：
`.provenance/slurm/20261003T120507Z_27a2533e1334_72b113e6_93925/results/`

- `refog-split-audit/{summary,refogs,first-split-nodes,examples}.json`、`refogs.tsv`
- `refog-pollution-audit.json`：全部RefOG重叠预测组的新增跨真值对、首次成员分区偏离、LCA和选择trace。
- `refog-candidate-coverage.json`：5个代表性节点的全部冻结候选。

工具：现有 `benchmarks/og_extraction/refog_split_audit.py`，新增
`benchmarks/og_extraction/refog_pollution_audit.py`。前者兼容旧报告布局，本次本地
`h5-report-compact.json`只封装实际 `h5-metrics/report.json`，未修改数值。
污染工具CLI：`--current CURRENT_RESULTS --baseline BASELINE_RESULTS --out REPORT_JSON`。

## 1. 剩余召回损失首次在哪一层发生

| 去向 | 原始真值对 | 官方权重对应百分点 |
|---|---:|---:|
| SSN跨连通组件，最终未保留 | 6,046 | **12.808652** |
| 同SSN组件，hierarchy分开后未重新合并 | 25,274 | **53.369350** |
| 最终保留 | 7,922 | **33.821998** |
| 合计 | 39,242 | 100 |

SSN项与1410777完全相同：35/70 RefOG跨组件，当前SSN的乐观关系保留上界仍为
87.191348%。新增召回收益发生于组件内部，而非SSN修复。
当前结构叶直接保留266对（1.090138个百分点），OG跨叶重新合并7,656对
（32.731861个百分点）。因此结构singleton不等于最终OG碎片。

**全部25,274个组件内丢失对没有共同的合格I/None祖先**；没有同一结构叶被OG再拆散、
没有共同合格祖先因SKIPPED_CONSUMED/排序而漏选的案例。
这是当前树与现有V1事件资格的衔接问题，不是已通过parity的OG实现偏离V1规则。
其首次物理分散节点事件：

| LCA事件 | 原始丢失对 | 当前损失百分点 | 1410777损失百分点 |
|---|---:|---:|---:|
| III-3 | 15,081 | **25.869326** | 40.716326 |
| III-1 | 8,229 | **22.995547** | 18.580321 |
| III-2 | 1,187 | 3.018766 | 2.047696 |
| II | 777 | 1.485710 | .255928 |

III-3显著减少，但binary带来的III-1/III-2/II仍会失去共同I祖先。
ARI=1只证明重复划分稳定，不证明该划分具有正确正交群边界。

## 2. 代表性低召回RefOG的具体路径

下表节点为该RefOG最终丢失对中**最浅的组件内首次分散节点**；同一RefOG不同蛋白对
有各自LCA，完整分布见JSON。跨组件的对更早在SSN阶段分散。

| RefOG | SSN组件数 | SSN丢失对 | hierarchy后丢失对 | 最浅节点component:node / depth | 分裂及事件 |
|---|---:|---:|---:|---|---|
| 021 | 12 | 1,421 | 6,114 | 0:0 / 0 | 69,642蛋白，470-way KWAY，III-3 |
| 019 | 2 | 12 | 65 | 3976:0 / 0 | 12蛋白，11-way FALLBACK_KWAY，III-3 |
| 007 | 1 | 0 | 1,449 | 0:28399 / 2 | 331蛋白，11-way FALLBACK_KWAY，III-3 |
| 011 | 7 | 2,160 | 1,787 | 153:0 / 0 | 89蛋白，3-way FALLBACK_KWAY，III-3 |
| 012 | 1 | 0 | 804 | 470:0 / 0 | 46蛋白，binary，III-1 |
| 023 | 1 | 0 | 2,092 | 0:37393 / 3 | 201蛋白，18-way FALLBACK_KWAY，III-3 |
| 058 | 2 | 225 | 845 | 49:0 / 0 | 180蛋白，3-way FALLBACK_KWAY，III-3 |
| 003 | 1 | 0 | 1,028 | 0:13894 / 1 | 1,169蛋白，38-way FALLBACK_KWAY，III-3 |
| 033 | 1 | 0 | 365 | 0:16872 / 1 | 114蛋白，binary，II |
| 006 | 4 | 138 | 829 | 7:341 / 2 | 182蛋白，8-way FALLBACK_KWAY，III-3 |
| 036 | 2 | 45 | 818 | 500:0 / 0 | 45蛋白，6-way FALLBACK_KWAY，III-3 |

007还有0:28731(depth5)的binary/II，使510对首次分散；不能只看其上游fallback。
021的大根仍在>=10000原kway作用域；后代0:84730的129蛋白12-way fallback另损失1,991对。

### 候选覆盖与准入冲突，不直接修改门槛

- **019/3976:0**：已测24点，2个原始binary均未合格，8个候选UNSTABLE；
  原门槛合格kway为2个，最终gamma1.12、ARI1的11-way。不是深度截断，也不是未执行搜索。
- **058/49:0**：24点原始count全>=3，未测得binary；15点通过原门槛。
  最终选gamma_min=.01的3-way。此证据仅证明有限日程没有采到binary，
  不证明区间内部或下界外存在合格解，也不批准修改下界。
- **036/500:0**：24点中10个原始binary，无binary合格；17点MAX_CHILD_FRACTION、
  5点NO_SPLIT+MAX_CHILD_FRACTION、1点UNSTABLE，唯一合格是6-way。
- **023/0:37544**：77蛋白/7-way/III-3为该RefOG最大组件内丢失LCA（1,303对）；
  24点中8个原始binary、零binary合格；12点比例失败、9点不分裂且比例失败、2点不稳定。
- **012/470:0**：7点中已有合格binary，gamma.08、ARI1、最大子群比例.5652。
  事件却是III-1；根分散520对，后代分别175/99/10，合计804。
  最终保留231/1035对，与旧树相同。这里继续寻找“任意binary”不是充分修复目标。

上述失败原因可重叠，原始binary count并不意味着每seed都binary，也不是精确无解证明。

## 3. 新增污染：完整预测组扫描，而非只看最佳组

相对1410777，按每个RefOG去除其低可信成员后，对**所有与可信成员重叠的预测组**
检查旧树分组不同、当前分组相同的“可信成员—组外成员”对。
总计7个RefOG、10个预测组、**146个新增跨真值对**（该统计不是官方precision分解）。
其中141对在OG选择阶段跨不同结构叶重新合并；5对已同属一个ONE_SPECIES结构叶。
“额外成员”是benchmark相对标签，不直接证明生物学非正交关系。

| RefOG | 新增跨真值对 | 当前共同OG来源component:source | 机制 |
|---|---:|---|---|
| 007 | 3 | 0:28764 | binary/I跨叶合并 |
| 011 | 5 | 153:133 | ONE_SPECIES六蛋白叶，None残余选择 |
| 015 | 23 | 943:0 | binary/I跨叶合并；新增OVERMERGE类别 |
| 020 | 48 | 362:12 | binary/I跨叶合并 |
| 021 | 19 | 958:2、0:84878、0:84999 | I跨叶合并，其中0:84999为KWAY |
| 023 | 37 | 0:37587(34对)、0:37703(3对) | binary/I跨叶合并 |
| 035 | 11 | 84:9 | binary/I跨叶合并 |

### 三个最佳组污染增加案例：首次拓扑偏离 → 最终错误合并

- **015**：固定component943根的同一31蛋白/5物种clade，旧node0为3-way/III-3，
  新node0为binary/I（gamma.05）。这就是首次成员分区偏离与新选择资格位置。
  额外 `ENSDARP00000135749` 与23个可信成员原本分组不同；当前仍在不同结构叶，
  部分对在node4/III-3已物理分散，却被根I最终重合并。根trace明确SELECTED。
  原组共31蛋白包含7个该RefOG低可信成员，诊断扣除后为23可信+1额外。
- **023**：首次成员分区偏离在固定343蛋白clade：旧0:30753(depth2)11-way/III-3，
  新0:37392(depth2)8-way FALLBACK_KWAY/III-3。两者事件相同但membership分区不同。
  下游新0:37587(depth12)的19蛋白/10物种binary/I被SELECTED，
  合并17可信与 `ENSCINP00000035305`、`FBpp0075495`，增加34个跨真值对。
  旧树两类成员的LCA为0:30952(depth12)15-way/III-3。
  新选择源19蛋白clade在旧树没有完全相同节点，不拿numeric ID作同节点对照。
- **035**：首次偏离在固定14蛋白clade：旧84:7(depth7)4-way/III-3，
  新84:7(depth7)binary/III-1。之后新84:9(depth9)12蛋白/10物种binary/I被SELECTED，
  11可信与 `WBGene00000020.1` 合并。旧这11对LCA在84:7/III-3；
  新source的12蛋白clade旧树不存在。

补充：**011**并非跨叶错误重合并。首次分区偏离在153:1(depth1)固定79蛋白clade：
旧3-way/III-3 → 新binary/III-1。后代153:133(depth4)聚成6个线虫蛋白的ONE_SPECIES叶，
其中可信 `WBGene00015545.1` 与5个组外成员由None残余组一起发布。
维持当前ONE_SPECIES停止语义时，OG没有再拆结构叶；此例需归于新拓扑下的同物种叶混合，
而不是I事件选择漏检或新增停止尺寸修改。

## 4. 结论与后续设计入口

1. SSN阶段仍有12.81个百分点不可跨组件恢复的损失，当前政策实验没有改变这一项。
2. 剩余53.37个百分点主要是“hierarchy分裂 → 没有合格共同I/None祖先”；
   III-3与binary/III-1均应独立处理，继续增深不能解释这些已resolved且已搜索的节点。
3. 新增污染表明“获得I资格”本身不是正确OG的保证；更早改变的成员拓扑、选择资格
   及ONE_SPECIES后代应一起审计。不是把全部binary自动替代fallback的证据。
4. 下一步应分别设计覆盖不足（019/058）、原门槛冲突（036/023）、
   合格binary但事件不适用（012）、新I污染（015/023/035）及同物种叶混合（011）的
   **冻结节点/完整后代对照**。当前证据不批准降低ARI、改变停止尺寸或改写OG/scorer。
   搜索日程或拓扑目标变更仍需显式ADR修订，纯增加预算也须先说明具体待覆盖区间。

本轮只有审计工具与文档，不修改生产科学逻辑，不提交新Slurm任务。
本地现有首次分散测试+新增污染fixture测试：**4 passed**；
真实数据运行的hash、成员唯一映射、source包含关系、组成员与导出一致性断言全部通过。

固定同mean SSN的原始V1优先四RefOG对照已完成，见 [V1原始hierarchy结果](v1-original-fixed-mean-hierarchy-results.md)。
