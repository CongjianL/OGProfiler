# 固定 hits 的逐层 SSN 差分

本诊断不重新运行 search、Leiden、MCL 或 OG 提取。
作业 1410698 对应修复前源码，其结论保存在独立结果报告中；当前脚本用于验证修复后的生产算法。

## 输入与参考

- 固定输入：已有运行的 `search/hits.parquet` 与 `input/proteins.parquet`。
- 两套算法共享物种顺序和物种内 V2 protein_id 顺序，明确排除原 FASTA
  顺序变化这个混杂因素；因此不是原有 OrthoFinder 图文件的简单重跑比较。
- 同源 hit 的 identity 使用 `(species_id, original_id)` 与共享 protein_id，
  不把不同物种同名 original_id 合并。
- 参考从实际 OrthoFinder 3.1.5 安装目录读取，复制相关源码到新运行目录，
  记录其 SHA256。执行原始 scnorm/WaterfallMethod 类与 sparse_max_row 定义，
  排除搜索及调度模块导入。
- 比较默认全量 OF `v2_scores=False` 分支，不覆盖增量 assign。
- 固定 hits 和蛋白表在比较前后重新计算 SHA256。

## 比较层

1. **最大 bitscore 去重**：相同输入的最大分数矩阵与 V2 实际去重算子
   作用于 bitscore 的诊断矩阵比较；修复后 V2 在 NBS 前进行最大 HSP 去重。
2. **NBS 样本/参数**：每个有向物种对记录原始/唯一 hits 数、两套拟合
   样本数、样本多重集差异及 `(a,b)`，包括退化和 too-few-hit 分支。
3. **完整方向 B**：比较所有非零归一化分数，不提前做方向阈值过滤。
4. **BH/RBH**：逐物种对比较支持集合。RBH 的比较范围是跨物种；
   同物种 BH 保留，OF 同物种 reciprocal BH 不参与跨物种 cutoff。
5. **每个蛋白 cutoff**：逐蛋白比较实际路径的最低 RBH/无 RBH fallback。
6. **方向 connect**：V2 调用实际生产 selector；独立记录的 cutoff 还必须
   重建出完全相同的方向选择集合，否则审计报错。
7. **方向 W 与无向对象**：分清以下对象，而不是把 OF 两个方向覆盖成一个值。

完整差异数量与 bounded samples 分开统计；浮点比较默认
`abs_error <= 1e-8 + 1e-6 * max(abs(left),abs(right))`。
同时记录 exact difference count 和容差外数量。

## 必须保持区分的图对象

### OrthoFinder → MCL

```text
B：完整方向归一化分数
C：方向阈值是否通过，数值为 0 或 1
W = (C + C.T) ⊙ B
```

W 一般不对称，两方向都通过时连接倍数为 2。
writer 将方向权重格式化为三位小数，MCL 消费这个矩阵。
报告分别统计完整精度与三位小数后的图支持和组件规模。

### 当前 V2 → Leiden

每个 canonical `(u,v)` 只有一个权重 Uuv；当前默认 `mean` 是
完整方向 W 的 `(Wuv+Wvu)/2`。igraph 使用 `directed=False`。
这仍是 OF 方向 W 的投影，而非 MCL 的输入对象。

### 仅供诊断的对象

- **V2 directional same-assembly**：将 OF 的 W 公式应用到 V2 的 B/C，
  比较上游差异，明确标注为诊断矩阵，不声称这是 V2 生产输出。
- **V2 symmetric lift**：把一条无向边的权重填入两个方向，只用于
  数值比较；统计 direction 数时不会误当作两条无向边。
- **OF mean projection**：`(W+W.T)/2`。
- **OF max projection**：`max(W,W.T)`。

后两者与 V2 比较时只取上三角。修复后生产默认使用 mean 投影；
max 保留为诊断候选。两者都不等价于方向 MCL 输入。

## 固定 B 的控制组

另将 V2 的完整 B 同时送入两套下游规则，记录：

```text
BH_sameB_control
RBH_sameB_control
cutoff_sameB_control
connect_sameB_control
```

这样可以区分 NBS 差异的下游传播与规则自身的边界差异。
修复前近似并列 RBH 的重复索引赋值影响，体现在同 B 的 cutoff/control 中；
修复后的 cutoff 精确保留这一 OF 边界语义，预期 control 差异为零。

## 输出

- `report.json`：完整报告、阶段状态、图对象/投影定义、依赖和输入校验。
- `progress.json`：阶段进度，不是完成标志。
- `normalization_pairs.tsv`：各方向物种对拟合样本、参数与警告。
- `stage_summary.tsv`：按方向物种对的完整差异统计。
- `stage_totals.tsv`：按阶段汇总，计数相加、最大误差取 max。
- `difference_samples.tsv`：有限样例，含两个端点的物种和原始 ID。
- `reference/orthofinder/`：本次真正执行的参考源码副本。

成功须同时满足：Slurm/进程正常退出、`completed=true`、
`input_integrity_verified=true`。这仅表示比较完成，不表示 SSN 相同。

## 执行

本地小型验证需要 SciPy；修复后的生产项目也使用 SciPy curve_fit，并已声明该依赖。
测试使用 `OF_SOURCE_ROOT` 指向实际参考源码目录，没有参考源码时显式 skip。

全量数据通过项目标准 wrapper 提交：

```bash
./dev/slurm-submit slurm/fixed_hits_ssn_audit.sh INPUT_RUN OF_SOURCE_ROOT
```

作业资源：1 CPU、128 GB、12 小时。库线程固定为 1；不启动 job array。
输入根与参考目录从参数传入，输出独占 `DEV_RUN_DIR/ssn-audit`。
项目 wrapper 保存 dirty diff、源文件 manifest、源码 snapshot、RUN_ID 和 JOB_ID。

初始选择是修改后 Open_Orthobench `v2_modified_rep1` 的已有 hits，
因为它直接对应此前 precision 崩溃的运行。
此次工作先量化 SSN 差异，不据此直接归因最终 OG precision。
