# 当前采用方案：冻结 soft42

按用户决定，结束本轮选择层评分优化，采用job1411104/current-soft42作为当前
Embleya 13基因组方案。不晋升绝对图惩罚或组件强度零模型，不继续提交优化作业。

## 配置与出处

- 完整配置：`presets/embleya-soft42.yaml`，逐项复制已完成soft42运行的run.yaml。
- 原运行源码：`405ce3e`；RUN_ID `20261004T081659Z_405ce3ea83d6_0366c702_12341`。
- `hierarchy.topology_policy=soft_binary_24_v2`，`hierarchy.max_depth=42`。
- `orthogroups.strategy=v1_compatible`，refinement=false；其余参数见完整配置。
- 此配置含原运行的线程/worker值；后续资源调整需显式记录，不能冒充原始运行。
- 仓库通用DEFAULT_CONFIG仍为kway/depth20；使用本方案须显式加载该配置。
  这不是将Embleya结论推广为所有数据集的通用默认。

## 已接受的表现

冻结OF assigned参考的一致性指标：precision **95.7133%**，recall **69.0770%**，
pair F1 **80.2424%**，macro best member F1 **91.1391%**。
完整输出覆盖128483蛋白；primary评估122954蛋白，OF unassigned5529不作主指标负例。
这是与参考方法的一致性，不是独立生物学正确性证明。用户认为当前表现已满足需求。

后续实验全部位于benchmarks，未替换生产OG算法，故无需reset Git历史或删除源码。
保留原始结果、诊断、失败实验与PR溯源，现有远端结果不覆盖、不重跑。
若未来重启研发，应另定目标与验收标准；本轮的待办里程碑暂不继续执行。

配置已用当前仓库源码的load_config验证，解析后与冻结run.yaml完全一致。
本地.venv中已有安装包版本较旧，会拒绝admission_policy；运行时应使用
当前仓库源码（如设置PYTHONPATH=src）或匹配版本的安装包，避免调用旧包。
