# 固定完整树切分：组件内强度零模型实验

算法身份 `fixed-tree-component-strength-null-v1`，只在benchmarks内实现，生产默认不变。
这是图目标对照，不是祖先OG模型；不以OF标签决定候选或评分。

## 在看本轮成绩前固定的目标

对一个SSN组件，m为所有无向边权之和，s_i为蛋白i在该组件内的加权度。
一个候选完整子树v的内部无向边权为w_v，节点强度K_v=Σ(i∈v)s_i。

```
keep(v) = w_v/m - (K_v/(2m))^2
best(v) = max(keep(v), sum(best(child)))
```

分辨率固定为1，无阈值、无搜索。零模型项包含平方中的自项；这正是本次固定定义，
不在实现中删除自项或根据子群数调整。m=0时所有候选分数为0、同分保留父节点，
明确记录零权组件；singleton保持覆盖。任意UNRESOLVED仍阻断。

相较绝对pair_penalty，它以组件内实际节点强度定义预期连接权重，均匀缩放全部
边权不改变目标；不需要子群内部密度，因此单点子群也有定义。
“局部”指各SSN组件独立归一化，不代表每个家族自适应校准。巨大component0仍可能
有尺度效应；本轮不能预先保证解决该问题，也不改变SSN、层级或Leiden目标。

完整切分互不嵌套，全蛋白覆盖。输出保留算法身份、resolution=1、pair_penalty=null、
组件总边权与节点强度；不得将这个1误称为旧pair_penalty=1或重用旧参数意义。
批量CLI使用互斥开关 `--strength-null`，原 `--penalties` 路径不变。

## 预定评估

冻结job1411104/current-soft42与同次OF参考。单个新候选，先预测后评估，
对比已冻结的soft42与job1411414的pair_penalty=0.1。不筛选新参数。
报告pair/B-cubed/macro、同/跨物种分层、误合并拆分及根/全叶退化。
仍为已暴露的13基因组探索，不是保留集验证；结果差也完整保留。

脚本 `slurm/embleya_strength_null.sh`：2CPU/8GB/2h，无数组，先fixture smoke。
新源码快照、新目录、不覆盖既有结果。数学fixture覆盖手算两团弱桥目标、
均匀缩放不变、全单点子群、零边、参考扰动不改变预测及新旧参数互斥。
