# SAC-SQ+ R399 — 独立 witness 下同输出合同配对实验与科研裁决

**完成日期：2026-11-11**  
**主稿：R398（冻结，不在本轮修改）**  
**研究问题：**在同一冻结边区间成本盒、同一独立生成的 witness pool、同一全部路径身份与精确有理 β-support 合同下，R390 exact continuous witness-interval cover 相对 R357 16-cell 多见证排除，实际节省多少完整认证成本？

## 1. 科研结论

**R399 成功补齐论文最直接的内部同保证性能证据。** 原始真实 L/U 输入上的11组工作负载全部由两种独立 witness 初始化后的认证器返回 `COMPLETE`；经过统一 R354 精确支持过滤，每一次运行的路径身份、分段数、闭支持上下端点都与 R392 冻结完整输出精确一致。

- **11/11 workload arms，135组配对、270次认证器执行，全部 COMPLETE 且精确支持一致。** n256的8组每组15次；n1024的3组每组5次。所有重复都从同一冻结 L/U 输入与新生成的 witness pool 开始。
- **认证器＋精确支持 CPU：** 11 arm 的每臂中位数加速几何平均 **2.633×**（旧16-cell / 新 exact interval）。
- **加入共同 witness 初始化 CPU 后：** 11 arm 的每臂中位数加速几何平均 **1.650×**；这一口径覆盖 **51 fixed-scenario DP calls + midpoint exact-k DP + 认证 + R354 支持**，但不包括上游昂贵模型拟合/源证据取得。
- **标签创建量：** 11 arm 跨臂几何平均旧/新 **20.90×**；极端 R343 n1024 **212,256→161**，Occupancy SAC n1024 **194,516→161**。
- SECOM n256 strict 的完整认证中位 CPU **0.0595→0.0621s**，新认证器略慢；这是同合同完整输出已满足但小问题中 hull 构造成本占比较高的实例。该盒加入初始化后新方法同样略慢（比率0.987×）。**所有工作负载的效率改进应按具体情况和汇总口径陈述，不宣称逐盒严格领先。**

这是比 R390 早期 archive-seeded n256 配对更有直接投稿价值的一轮：现在完整输出、独立候选生成及新证书效率使用同一配对协议。它仍属于**自有两种同骨架算法的受控对照**，不能称为 Tropical Dijkstra 原作者的官方实现或独立外部同合同算法基线。

## 2. 实验设计及可复现方法

**输入：**R378 完成严格真实取证后归档的11个 `(L,U)` 独立边区间成本矩阵，`h=8`，`β∈[100000000000,1200000000000]`（1e12整数缩放）；8个 n256 arm、3个 n1024 arm。这11个 arm 包括同一数据来源的不同成本盒与控制器输出，**不是11个彼此独立数据集**。

**独立 witness：**仅读取 L/U。令 `M=L+floor((U-L)/2)`，在每个L/M/U场景对17个等距惩罚点进行固定成本最短路；令种子分段数区间为 `[kmin,kmax]`，再从 M 场景 exact-k DP 恢复该范围每种 k 的路径，去重为 witness pool。每盒独立重算5次，确认文本与 R392 的事后存档 seed 完全一致。最终路径 archive **不参与任何候选生成、C++认证器或 R354 的输入**。

**同骨架认证：**
1. `baseline16`：R378 `aug_tropical_multiwitness.cpp`，16个闭惩罚 cell，单 witness 在一个 cell 两端严格有效才能覆盖该 cell；保留 same-node strict dominance。
2. `exact_interval`：R390 `exact_interval_cover_hull.cpp`，count-indexed suffix lower hull + exact rational strict open witness intervals 并集覆盖；同样保留 same-node dominance。
3. 两者均设置 `max_created=500000`，时间上限120秒，均运行 **R354 `exact_best_scenario_support.cpp`** 并在 `COMPLETE` 后对每条候选计算精确闭 support。

**完全相同输出：**本轮对每次运行检查 `PATH` 身份、`K`、`LO`、`HI` 一一匹配 R392 `exact_support.txt`；支持过滤后的多余 sink 候选不算最终身份。结果检查只在运行结束后进行；与 R392 事后归档值比较不是候选生成时泄漏最终输出。

**CPU：**每个独立程序调用的子进程 user+system CPU via Python `resource.getrusage(RUSAGE_CHILDREN)` 差量；初始化取Python进程 `time.process_time()`。两认证器新进程交替顺序调用，统计各 arm 中位数；R354 同样独立子进程。完整成本定义为 `median(shared_seed + exact-k CPU) + median(enumeration + support CPU)`。**这个“完整”只指冻结成本盒之后的 witness+认证+支持四阶段，不包含生成 L/U 所需的 ridge-logistic 来源拟合。**

**运行环境：**Linux x86_64（Intel Xeon Platinum 8573C 容器），同进程宿主、相同 g++ `-O3 -std=c++17` 编译设置；NumPy 固定实现、OPENBLAS_NUM_THREADS/OMP_NUM_THREADS=1。CPU计时不跨机器混算。n256每盒15次，n1024每盒5次。并未以未来投稿期刊集群成本推断本地延迟。

## 3. 全11盒原始量化结果

表中所有CPU以秒为单位：`旧→新` 为 16-cell→exact-interval；“认证CPU”含C++ enumerator与独立R354 exact support，不含初始化；“全CPU”比率含共同初始化。路径数量为精确支持筛选之后的可能身份数。

| 工作负载 | 独立 seeds→witness pool | 创建labels 16-cell→exact | 支持路径 | 认证CPU旧→新 (s) | 认证加速 | 四阶段加速 |
|---|---:|---:|---:|---:|---:|---:|
| R343 Occupancy n1024 | 8→14 | 212,256→161 | 14 | 23.013→1.270 | 18.12× | 11.32× |
| R372 baseline n256 | 6→7 | 1,475→37 | 7 | 0.048→0.016 | 2.91× | 1.21× |
| R372 certificate-debt n256 | 6→7 | 2,769→37 | 7 | 0.079→0.015 | 5.13× | 1.42× |
| Occupancy n1024 strict | 8→14 | 11,516→161 | 14 | 2.516→1.254 | 2.01× | 1.57× |
| Occupancy n1024 SAC | 8→14 | 194,516→161 | 14 | 11.628→1.324 | 8.78× | 5.30× |
| EEG n256 strict | 9→15 | 242→136 | 15 | 0.050→0.034 | 1.47× | 1.09× |
| EEG n256 SAC | 9→15 | 462→136 | 15 | 0.058→0.033 | 1.74× | 1.13× |
| Occupancy n256 strict | 5→5 | 104→22 | 5 | 0.025→0.014 | 1.80× | 1.09× |
| Occupancy n256 SAC | 5→5 | 197→22 | 5 | 0.024→0.015 | 1.63× | 1.07× |
| SECOM n256 strict | 10→21 | 254→179 | 18 | 0.059→0.062 | 0.96× | 0.99× |
| SECOM n256 SAC | 10→21 | 483→179 | 18 | 0.092→0.074 | 1.23× | 1.08× |

**统计汇总：**11盒全部同支持输出；各臂CPU比率几何平均认证 2.633×、四阶段 1.650×。由于部分工作负载共享同一数据来源或同一任务下不同成本盒，汇总只作为**描述性 arm-level aggregate**，不当作独立样本显著性检验。

**分项CPU和逐次观测**详见 `R399_PAIRED_SUMMARY.csv`、`runs/R399_SUMMARY.json`、`runs/R399_RAW_RUNS.json`、每盒每轮 `.log` 与 `.support.txt`。独立候选文本在每盒 `independent_witnesses.txt`；全部程序源码见 `src/`。

## 4. 与历史结果的关系

R390的1.43–2.17× n256局部 CPU 对比使用了历史 archive witness；R392的11/11 COMPLETE 使用从输入生成的 candidate pool。R399将二者合并到统一同机协议，以更强的证据证明**新 exact continuous witness-interval certificate 在真实独立初始化时仍提供计算收益**。

R399 不是原ECS stagewise CPU的重新测试；原四任务4.53–10.48× CPU仍是旧在线取证 + 历史16-cell R357 + R354的分阶段计时，并非 R399 四阶段的同机时钟。R399 的11臂认证 2.633×和四阶段 1.650×不可与ECS倍率相乘；R343的18.1×认证倍率也不能当作所有输入的全流程加速。

一个重要的技术解释：即使最终可能路径集合相同，16-cell只有当某一 witness 在格点两端严格胜出时才能覆盖该cell；精确开区间可以在 cell 内切换 witness，并进行完整的有理覆盖。实际创建标签数因此下降，但对于小成本盒，精确 hull 的预处理也有固定开销，SECOM strict的局部负例体现了这一权衡。

## 5. R400论文采用建议与停止裁决

**建议立即进入 R400 结果整合，不再增加泛化数据集、更多 β 细网格或 SKAB 大预算。** R399在最优投入/产出比上完成了该阶段必需实验。论文主论证现可直接为：

1. 连续多 witness prefix cover 是核心可证数学机制；同β的可能路径身份与精确 support 构成新的完整输出合同。
2. L/M/U + exact-k 初始化独立于历史 archive；R392/R393 确认标准真实盒、同count指数路径、singleton transition的输出正确性。
3. **R399**在11个历史严格来源成本盒上，在独立witness、相同输出、同机配对中达到2.63×认证CPU以及1.65×四阶段CPU几何平均优势；最强实例与小盒例外一并在相应表中给出。
4. 原真实 GLM 取证/共享源证据仍作为外部参考的阶段价值。明确R399本地冻结盒与早期原ECS拟合流程的界面。

R398 原稿及五幅图在本轮未更改；R400整合时优先替换 RQ2实验表、更新摘要与结论一到两句。**保留成熟图形**，除非数值口径需要调整，优先用新表替换/补充旧archive CPU表，避免无意义重画五张图。没有必要开发另一个全新外部算法基线来重复证明完整区间合同。

## 6. 可直接复现

解压本研究包后，在具备 `g++`、Python3和NumPy的Linux环境中执行：

```bash
g++ -O3 -std=c++17 src/aug_tropical_multiwitness.cpp -o baseline16
g++ -O3 -std=c++17 src/exact_interval_cover_hull.cpp -o exact_interval
g++ -O3 -std=c++17 src/exact_best_scenario_support.cpp -o exact_support
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 run_r399.py
```

源文件会作为 `src/` 打包；`run_r399.py` 从 `r392/source/`导入独立候选生成函数，读取 `../R392-R393/real_boxes/inputs/`共享冻结二进制，仅在生成后对比 `r392/results/`保存的历史完整支持。旧实现与新实现共享相同认证输出和R354过滤程序。脚本将重现 `runs/` 中 JSON、逐次日志与总结。

## 7. 证据边界

- 同成本盒与同输出合同下的两种**内部认证机制**的公平配对，不是与第三方 Tropical Dijkstra 官方实现的对照。
- R399没有从原始UCI观测文件重新训练 ridge-logistic，也没有重跑整轮高精度源取证。
- 对同一数据来源下多个 arm 的描述性汇总不作统计独立性推断。
- CPU口径清楚区分局部认证、冻结盒四阶段与原ECS全阶段。
- 原先 SKAB n768 输出丰富、未能在预算内COMPLETE的情况未被本轮实验改变。

**最终科研裁决：** R399满足用户批准的实验合同，为论文补齐一条独立初始化、完整输出、同机配对的关键性能证据。下一步是 R400论文冻结与选刊投稿。
