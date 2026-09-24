# 几何可见域：展示图片与动图

这些图根据原报告第 6.3 节及图 12 的思想重新制作，展示“角度局部极值 → 同一射线双侧序号比较 → 全局候选比较 → 完整可见域分割”。所有切点、交点和盲区都来自下述实际几何的数值计算，没有手工放置虚构交点。

**范围：** 这是用于解释算法的光滑道路示例（illustrative geometry），不是 Monza/Nürburgring 的运行结果，也不是原始离散程序的性能或正确性认证。原报告图片未被重绘为新的实验数据。

## 展示文件

| 文件 | 用途 |
|---|---|
| [01_tangent_events.png](figures/01_tangent_events.png) / [SVG](figures/01_tangent_events.svg) | 空间道路和左右边界方位角曲线对应起来；标出通过前向、内侧判据的三个切向事件。 |
| [02_two_stage_selection.png](figures/02_two_stage_selection.png) / [SVG](figures/02_two_stage_selection.svg) | 在同一实际道路上展示两层比较及数值表格，是推荐放在仓库首页的总览图。 |
| [03_complete_visibility.png](figures/03_complete_visibility.png) / [SVG](figures/03_complete_visibility.svg) | 对比只保留整体前沿与保留完整分割；说明最远可见位置不能替代局部盲区信息。 |
| [01_boundary_angle_scan.gif](animations/01_boundary_angle_scan.gif) | 同时沿两边界的共同道路序号扫描，空间视线与角度曲线同步变化。 |
| [02_two_stage_story.gif](animations/02_two_stage_story.gif) | 六步解释三个事件的分类、候选选择及最终局部盲区保留。 |

PNG 为 240 dpi；SVG 为矢量图。英文标签方便后续放进英文论文；可调整生成脚本的颜色、文字和尺寸。GIF 用于 README / 演示，不用于表达计算性能。图中的“two levels”指两个组织层次，不是整个算法只有两次标量比较。

## 实际使用的道路

中心线是半径 `R = 10 m`、共同弧长 `s ∈ [0,18] m` 的圆弧：

```text
c(s) = (R sin(s/R), R - R cos(s/R))
N(s) = (-sin(s/R), cos(s/R))
B_L(s) = c(s) + w_L(s) N(s)
B_R(s) = c(s) - w_R(s) N(s)
P = c(0) = (0,0)
```

宽度为正的高斯函数叠加，精确定义位于 [`scripts/geometry_model.py`](../scripts/geometry_model.py)。全窗口转角小于 π，所有半径为正，各道路横截面是不同极角的径向线段，因此该几何具有全局单射的 Frenet 映射；不存在道路自交或内洞。视线不允许穿越道路外部，没有独立内部遮挡物。

| 事件 | 切点序号 | 同侧首个后续命中 | 异侧首个后续命中 | 第一层分类 |
|---|---:|---:|---:|---|
| T1 | 1.896534 | 4.875193 | 13.084489 | 同侧盲区 |
| T2 | 5.996138 | 不存在于窗口内 | 10.306371 | 整体前沿候选 |
| T3 | 10.848432 | 不存在于窗口内 | 14.954958 | 整体前沿候选 |

因此第二层选择 `min(10.306371, 14.954958) = 10.306371 m`。T1 的同侧命中值不进入这次最小值比较。T3 位于更远处，保留为一个冗余前沿候选，用于说明扩大候选集不必改变已经找到的首个前沿。图中没有把每个事件都假设为从 P 可见。

完整可见区域由 T2 的近侧区域扣除 T1 产生的实际边界弧盲区得到。盲区有面积，边界是实际道路弧与分割线；不是三点共线的 `P–T–U` 三角形。

## 可复现与核验

依赖：Python 3、NumPy、SciPy、Matplotlib、Shapely、Pillow。

从仓库根目录运行：

```bash
python3 scripts/generate_idea_visuals.py --out docs/figures --animations docs/animations
```

运行时先进行几何核验，通过后再输出图片与动图。核验结果记录于 [`geometry_checks.json`](figures/geometry_checks.json)：

- 检查道路多边形有效性、切向与对齐残差；
- 检查每个分割线的内部点处于道路内，且道路序号严格向前增加；
- 固定随机种子，在 4,000 个道路点上比较事件分割与独立的整条视线径向净空判据；后者在视线参数上密集括区间，再对局部极小值数值精化。

这些数值检查是展示几何的一致性核对，不能替代数学证明；数值根枚举也没有被包装为形式化的完备性证书。该脚本不调用原库中的 OCP 或 FOV 缓存，不证明原库所有阈值与分支正确，不包含实车结果，也不声称速度或安全性能提升。
