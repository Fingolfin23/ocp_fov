# 代码与资料索引

当前代码库保留几何算法、规划接口、独立实验、可复现的展示材料和数学证明。重复备份、空占位、系统缓存及无消费路径的生成副本已清理；具有不同模型或实验参数的版本仍予以保留。

## 阅读顺序

1. [首页](../README.md)：通过图片和动图理解“角度极值—双侧命中比较—跨候选比较—可见域分割”。
2. [数学证明 PDF](proofs/geometric_visibility_proofs.pdf)：了解模型、假设及几何结论；[Markdown 源稿](proofs/geometric_visibility_proofs.md)可用于进一步修订。
3. `examples/inspect_geometry.py`：检查一个真实赛道观察位置的计算结果。
4. `fov_extended.py`：原有几何实现；再看 `fov_to_frenet_bounds.py` 的规划接口。
5. [实现说明](implementation_notes.md)：区分精确数学构造、采样实现和规划实验目前的状态。

## 算法与实验入口

| 类别 | 文件 | 职责与入口 |
| --- | --- | --- |
| 核心几何 | `fov_extended.py` | 生成采样边界、检测方位角极值、比较同射线同侧／异侧后续命中、选择全局前沿、输出局部事件。便捷入口为 `compute_fov_at_s0`。 |
| 缓存生成 | `precompute_fov_cache.py` | 读取 CSV，生成完整 JSON 与可选使用的紧凑 NPZ 输出；仓库只提交当前使用的 JSON。 |
| Frenet 规划接口 | `fov_to_frenet_bounds.py` | 将缓存前沿与未知区域映射为横向边界，提供 CasADi Opti 硬／软约束接口。 |
| 单点可视化 | `plot_fov.py` | 选择一个缓存观察位置，展示切点、射线及前沿。 |
| 原有赛道图与动图 | `viz_fov.py` | 总览、跟随视角、GIF／MP4；通过 `--track` 显式指定轨迹数据。 |
| Monza 全圈展示 | `scripts/generate_monza_replay.py`、`scripts/monza_replay_geometry.py` | 周期闭合赛道，全局地图、跟随镜头与前沿距离曲线；区分原采样算法输出与独立参考填色。 |
| 原理图与教学动图 | `scripts/generate_idea_visuals.py`、`scripts/geometry_model.py` | 依据报告中的几何构造重新计算示意场景，导出 PNG、SVG 和 GIF。 |
| FOV 规划实验 | `OPT_LAPTIME_FULLFOV.py` | 使用本库轨迹和缓存的圈速优化实验；运行结果仍需结合实现说明解释。 |
| 圈速基线 | `OPT_LAPTIME.py` | 不接入完整 FOV 规划带的优化实验。 |
| 障碍物扩展 | `OCP.py`、`solve_min_time_with_fov.py` | 不同的障碍物规划入口，引用本库未收录的 `obstacle_avoidance` 包和障碍物 CSV。 |
| Opti 版本 | `solve_min_time_opti_with_fov.py` | 导入前述求解器以复用模型，也需要外部障碍物依赖。 |
| 轮胎模型实验 | `compare_tire_models.py` | 线性与 Pacejka 轮胎模型对比，有命令行参数；不是几何算法验证程序。 |

## 保留的历史实验

| 文件 | 保留原因 |
| --- | --- |
| `OPT_LAPTIME_.py` | 与主圈速基线的轮胎表达式、载荷使用、赛道切片和迭代参数不同。与其完全相同的 `OPT_LAPTIME copy 2.py` 已移除。 |
| `OPT_LAPTIME_FULLFOV copy.py` | 与主版本的侧向力增益、求解迭代上限、制动参数不同，属于独立参数实验。 |
| `compute_fov(banned).ipynb` | 早期中心线射线与滤波探索，有独立源码；已清除输出图像和执行计数。它不是当前可直接运行的算法入口。 |

未完成的 `fov(banned).py`、重复绘图片段 `test.py` 和两个空文件 `solve_min_time_mpc_fov.py`、`viz_ocp_fov.py` 已移除。其旧内容仍可从清理前的 [Git 提交 bfc877f](https://github.com/Fingolfin23/ocp_fov/tree/bfc877f4cdf556896dc1bae3ae2dd884b0817ad5) 查看。此次清理未改变现有几何函数或求解器参数。

## 数据与展示材料

- `Monza.csv`、`Nuerburgring.csv`：中心线与道路宽度；列为 `# x_m,y_m,w_tr_right_m,w_tr_left_m`。
- `fov_cache.json`：原有 Monza 示例缓存，共 1,158 个观察位置，间距 5 m，前向查询窗口 200 m。
- `fov_cache.npz`：可由缓存生成器重建；与 JSON 对应数组相同且没有当前读取入口，因此不再入库。
- `fov_follow.gif`、`fov_follow_b.gif`：两份原有展示，保留原位置并在 README 中直接显示。
- `docs/figures/`：首页 PNG、论文用 SVG、Monza 海报及数值检查记录。
- `docs/animations/`：教学动图与 Monza 全圈 GIF／高清 MP4。
- `docs/proofs/`：英文证明 PDF、可编辑源稿、三份扩展推导及解析例子的验证材料。

旧单点截图 `monza_930.png` 已移除，可用 `plot_fov.py` 从保留的 JSON 重建。证明构建所产生的 HTML 和排版日志也只保留为本地生成文件。重建方法见 [构建说明](BUILDING.md)。

## 依赖分层

| 文件 | 用途 |
| --- | --- |
| `requirements.txt` | 核心几何、数据加载及现有可视化所需包 |
| `requirements-ocp.txt` | 在基础依赖上加入 CasADi；不包含缺失的外部障碍物模块与数据 |
| `requirements-docs.txt` | 基础依赖、几何参考检查所需 Shapely，以及全圈视频编码依赖 |
| `docs/package.json` | 可选 PDF 重建依赖；不参与几何算法运行 |

当前没有覆盖所有历史求解器的统一版本锁和完整兼容性验证。原有 `viz_fov.py` 的 MP4 导出需要系统 FFmpeg；新全圈脚本使用 `imageio-ffmpeg` 提供的编码器。
