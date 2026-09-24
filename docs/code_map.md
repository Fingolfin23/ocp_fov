# 代码与资料索引

本次整理保留原有脚本、CSV 和缓存的根目录路径，通过说明文档和示例建立清晰入口。现有导入和默认数据路径依赖这种布局；本次没有借整理目录改变核心算法、车辆模型或历史实验参数。

## 阅读顺序

1. [首页](../README.md)：通过图片和动图理解“角度极值—双侧命中比较—跨候选比较—可见域分割”。
2. [数学证明 PDF](proofs/geometric_visibility_proofs.pdf)：了解模型、假设及几何结论；[Markdown 源稿](proofs/geometric_visibility_proofs.md)可用于进一步修订。
3. `examples/inspect_geometry.py`：检查一个真实赛道观察位置的计算结果。
4. `fov_extended.py`：对应原有几何实现；再看 `fov_to_frenet_bounds.py` 的规划接口。
5. [实现说明](implementation_notes.md)：区分精确数学构造、采样实现和规划实验目前的状态。

## 原有代码

| 类别 | 文件 | 职责与入口 |
| --- | --- | --- |
| 核心几何 | `fov_extended.py` | 生成采样边界、检测方位角极值、比较同射线同侧／异侧后续命中、选择全局前沿、输出局部事件。便捷入口为 `compute_fov_at_s0`。 |
| 缓存生成 | `precompute_fov_cache.py` | 读取 CSV，为采样观察位置生成 `fov_cache.json` 和紧凑 NPZ；有命令行参数。 |
| Frenet 规划接口 | `fov_to_frenet_bounds.py` | 将缓存前沿与未知区域映射为横向边界，提供 CasADi Opti 硬／软约束接口。 |
| 单点可视化 | `plot_fov.py` | 选择一个缓存观察位置，展示切点、射线及前沿。 |
| 赛道图与动图 | `viz_fov.py` | 总览、跟随视角、GIF／MP4；通过 `--track` 显式指定轨迹数据。 |
| FOV 规划实验 | `OPT_LAPTIME_FULLFOV.py` | 使用本库轨迹和缓存的圈速优化实验；运行结果仍需结合实现说明解释。 |
| 圈速基线 | `OPT_LAPTIME.py` | 不接入完整 FOV 规划带的优化实验。 |
| 障碍物扩展 | `OCP.py`、`solve_min_time_with_fov.py` | 引用了本库未收录的 `obstacle_avoidance` 包和障碍物 CSV。 |
| Opti 版本 | `solve_min_time_opti_with_fov.py` | 通过导入前述求解器复用模型，因此也需要外部障碍物依赖。 |
| 轮胎模型实验 | `compare_tire_models.py` | 线性与 Pacejka 轮胎模型对比，有命令行参数；不是几何算法验证程序。 |

## 历史文件

| 文件 | 当前状态 |
| --- | --- |
| `OPT_LAPTIME_.py`、`OPT_LAPTIME copy 2.py` | 两份内容完全相同的历史版本，保留以便追溯。 |
| `OPT_LAPTIME_FULLFOV copy.py` | 与主版本的侧向力增益、求解迭代上限、制动参数不同，不应按重复文件直接删除。 |
| `fov(banned).py`、`compute_fov(banned).ipynb` | 早期 FOV 探索；不是推荐入口。 |
| `test.py` | 历史绘图探索，引用外部相对路径；不是自动化回归测试集。 |
| `solve_min_time_mpc_fov.py`、`viz_ocp_fov.py` | 空文件占位，尚未实现相应功能。 |

## 数据与原有输出

- `Monza.csv`、`Nuerburgring.csv`：中心线与道路宽度；列为 `# x_m,y_m,w_tr_right_m,w_tr_left_m`。
- `fov_cache.json`：原有 Monza 示例缓存，共 1,158 个观察位置，观察位置间距 5 m，前向查询窗口 200 m。
- `fov_cache.npz`：紧凑缓存；现有 Python 规划入口主要读取 JSON，未发现已接入的 NPZ 消费路径。
- `monza_930.png`、`fov_follow.gif`、`fov_follow_b.gif`：原有展示输出，保留原位置。

新说明图位于 `docs/figures/`，动图位于 `docs/animations/`，生成入口为 `scripts/generate_idea_visuals.py`。它们依据报告第 6.3 节及图 12 的几何构造重新计算说明场景，不作为原有赛道实验结果。

## 依赖分层

| 文件 | 用途 |
| --- | --- |
| `requirements.txt` | 核心几何、数据加载及现有可视化所需包 |
| `requirements-ocp.txt` | 在基础依赖上加入 CasADi；不包含缺失的外部障碍物模块与数据 |
| `requirements-docs.txt` | 基础依赖与说明图所需 Shapely |
| `docs/package.json` | 可选 PDF 重建依赖；不参与几何算法运行 |

Python 依赖暂未锁定版本：当前没有覆盖所有历史求解器的完整兼容性验证。GIF 导出使用 Pillow；原有 MP4 导出还需要 FFmpeg。
