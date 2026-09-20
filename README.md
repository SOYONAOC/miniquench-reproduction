# Mini-quenching 复现审计（独立 CPU 项目）

固定论文：Judah Luberto & Steven R. Furlanetto, *A model to mini-quench early galaxies by balancing stellar feedback and gas accretion*, [arXiv:2609.19265v1](https://arxiv.org/abs/2609.19265v1)。

**本轮已实际实现、运行和测试，但没有通过完整论文复现验收。** 图 1 的六条参数化 IMF 曲线与原图吻合；图 4 是定量有偏差的近似；图 2/5/8 的完整反馈复现受热学闭合和启动边界问题阻碍。`sn_only` 是独立调试对照，不能称为完整模型。先读 [中文报告](reproduction_report.md)，审阅版幻灯片位于 `slides/reproduction_review.pdf`：16 页（含图 2 偏差诊断），所有比较页均为左侧原图、右侧复现，多面板逐项拆页。

## 外部审阅入口

先读 [给 GPT Pro 的阅读入口](REVIEW_START_HERE.md)，再看 [图 2 最新失败诊断](docs/figure2_failure_audit.md)。对照 PDF 均为左侧原图、右侧计算。仓库包含真实结果及失败记录；虚拟环境、可重建缓存和第三方 Git 历史不上传。

## 已测试的运行方式

在 `miniquench_reproduction/` 目录执行；使用自己的 `.venv`，不导入 AuroraLF。下面这些命令已经实际执行。第一次安装使用 `uv venv --seed miniquench_reproduction/.venv`（从父工作目录）；实际成功的版本记录在 `requirements.lock.txt` 和 `outputs/logs/environment.txt`。

```bash
.venv/bin/python -m pip install -e ".[dev]"
.venv/bin/python -m pytest -q
.venv/bin/python scripts/prepare_sources.py
.venv/bin/python scripts/extract_reference.py
.venv/bin/python scripts/reproduce.py --figures 2 --profile quick
.venv/bin/python scripts/reproduce.py --figures 1 2 4 5 8 --profile full --workers 3
.venv/bin/python scripts/reproduce.py --figures 2 5 8 --profile full --model sn_only --workers 3
.venv/bin/python scripts/validate.py
.venv/bin/python scripts/audit_figure4.py
.venv/bin/python scripts/compare_paper.py
.venv/bin/python scripts/render_results.py
.venv/bin/python scripts/prepare_comparison_assets.py
```

源文件已随项目保存在 `external_data/`，Nebrin 原始代码在 `third_party/`。准备脚本仅从这些真实资料提取表格和校验值，不联网生成替代资料。图形使用随项目保存的 ApJ 样式，字体明确为 Matplotlib 自带 DejaVu Serif，避免依赖宿主未安装的 Times New Roman。

`quick` 使用 9 个质量点、0.2 Myr 历史格；`full` 使用 31 个质量点、0.05 Myr 历史格，并在可判定的启动临界区加密。最大步长和输出间隔另行配置。最多 4 个工作进程，每进程 BLAS 单线程。`reproduce.py` 默认 `full_approx`，默认回落分支是 `until_return`（见下）；输出 CSV/JSON 包含失败与未返回状态，命令完成不意味着科学验收通过。

## 两层约定与边界

- 直接物理输入：论文自己的宇宙学、等温球、增长率、IMF、寿命、SN 能量及 Compton 冷却。`fstar`、观测映射效率、单次形成质量均独立。
- 未被作者唯一指定的部分：恒星半径、冷却壳层厚度/密度/中性度、辐射源空间分布等。`full_approx` 包含 SN、UV、Lyα 三项，调用 Nebrin 的真实拟合和 Schure 的真实冷却表，但这些闭合属于明确标记的推断。
- 低层 `Config` 默认 `loading='outward_only'`，对应只在 v>0 时加载的字面解释；CLI 扫描默认 `--loading until_return`，作为数值上较容易审计的独立分支。该选择未被确认为作者实现，两者的敏感性结果均保留。
- SFR 只在实际可解析的启动后关闭，100 Myr 是启动观察窗口。形成历史与自适应积分的 RHS 分离；关闭 SFR 后老恒星继续发光和产生 SN。
- `completed_cycle`：确实返回移动 R0；`not_launched_within_window`：未启动；`launched_not_returned_by_tmax`：只有时长下限；`numerical_failure`：显式错误，缺值不填零。
- `validate.py` 的收敛通过仅适用于它列明的代表轨道和 SN-only 临界质量，不涵盖完整反馈的失败点。

## 输出

- `docs/equation_map.md`、`parameters.md`、`ambiguities.md`：出处、参数与未决约定。
- `data/manifest.json`：原始论文、源包、外部冷却/反馈资料及 SHA256。
- `data/reference/`：作者矢量 PDF 的曲线路径坐标及轴标定；并非作者求解器数组。
- `configs/`：各图和实际扫描配置。
- `data/curves/summary_full_*.json`、`scan_full_*.csv`：所有运行状态、事件、配置和耗时。
- `data/curves/fig2_*_full_*.csv`：半径、速度、压力、SFR、质量、各力、反馈与温度的完整时间序列。
- `data/cache/`：按物理配置和代码指纹隔离的核/轨道缓存；允许断点续跑，不混用旧实现。
- `data/convergence.csv`、`threshold_convergence.json`、`sensitivity.json`、`paper_comparison.json`：数值误差和论文残差分别存储。
- `outputs/figures/`：当前 PNG/PDF；`outputs/logs/`：实际运行日志。

最初复现阶段没有提交或推送；之后按用户明确要求整理为独立 GitHub 仓库供外部审阅。没有申请 GPU、集群或付费资源，没有接入 AuroraLF，也没有下载模拟快照。其余论文图的状态逐项列在报告中。
