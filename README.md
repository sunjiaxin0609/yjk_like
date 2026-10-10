# yjk_like

对标盈建科（YJK）的结构计算设计软件：**C++17 内核 + Qt6 桌面端**，
零外部数值库依赖（仅标准库），自实现稀疏直接法、预条件 CG 与广义特征值求解。

> 算法细节见 [`NOTES.md`](NOTES.md)（数学层 / 单元层 / 模型层 / 静力 / 模态 /
> 反应谱 / 组合 / 后处理 / 混凝土验算的公式与流程逐条对照代码）。

---

## 功能总览

| 层 | 内容 | 状态 |
| --- | --- | --- |
| 数学 | 稀疏 LDLᵀ 直接法（Gilbert–Peierls）· 预条件 CG · 子空间迭代特征值 | ✅ |
| 单元 | 3D 梁（Timoshenko 12 DOF，端部释放/铰）· 壳（膜 + MITC4 板弯曲）· 弹簧 | ✅ |
| 模型 | 材料/截面库 · 轴网建模（柱/梁/板/墙自动生成）· 刚性楼板 MPC · 45° 导荷 · 剪力墙 | ✅ |
| 分析 | 线性静力 · 模态（振型数自动）· 反应谱（GB 50011，SRSS/CQC）· 荷载组合与包络 | ✅ |
| 验算 | 混凝土正截面承载力（受弯/受压/受扭，GB 50010 第 6 章） | ✅ P4-T1 |
| 界面 | 命令行 yjk · HTML 后处理界面 · Qt6 桌面端（交互建模、撤销/重做、另存 .yjk） | ✅ |

内核 22 个测试目标 / 870+ 断言，GCC 16.2 与 GCC 13.1 双编译器全绿（2026-10-10）。

---

## 快速上手

### 命令行

```bash
yjk run examples/frame3x4.yjk -o out/frame3x4
```

输出文件（静态分析 3 个；脚本含动态命令时多出 2 个）：

| 文件 | 内容 |
| --- | --- |
| `out/frame3x4.html` | 浏览器直接打开：三维变形动画、应力云图、内力图、楼层表（离线单文件） |
| `out/frame3x4.json` | 结果数据（自带 `_schema` 说明字段含义） |
| `out/frame3x4.txt` | 文本报告：位移包络、楼层指标、构件应力 TOP 10、平衡校核 |
| `out/frame3x4.dynamic.json` | 动态分析数据：模态/反应谱/组合内力（含模态命令时） |
| `out/frame3x4.dynamic.txt` | 动态分析文本报告：频率/周期/参与质量、逐层剪重比、组合控制内力 |

动态分析（模态 + 反应谱 + 荷载组合）示例：

```bash
yjk run examples/frame3x4_dynamic.yjk -o out/frame3x4_dynamic
```

### 桌面端

```bash
scripts/build-gui.sh --run     # 构建并启动
```

界面里打开 `.yjk` 脚本，按 `F5` 即开始计算。

---

## 构建与测试

### 内核（无外部依赖）

```bash
cmake -B build -G Ninja -DCMAKE_MAKE_PROGRAM=<ninja路径> \
      -DCMAKE_CXX_COMPILER=<g++路径> -DYJK_BUILD_TESTS=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

需要 CMake ≥ 3.16 与 C++17 编译器。常用选项：
`YJK_BUILD_TESTS`（默认 ON）、`YJK_BUILD_GUI`（默认 OFF）、`YJK_BUILD_SHARED`、`YJK_MSVC`。

### Qt6 桌面端

```bash
scripts/build-gui.sh            # 配置 + 编译 + 测试
scripts/build-gui.sh --deploy   # 打可独立拷贝的 dist/（已收齐运行库，拷走即可运行）
scripts/build-gui.sh --run      # 构建后直接启动
```

脚本自动使用 **Qt 自带的 MinGW** 编译 —— Qt 与 MinGW 版本必须配套
（6.8.3 对应 GCC 13.1.0），用别的编译器会"能编过、一运行就崩"。
首次安装 Qt（约 300 MB）：

```bash
pip install aqtinstall
python -m aqt install-qt   windows desktop 6.8.3 win64_mingw \
       -O C:/toolchain/Qt --archives qtbase
python -m aqt install-tool windows desktop tools_mingw1310 \
       qt.tools.win64_mingw1310 -O C:/toolchain/Qt
```

Qt 装在别处时用环境变量指定：`QTDIR=... QTMINGW=... scripts/build-gui.sh`。

---

## 目录结构

```
include/yjk/math/       数学层：SparseMatrix / Solver / Ordering / EigenSolver
include/yjk/element/    单元层：BeamElement3D / ShellElement4 / SpringElement3D
include/yjk/model/      模型层：Section / Material / Model（MPC）/ GridMesh（轴网·导荷）
include/yjk/analysis/   分析层：StaticAnalysis / ModalAnalysis / ResponseSpectrum / LoadCombination
include/yjk/post/       后处理层：PostProcessor / ResultExport / ViewerTemplate
include/yjk/check/      规范验算：RC（混凝土正截面，GB 50010）
include/yjk/interact/   交互建模纯函数：BuildOps / PickCore / CommandStack
include/yjk/io/         建模脚本：ModelScript（.yjk）
src/                    对应实现
apps/yjk.cpp            命令行入口（run / check）
apps/ui/                Qt 桌面端
scripts/build-gui.sh    Qt 构建脚本
tests/                  内核测试目标
examples/               frame3x4.yjk · frame3x4_dynamic.yjk · frame_wall.yjk 等
```

---

## 文档与验证

- [`NOTES.md`](NOTES.md)：算法与计算过程详解（公式、参数、流程逐条对照代码）；
  已验证精度锚点（悬臂梁/简支板/模态频率/剪重比等对解析解）见其 §15。
- 设计约定：内部统一 `kN-m-N`（长度 m、力 kN、应力 kPa）；内核不依赖 Qt；
  绝不静默出错（奇异/不收敛/精度不足必须给出可读诊断）；未验证的代码等于零
  —— 每个单元/算法都必须对标解析解或独立实现。

## 已知限制

| 限制 | 影响 |
| --- | --- |
| 无实体单元、墙元（平面内刚域） | 实体段只能靠细分的壳元近似 |
| 刚性楼板假定仅单塔 | 主节点每层 3 个自由度，多塔/大底盘需分块 |
| 无自动重编号 | 数万自由度求解慢；依赖"按楼层编号"的建模约定 |
| P-Δ · 温度 · 施工模拟 | 二阶效应 / 温度工况 / 竖向变形差不可用 |
| 暂未第三方交叉验证 | 验证均为简单构件对标解析解 |
| 导入导出仅 `.yjk` | YJK/PKPM、DXF、IFC 未接入 |

## 开发路线

- 🔵 **P4 规范验算与配筋**（进行中）：斜截面 / 裂缝 / 挠度 / 实配钢筋 / 钢结构 / 抗震构造
- 🔵 **P5 工程化**：第三方交叉验证 · 求解器性能（AMD 重编号）· 2D 平面图 · 导入导出

---

## 使用说明

### 命令行 yjk

```
yjk run   <模型.yjk> [-o 输出前缀]    建模 → 校核 → 求解 → 后处理 → 导出
yjk check <模型.yjk>                  只建模与校核，不求解
yjk help
```

| 退出码 | 含义 |
| --- | --- |
| 0 | 成功 |
| 1 | 脚本 / 建模错误 |
| 2 | 校核未通过 |
| 3 | 求解失败 |
| 4 | 用法错误 |

输出前缀缺省为脚本内 `out` 命令指定的值，否则为 `yjk_out`。

### `.yjk` 脚本速查

建模命令（完整命令集见 `src/io/ModelScript.cpp`）：

| 命令 | 作用 |
| --- | --- |
| `material C30 concrete 30` | 定义材料（名 / 类型 / 强度） |
| `section COL rect 0.5 0.5` | 定义截面（柱 500×500） |
| `grid.axisX / axisY / story` | 轴网：X/Y 轴线坐标与楼面标高 |
| `grid.column / beam / slabthick` | 沿轴网生成柱/梁、板厚 |
| `grid.divX / divY` | 跨内剖分数（梁板共用、节点重合） |
| `grid.slab off` · `slabload.to beam` | 不建板 + 楼面荷载 45° 线导荷（真实工程常规做法） |
| `grid.wallx / wally / wallthick` | 沿轴线生成剪力墙 |
| `support base fixed` | 支座（`support base/node/story <目标号> <fixed\|pinned\|roller>`，见 §14.2） |
| `selfweight on` | 梁柱自重 |
| `slabload 4` | 楼面均布荷载 kPa（正数向下） |
| `nodeload top 15 0 0` | 节点荷载（位置 / Fx / Fy / Fz） |
| `out frame3x4` | 输出前缀 |

动态分析命令：

| 命令 | 作用 |
| --- | --- |
| `case D` | 切换后续荷载所属工况（缺省 D） |
| `modal auto` | 模态分析，振型数 `max(9, 层数×3)` |
| `spectrum 0.16 0.40 cqc` | 反应谱：αmax、Tg、组合方式（第 4 参为竖向系数 k，默认 0.65） |
| `gmass D L 0.5` | 重力代表值（恒 + ψ·活）转节点质量 |
| `combo basicvar B1 D L W` | 荷载组合（standard / basicvar / basicperm / seismic） |

### Qt 桌面端

| 操作 | 效果 |
| --- | --- |
| 左键拖动 | 旋转（轨道相机） |
| 中键 / Shift+左键 拖动 | 平移 |
| 滚轮 | 以光标为锚点缩放 |
| 左键单击 | 选中构件（属性面板 + 模型树同步） |
| `1` `2` `3` `4` | 轴测 / 正立面 / 侧立面 / 俯视 |
| `Home` / `F` | 缩放到全览 |
| `O` | 正交 / 透视切换 |
| 空格 | 变形动画开 / 关 |
| `F5` | 开始计算 |
| `Ctrl+Z` / `Ctrl+Y` | 撤销 / 重做 |

界面交互建模（拾取、增删构件、属性编辑、另存为 `.yjk`）与脚本走同一条内核链路，
导出后逐字段一致。

无头渲染（界面层自检 / 出图）：

```bash
yjk_ui --script examples/frame3x4.yjk --shot out.png \
       --preset iso --contour disp --window
```