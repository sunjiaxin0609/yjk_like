# ============================================================
#  yjk_like  构建配置
#
#  依赖：
#    · C++17 编译器（MSVC / MinGW-w64 / Clang）
#    · CMake ≥ 3.16
#    · Qt6      —— 仅 GUI 需要（选项 YJK_BUILD_GUI，默认 OFF）
#
#  快速开始：
#    cmake -B build -DCMAKE_BUILD_TYPE=Release
#    cmake --build build
#    ctest --test-dir build --output-on-failure
# ============================================================

cmake_minimum_required(VERSION 3.16)
project(yjk_like VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---- 选项 ----
option(YJK_BUILD_TESTS "构建单元测试" ON)
option(YJK_BUILD_GUI   "构建 Qt 界面（需要 Qt6）" OFF)
option(YJK_BUILD_SHARED"构建共享库"       OFF)
option(YJK_ASAN        "启用地址消毒器（调试用）" OFF)
option(YJK_MSVC        "生成 Visual Studio 工程（Windows）" ON)

if(YJK_ASAN)
  add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
  add_link_options(-fsanitize=address,undefined)
endif()

if(MSVC)
  add_compile_options(/W4 /permissive- /utf-8 /Zc:__cplusplus)
  add_compile_definitions(_USE_MATH_DEFINES NOMINMAX WIN32_LEAN_AND_MEAN)
  # 结构分析用 double 精度即可，但 /fp:precise 保证跨平台一致性
  add_compile_options(/fp:precise)
else()
  add_compile_options(-Wall -Wextra -Wpedantic -Wshadow)
endif()

# ============================================================
#  yjk_core —— 有限元内核（零外部依赖，纯标准库）
#  界面层（Qt）只依赖这一层，不反向依赖。
# ============================================================
set(YJK_CORE_SOURCES
  src/math/Math.cpp
  src/math/EigenSolver.cpp
  src/model/Model.cpp
  src/model/GridMesh.cpp
src/analysis/StaticAnalysis.cpp
  src/analysis/ModalAnalysis.cpp
  src/analysis/ResponseSpectrum.cpp
  src/analysis/LoadCombination.cpp
  src/post/PostProcessor.cpp
  src/post/ResultExport.cpp
src/io/ModelScript.cpp
  src/io/YjkWriter.cpp
  src/interact/PickCore.cpp
  src/interact/BuildOps.cpp
)

add_library(yjk_core STATIC ${YJK_CORE_SOURCES})
target_include_directories(yjk_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
set_target_properties(yjk_core PROPERTIES POSITION_INDEPENDENT_CODE ON)

# ============================================================
#  测试
# ============================================================
if(YJK_BUILD_TESTS)
  enable_testing()

# 数学层：所有稀疏求解器必须与解析解/稠密解严格一致
  add_executable(test_math tests/test_math.cpp)
  target_link_libraries(test_math PRIVATE yjk_core)
  add_test(NAME math COMMAND test_math)

  # 广义特征值：Kφ=λMφ 与解析解（弹簧链闭式解）、稠密雅可比互验、残差判据
  add_executable(test_eigen tests/test_eigen.cpp)
  target_link_libraries(test_eigen PRIVATE yjk_core)
  add_test(NAME eigen COMMAND test_eigen)

  # 截面与材料：几何性质必须对标解析解与数值积分
  #   这一层是单元层的物理输入，算错不报错但内力全错
  add_executable(test_section tests/test_section.cpp)
  target_link_libraries(test_section PRIVATE yjk_core)
  add_test(NAME section COMMAND test_section)

  # 三维梁单元：局部坐标、刚度解析解、等效节点荷载、内力恢复、整体结构
  add_executable(test_beam tests/test_beam.cpp)
  target_link_libraries(test_beam PRIVATE yjk_core)
  add_test(NAME beam COMMAND test_beam)

  # 四节点壳单元：几何、膜（面内）已验证；板（面外）待 MITC4
  add_executable(test_shell tests/test_shell.cpp)
  target_link_libraries(test_shell PRIVATE yjk_core)
  add_test(NAME shell COMMAND test_shell)

  # 单元层：其他单元类型
  add_executable(test_element tests/test_element.cpp)
  target_link_libraries(test_element PRIVATE yjk_core)
  add_test(NAME element COMMAND test_element)

  # 模型层 + 分析算法：整体结构静力/模态
  add_executable(test_analysis tests/test_analysis.cpp)
  target_link_libraries(test_analysis PRIVATE yjk_core)
  add_test(NAME analysis COMMAND test_analysis)

  # 模型层端到端：节点/单元/荷载/约束 → 自由度编号 → 静力求解 → 内力与反力
  #   判据全部对标解析解（悬臂梁 PL³/3EI、简支梁反力 wL/2、门式刚架残差）
  add_executable(test_model tests/test_model.cpp)
  target_link_libraries(test_model PRIVATE yjk_core)
  add_test(NAME model COMMAND test_model)

  # 轴网建模：GridMesh 自动生成节点/柱/梁/板
  #   最关键的一条判据：板边节点必须与梁节点重合（否则荷载传不到梁上，且不报错）
  add_executable(test_grid tests/test_grid.cpp)
  target_link_libraries(test_grid PRIVATE yjk_core)
  add_test(NAME grid COMMAND test_grid)

  # 后处理 + 导出 + 建模脚本
  #   判据：内力图逐点对照解析解（含跨中抛物线与符号约定）、
  #         应力 σ=N/A±M·h/I 可手算复核、剪重比分母含面荷载、
  #         导出的 JSON 可被解析（NaN 会让界面白屏）
  add_executable(test_post tests/test_post.cpp)
  target_link_libraries(test_post PRIVATE yjk_core)
  add_test(NAME post COMMAND test_post)

# 刚性楼板（MPC）+ 板荷载导荷 —— "不建板也能算"这条路
  add_executable(test_floor tests/test_floor.cpp)
  target_link_libraries(test_floor PRIVATE yjk_core)
  add_test(NAME floor COMMAND test_floor)

  # 模态分析：ω/T/振型对标解析解（单自由度、两自由度链）、质量归一、刚性楼板
  #   质量折减与 MPC 刚度同源展开表；X 平动与 RZ 扭转简并时从属位移仍精确
  #   满足 ux_i = UX − Θ·dy_i —— 这是后续反应谱（T4）参与质量校验的基石
  add_executable(test_modal tests/test_modal.cpp)
  target_link_libraries(test_modal PRIVATE yjk_core)
  add_test(NAME modal COMMAND test_modal)

  # 设计反应谱与振型分解反应谱法：α(T) 四段曲线解析值、单自由度
  #   V = α·m*·g、两自由度逐模态 V_j = α_j·γ_j²·g 与 SRSS/CQC 组合、
  #   参与质量 ≥ 90% 检查（含方向无可动质量豁免与 checkMass 开关）
  add_executable(test_spectrum tests/test_spectrum.cpp)
  target_link_libraries(test_spectrum PRIVATE yjk_core)
  add_test(NAME spectrum COMMAND test_spectrum)

  # 荷载组合：标准/基本/地震组合系数与逐分量线性组合
  #   判据：用例工况结果手算组合值逐点对照（位移/内力/层剪力/包络极值）
  add_executable(test_combination tests/test_combination.cpp)
  target_link_libraries(test_combination PRIVATE yjk_core)
  add_test(NAME combination COMMAND test_combination)

  # P1 抗震精度补全：恒+0.5活转质量（GB 50011 5.1.3，Σweight/g + 自重不重复）、
  #   modal 自动振型数 max(9, 层数×3)、spectrum cqc 第三参透传 + 竖向系数
  #   线性比例锚点 + 手写 ρ_jk 交叉项 √(ΣΣρV_jV_k)、两层网格层间位移角/剪重比
  add_executable(test_p1_seismic tests/test_p1_seismic.cpp)
  target_link_libraries(test_p1_seismic PRIVATE yjk_core)
  add_test(NAME p1_seismic COMMAND test_p1_seismic)

  # 端部释放（铰）+ 弹簧单元：新 P1 特性
  #   test_hinge —— 静力凝聚解析解（简支含剪切 PL³/48EI·(1+Φ)、悬臂、端弯矩为 0）、
  #                 刚体模式、45° 倾斜、绕轴旋转（弯曲平面随 up 旋转 → 弱轴弯曲）
  #   test_spring—— 轴向/横向串联等效刚度、并联分配、45° 投影刚度、
  #                 k=0 释放语义（无荷载自动锁定、带荷载拦截）、平动刚体零内力、局部对角线
  add_executable(test_hinge tests/test_hinge.cpp)
  target_link_libraries(test_hinge PRIVATE yjk_core)
  add_test(NAME hinge COMMAND test_hinge)
  add_executable(test_spring tests/test_spring.cpp)
  target_link_libraries(test_spring PRIVATE yjk_core)
  add_test(NAME spring COMMAND test_spring)

# 命令行冒烟：用示例模型跑完整条链路并写出三个结果文件。
  # 这条测试保护的是"各模块能不能接起来"，而不是数值正确性 ——
  # 后者由上面各套件负责。
  add_test(NAME cli_smoke
           COMMAND yjk run ${CMAKE_CURRENT_SOURCE_DIR}/examples/frame3x4.yjk
                   -o ${CMAKE_CURRENT_BINARY_DIR}/cli_smoke)

  # 剪力墙 P2：墙进入轴网/网格/脚本/自重/楼层剪力 整条链路
  #   test_wall —— 节点重合（墙与梁板共用 nodeAt）、悬臂墙解析解（Timoshenko）、
  #                双墙并联（刚度翻倍/反力平分）、墙底剪力对账（层剪力 =
  #                基底水平反力）、墙模态（质量与基频）、端到端 frame_wall.yjk
  add_executable(test_wall tests/test_wall.cpp)
  target_link_libraries(test_wall PRIVATE yjk_core)
  add_test(NAME wall COMMAND test_wall
           ${CMAKE_CURRENT_SOURCE_DIR}/examples/frame_wall.yjk)

  # P3 交互建模 T1：内核可变接口（addColumn/addWall/addNodeLoad/remove/clear
  #   与按 id 索引查询）—— 脚本层与交互层共用同一套接口。
  #   核心验收：Model 建模型 → YjkWriter 序列化 → ModelScript 读回，
  #   节点/单元/荷载集合逐字段一致。
  add_executable(test_interact tests/test_interact.cpp)
  target_link_libraries(test_interact PRIVATE yjk_core)
  add_test(NAME interact COMMAND test_interact)

  # P3-T2 射线拾取核心：yjk::interact 纯函数（无 Qt 依赖），
  # View3D 与脚本层共用同一套命中逻辑，单测覆盖节点/梁/壳/遮挡/tMax。
  add_executable(test_pick tests/test_pick.cpp)
  target_link_libraries(test_pick PRIVATE yjk_core)
  add_test(NAME pick COMMAND test_pick)
endif()

# ============================================================
#  命令行工具 yjk
#
#  这是目前的主要入口：yjk run 模型.yjk 一步走到结果文件。
#  它不依赖 Qt，可以在 CI / 服务器上跑。
# ============================================================
add_executable(yjk apps/yjk.cpp)
target_link_libraries(yjk PRIVATE yjk_core)

# ============================================================
#  Qt 界面
#
#  用 -DYJK_BUILD_GUI=ON 打开。需要 Qt6 的 Core / Gui / Widgets。
#
#  【为什么必须用【Qt 自带的 MinGW】而不是系统里的编译器】
#  Qt 的 mingw_64 包是用特定版本的 MinGW-w64（6.8 对应 GCC 13.1.0，
#  x86_64-posix-seh）构建的。用别的 GCC 编译本程序时，
#  链接阶段会和 Qt 的 DLL 共享 libstdc++-6.dll / libgcc_s_seh-1.dll，
#  版本不一致会出现"能编过、一运行就崩"或异常到处乱飞。
#  正确做法：Qt 在哪装的，就用它 Tools/mingw*_64 里的编译器编本程序。
#  aqt 安装时对应 tools_mingw1310。
# ============================================================
if(YJK_BUILD_GUI)
  find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets)
  set(CMAKE_AUTOMOC ON)

  add_executable(yjk_ui
    apps/ui/main.cpp
    apps/ui/MainWindow.cpp
    apps/ui/MainWindow.h
    apps/ui/View3D.cpp
    apps/ui/View3D.h
    apps/ui/SceneData.cpp
    apps/ui/SceneData.h
    apps/ui/ModelTreePanel.cpp
    apps/ui/ModelTreePanel.h
    apps/ui/PropertyPanel.cpp
    apps/ui/PropertyPanel.h
    apps/ui/ResultPanel.cpp
    apps/ui/ResultPanel.h
    apps/ui/Theme.cpp
    apps/ui/Theme.h
    apps/ui/Camera.h
    apps/ui/ColorMap.h
  )
  target_include_directories(yjk_ui PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/apps)
  target_link_libraries(yjk_ui PRIVATE yjk_core Qt6::Core Qt6::Gui Qt6::Widgets)

  # 把示例脚本复制到可执行文件旁边：
  # 界面里的"开始计算"在没选文件时会去找它，让程序开箱可用。
  add_custom_command(TARGET yjk_ui POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_CURRENT_SOURCE_DIR}/examples
            $<TARGET_FILE_DIR:yjk_ui>/examples)

  # 界面自检：无头渲染一张图。保护的是"渲染管线没被改坏"，
  # 而不是数值正确性（那由内核各测试套件负责）。
  if(YJK_BUILD_TESTS)
    add_test(NAME ui_render
             COMMAND yjk_ui --script ${CMAKE_CURRENT_SOURCE_DIR}/examples/frame3x4.yjk
                            --shot ${CMAKE_CURRENT_BINARY_DIR}/ui_shot.png
                            --preset iso --contour disp --size 1200x760)
    set_tests_properties(ui_render PROPERTIES
      ENVIRONMENT "QT_QPA_PLATFORM=windows"
      TIMEOUT 90)
  endif()
endif()

# ============================================================
#  安装
# ============================================================
install(TARGETS yjk_core ARCHIVE DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)
