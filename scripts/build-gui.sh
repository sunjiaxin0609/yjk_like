#!/usr/bin/env bash
# =============================================================================
#  scripts/build-gui.sh  ——  构建 Qt 桌面端
#
#  【为什么要脚本】配置 GUI 构建的 cmake 命令行有 7 个参数（Qt 路径、
#  MinGW 路径、两个编译器、make 程序、GENERATOR……），少一个就会得到
#  难懂的报错。而且 Qt 与 MinGW 【必须配套】—— 用系统里别的 GCC 编，
#  链接能过但一运行就和 Qt 的 DLL 抢 libstdc++，异常乱飞。
#  这些知识不该靠人记，应该固化在脚本里。
#
#  用法：
#    scripts/build-gui.sh            配置 + 编译
#    scripts/build-gui.sh --deploy   再打一个可独立拷贝的 dist/ 目录
#    scripts/build-gui.sh --run      编译后直接启动
#    scripts/build-gui.sh --shot out.png [--preset iso] [--contour disp]
# =============================================================================
set -euo pipefail

# ---- 路径风格：同一个路径在这里有两种写法，不能混 ----
#
# 【这是个真实踩过的坑】Git Bash 内部按 Unix 风格（/c/...）解析 PATH，
# 而传给 cmake 的 -D 参数必须是 Windows 风格（C:/...）。
# 早期版本把 Windows 风格写进了 PATH，结果 PATH 条目失效 ——
# 症状是 cmake 报 "Test run of moc executable failed / 0xc0000135"，
# 也就是 moc 找不到自己的 DLL。报错信息完全指向 Qt，
# 与真正的原因（PATH 写法）毫无关系，极难定位。
QTDIR_U="${QTDIR:-/c/toolchain/Qt/6.8.3/mingw_64}"
QTMINGW_U="${QTMINGW:-/c/toolchain/Qt/Tools/mingw1310_64}"
# 转成 cmake 需要的 Windows 风格
QTDIR_W="$(cygpath -m "$QTDIR_U" 2>/dev/null || echo "$QTDIR_U")"
QTMINGW_W="$(cygpath -m "$QTMINGW_U" 2>/dev/null || echo "$QTMINGW_U")"
CMAKEROOT="${CMAKEROOT:-/c/toolchain/w64devkit/bin}"
BUILDDIR="${BUILDDIR:-build-gui}"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

CMAKE="$CMAKEROOT/cmake.exe"
CTEST="$CMAKEROOT/ctest.exe"

# ---- 前置检查 ----
# 这三个路径缺一个，cmake 会给出"找不到 Qt6Config.cmake"这类误导性的报错，
# 不如在这里直接说清楚缺什么、去哪儿装。
for p in "$QTDIR_U" "$QTMINGW_U" "$CMAKE"; do
  if [ ! -e "$p" ]; then
    cat >&2 <<EOF
错误：找不到 $p

Qt 6 用 aqtinstall 安装（见 README「Qt 桌面端」一节）：
  pip install aqtinstall
  python -m aqt install-qt windows desktop 6.8.3 win64_mingw \\
      -O C:/toolchain/Qt --archives qtbase
  python -m aqt install-tool windows desktop tools_mingw1310 \\
      qt.tools.win64_mingw1310 -O C:/toolchain/Qt

若装在别处，用环境变量指定：QTDIR=... QTMINGW=... scripts/build-gui.sh
EOF
    exit 1
  fi
done

export PATH="$QTMINGW_U/bin:$QTDIR_U/bin:$PATH"

echo "== 配置 =="
"$CMAKE" -B "$BUILDDIR" -G "MinGW Makefiles" \
  -DCMAKE_BUILD_TYPE=Release \
  -DYJK_BUILD_GUI=ON \
  -DYJK_BUILD_TESTS="${YJK_TESTS:-ON}" \
  -DCMAKE_PREFIX_PATH="$QTDIR_W" \
  -DCMAKE_CXX_COMPILER="$QTMINGW_W/bin/g++.exe" \
  -DCMAKE_MAKE_PROGRAM="$QTMINGW_W/bin/mingw32-make.exe" \
  2>&1 | tail -3

echo "== 编译 =="
"$CMAKE" --build "$BUILDDIR" -j"$(nproc 2>/dev/null || echo 4)"

if [ "${YJK_TESTS:-ON}" = "ON" ]; then
  echo "== 测试 =="
  "$CTEST" --test-dir "$BUILDDIR" --output-on-failure 2>&1 | tail -6
fi

if [ "${1:-}" = "--deploy" ]; then
  echo "== 打包 dist/ =="
  rm -rf dist && mkdir -p dist
  cp "$BUILDDIR/yjk_ui.exe" dist/
  cp -r examples dist/
  windeployqt --release --no-translations --no-system-d3d-compiler \
              --no-opengl-sw --compiler-runtime dist/yjk_ui.exe >/dev/null
  echo "完成：dist/ （$(du -sh dist | cut -f1)，可直接拷到别的机器运行）"
fi

if [ "${1:-}" = "--run" ]; then
  echo "== 启动 =="
  exec "./$BUILDDIR/yjk_ui.exe" ${2:+"$2"}
fi

if [ "${1:-}" = "--shot" ]; then
  shift
  echo "== 截图 =="
  "./$BUILDDIR/yjk_ui.exe" --script examples/frame3x4.yjk --shot "$@" --window
fi

echo "完成。可执行文件：$BUILDDIR/yjk_ui.exe"
