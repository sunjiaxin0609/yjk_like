// =============================================================================
//  yjk/math/Version.h  ——  版本与诊断接口
// =============================================================================
#pragma once

#include "yjk/math/Solver.h"

namespace yjk {

// 版本号（0.1.x = 数学层；0.2.x = 单元层；0.3.x = 静力与模态）
const char* version();

// 编译信息：C++ 标准、编译器、字长
const char* buildInfo();

// 把求解报告转成多行可读诊断（工程师排查问题的入口）
const char* diagnose(const SolveReport& rep);

}  // namespace yjk
