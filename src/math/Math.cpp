// =============================================================================
//  src/math/Math.cpp  ——  数学层的显式实例化与工具函数
//
//  yjk_core 主体是头文件实现（模板为主），此文件仅提供：
//   1. 版本信息
//   2. 非模板的诊断/统计工具
//   3. 显式实例化，减少其他翻译单元的编译开销
// =============================================================================
#include "yjk/math/Types.h"
#include "yjk/math/SparseMatrix.h"
#include "yjk/math/Ordering.h"
#include "yjk/math/Solver.h"

#include <cstdio>
#include <cstring>

namespace yjk {

const char* version() { return "0.1.0"; }
const char* buildInfo() {
  static char buf[256];
  std::snprintf(buf, sizeof(buf),
                "yjk_like %s | C++%ld | %s %d-bit | sizeof(double)=%zu",
                version(), static_cast<long>(__cplusplus),
#if defined(_MSC_VER)
                "MSVC", (int)(sizeof(void*) * 8),
#elif defined(__GNUC__)
                "GCC", (int)(sizeof(void*) * 8),
#else
                "unknown", (int)(sizeof(void*) * 8),
#endif
                sizeof(double));
  return buf;
}

// 求解器综合诊断：把一份求解报告转成工程师能读懂的多行文本。
// 结构软件的报错必须能直接指导排查 —— "奇异"要告诉用户是哪个构件、
// 大概什么原因。这一点比算得快更重要。
const char* diagnose(const SolveReport& rep) {
  static char buf[512];
  std::snprintf(buf, sizeof(buf),
                "自由度 n = %d, 非零元 = %zu, 带宽 = %d, 填充 = %zu (%.1f%%), "
                "最小主元 = %.4e, 耗时 = %.4f s\n诊断：%s",
                rep.n, rep.nnzA, rep.bandwidth, rep.fill,
                rep.nnzA ? 100.0 * static_cast<double>(rep.fill) / static_cast<double>(rep.nnzA) : 0.0,
                rep.minPivot, rep.seconds, rep.message().c_str());
  return buf;
}

}  // namespace yjk
