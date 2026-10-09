// =============================================================================
//  apps/ui/Theme.cpp
// =============================================================================
#include "Theme.h"

#include <QFontDatabase>
#include <QLocale>
#include <cmath>

namespace ui {

QFont monoFont(int pt) {
  QFont f = QFontDatabase::systemFont(QFontDatabase::FixedFont);
  f.setPointSize(pt);
  f.setStyleHint(QFont::Monospace);
  return f;
}

QFont uiFont(int pt) {
  QFont f = QFontDatabase::systemFont(QFontDatabase::GeneralFont);
  f.setPointSize(pt);
  return f;
}

// 数值格式化。
//
// 【为什么不用 QString::number(v, 'g', 4)】
// 'g' 在 1000 ~ 999999 这一段会输出科学计数法（1.234e+05），
// 而结构工程师读的是 kN、kN·m 这类工程量，e+05 完全不可读。
// 这里改成：小量走定点、中量走定点、超大量才用科学计数法，
// 并做四舍五入以避免 -0.000 这种负零显示。
QString fmtNum(double v, int sig) {
  if (!std::isfinite(v)) return "—";
  const double a = std::abs(v);
  // 【阈值只能是 0，不能是"很小就当 0"】
  // 曾经用 1e-12 作为"视为零"，结果残差 1.26e-13 显示成 "0" ——
  // 而残差恰恰是"这个解可信"的证据，显示成 0 反而像没算。
  // 结构工程里任何有量纲的量的有效下限都在 1e-9 以上（m / kN / kPa），
  // 比它更小的只可能是无量纲的收敛指标，必须如实显示。
  if (a == 0.0) return "0";
  if (a >= 1e7 || a < 1e-4) {
    return QString::number(v, 'e', 2);
  }
  // 有效数字 → 小数位数
  const int intDigits = static_cast<int>(std::floor(std::log10(a))) + 1;
  int dec = sig - intDigits;
  dec = std::clamp(dec, 0, 6);
  QString s = QString::number(v, 'f', dec);
  // 消除 "-0.000"：负零在读数上是噪音，而且排序时会出现在 0 之前。
  if (v < 0.0 && s.toDouble() == 0.0) s.remove('-');
  return s;
}

}  // namespace ui
